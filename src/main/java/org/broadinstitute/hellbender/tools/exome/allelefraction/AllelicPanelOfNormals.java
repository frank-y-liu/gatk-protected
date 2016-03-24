package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import static java.lang.Math.sqrt;

/**
 * Represents the panel of normals used for allele-bias correction.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPanelOfNormals {
    private static final Logger logger = LogManager.getLogger(AllelicPanelOfNormals.class);

    public static final AllelicPanelOfNormals EMPTY_PON = new AllelicPanelOfNormals();

    private final Map<SimpleInterval, HyperparameterValues> siteToHyperparameterPairMap = new HashMap<>();
    private HyperparameterValues mleHyperparameterValues;

    public AllelicPanelOfNormals() {}

    public AllelicPanelOfNormals(final File inputFile) {
        Utils.nonNull(inputFile);
        Utils.regularReadableUserFile(inputFile);

        final AllelicCountCollection counts = new AllelicCountCollection(inputFile);
        mleHyperparameterValues = calculateMLEHyperparameterValues(counts);
        initializeSiteToHyperparameterPairMap(counts);
    }

    public double getAlpha(final SimpleInterval site) {
        Utils.nonNull(site);
        return siteToHyperparameterPairMap.getOrDefault(site, mleHyperparameterValues).alpha;
    }

    public double getBeta(final SimpleInterval site) {
        Utils.nonNull(site);
        return siteToHyperparameterPairMap.getOrDefault(site, mleHyperparameterValues).beta;
    }

    private class HyperparameterValues {
        private final double alpha;
        private final double beta;

        private HyperparameterValues(final double alpha, final double beta) {
            this.alpha = alpha;
            this.beta = beta;
        }

        private HyperparameterValues(final int a, final int r) {
            final double f = 0.5;
            final int n = a + r;
            final double alpha = mleHyperparameterValues.alpha;
            final double beta = mleHyperparameterValues.beta;
            final double w = (1 - f) * (a - alpha + 1) + beta * f;
            final double lambda0 = (sqrt(w * w + 4 * beta * f * (1 - f) * (r + alpha - 1)) - w) / (2 * beta * (1 - f));
            final double y = (1 - f)/(f + (1 - f) * lambda0);
            final double kappa = n * y * y - (r + alpha - 1) / (lambda0 * lambda0);
            this.alpha = 1 - kappa * lambda0 * lambda0;
            this.beta = -kappa * lambda0;
        }
    }

    /**
     * Find MLE hyperparameter values from counts in panel of normals.  See analogous code in {@link AlleleFractionInitializer}.
     */
    private HyperparameterValues calculateMLEHyperparameterValues(final AllelicCountCollection counts) {
        double meanBias = AlleleFractionInitializer.INITIAL_MEAN_BIAS;
        double biasVariance = AlleleFractionInitializer.INITIAL_BIAS_VARIANCE;
        double previousIterationLogLikelihood;
        double nextIterationLogLikelihood = Double.NEGATIVE_INFINITY;
        int iteration = 1;
        logger.info(String.format("Initializing MLE hyperparameter values for allelic panel of normals.  Iterating until log likelihood converges to within %.3f.",
                AlleleFractionInitializer.LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD));
        do {
            previousIterationLogLikelihood = nextIterationLogLikelihood;
            meanBias = estimateMeanBias(meanBias, biasVariance, counts);
            biasVariance = estimateBiasVariance(meanBias, biasVariance, counts);
            nextIterationLogLikelihood = logLikelihood(meanBias, biasVariance, counts);
            logger.info(String.format("Iteration %d, model log likelihood = %.3f.", iteration, nextIterationLogLikelihood));
            iteration++;
        } while (iteration < AlleleFractionInitializer.MAX_ITERATIONS &&
                nextIterationLogLikelihood - previousIterationLogLikelihood > AlleleFractionInitializer.LOG_LIKELIHOOD_CONVERGENCE_THRESHOLD);

        final double alpha = alpha(meanBias, biasVariance);
        final double beta = beta(meanBias, biasVariance);
        logger.info("MLE hyperparameter values for allelic panel of normals found:");
        logger.info("alpha = " + alpha);
        logger.info("beta = " + beta);
        return new HyperparameterValues(alpha, beta);
    }

    private double estimateMeanBias(final double meanBias, final double biasVariance, final AllelicCountCollection counts) {
        final UnivariateObjectiveFunction objective = new UnivariateObjectiveFunction(proposedMeanBias ->
                logLikelihood(proposedMeanBias, biasVariance, counts));
        final SearchInterval searchInterval = new SearchInterval(0.0, AlleleFractionInitializer.MAX_REASONABLE_MEAN_BIAS, meanBias);
        return AlleleFractionInitializer.OPTIMIZER.optimize(objective, GoalType.MAXIMIZE, searchInterval, AlleleFractionInitializer.BRENT_MAX_EVAL).getPoint();
    }

    private double estimateBiasVariance(final double meanBias, final double biasVariance, final AllelicCountCollection counts) {
        final UnivariateObjectiveFunction objective = new UnivariateObjectiveFunction(proposedBiasVariance ->
                logLikelihood(meanBias, proposedBiasVariance, counts));
        final SearchInterval searchInterval = new SearchInterval(0.0, AlleleFractionInitializer.MAX_REASONABLE_BIAS_VARIANCE, biasVariance);
        return AlleleFractionInitializer.OPTIMIZER.optimize(objective, GoalType.MAXIMIZE, searchInterval, AlleleFractionInitializer.BRENT_MAX_EVAL).getPoint();
    }

    private double logLikelihood(final double meanBias, final double biasVariance, final AllelicCountCollection counts) {
        final double alpha = alpha(meanBias, biasVariance);
        final double beta = beta(meanBias, biasVariance);
        return counts.getCounts().stream().mapToDouble(c -> AlleleFractionLikelihoods.logPhi(alpha, beta, 0.5, c.getAltReadCount(), c.getRefReadCount())).sum();
    }

    private void initializeSiteToHyperparameterPairMap(final AllelicCountCollection counts) {
        logger.info("Initializing allelic panel of normals...");
        for (final AllelicCount count : counts.getCounts()) {
            final SimpleInterval site = count.getInterval();
            final HyperparameterValues hyperparameterValues = new HyperparameterValues(count.getAltReadCount(), count.getRefReadCount());
            if (siteToHyperparameterPairMap.containsKey(site)) {
                throw new UserException.BadInput("Input file for allelic panel of normals contains duplicate sites.");
            } else {
                siteToHyperparameterPairMap.put(site, hyperparameterValues);
            }
        }
        logger.info("Allelic panel of normals initialized.");
    }

    private double alpha(final double meanBias, final double biasVariance) {
        return meanBias * meanBias / biasVariance;
    }

    private double beta(final double meanBias, final double biasVariance) {
        return meanBias / biasVariance;
    }
}
