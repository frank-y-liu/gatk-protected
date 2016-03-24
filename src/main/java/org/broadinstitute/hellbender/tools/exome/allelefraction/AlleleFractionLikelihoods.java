package org.broadinstitute.hellbender.tools.exome.allelefraction;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.special.Gamma;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collection;
import java.util.function.Function;
import java.util.stream.IntStream;

import static java.lang.Math.log;
import static java.lang.Math.sqrt;
import static org.broadinstitute.hellbender.utils.MathUtils.log10Factorial;
import static org.broadinstitute.hellbender.utils.MathUtils.log10ToLog;

/**
 * Contains likelihood methods for the allele-fraction model.
 * See docs/CNVs/CNV-methods.pdf for a thorough description of the model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionLikelihoods {
    private AlleleFractionLikelihoods() {}

    /**
     * Compute the log-likelihood of a alt reads and r ref reads given minor fraction f and gamma hyperparameters
     * (specifying the distribution on allelic biases) mu (mean) and beta (rate = mean/variance) and given
     * an alt minor, ref minor, or outlier indicator state.  Note that this is a partially collapsed log-likelihood in that the
     * latent variable corresponding to the allelic bias at this site has been marginalized out but the indicator
     * variable has not been marginalized out.
     * <p>
     * See docs/CNVs/CNV-methods.pdf for derivation.
     * <p>
     * Finally, note that this is a static method and does not get mu, beta, and minorFraction from an AlleleFractionState object
     * We need such functionality because MCMC evaluates the likelihood under proposed parameter changes.
     *
     * @param state allele fraction state
     * @param segment index of segment containijng this het site
     * @param count AllelicCount of alt and ref reads
     * @param indicator the hidden state (alt minor / ref minor / outlier)
     *
     * @return if indicator == ALT_MINOR:
     * <p>
     * log { [beta^alpha / Gamma(alpha)][(1-pi)/2] * int_{0 to infty} f^a * (1-f)^r * lambda^(alpha + r - 1) * exp(-beta*lambda)/(f + (1-f)*lambda)^n d lambda }
     * <p>
     * if indicator == REF_MINOR, same as ALT_MINOR but with f <--> 1 - f
     * <p>
     * if indicator == OUTLIER log {pi * a!r!/(n+1)!}
     * <p>
     * where alpha = mu*beta and n = a + r
     */
    public static double hetLogLikelihood(final AlleleFractionState state, final int segment, final AllelicCount count, final AlleleFractionIndicator indicator) {
        return hetLogLikelihood(state, segment, count, indicator, AllelicPanelOfNormals.EMPTY_PON);
    }

    public static double hetLogLikelihood(final AlleleFractionState state, final int segment, final AllelicCount count, final AlleleFractionIndicator indicator,
                                          final AllelicPanelOfNormals allelicPON) {
        final SimpleInterval site = count.getInterval();
        final double alpha;
        final double beta;
        if (allelicPON.equals(AllelicPanelOfNormals.EMPTY_PON)) {
            beta = state.meanBias() / state.biasVariance();
            alpha = state.meanBias() * beta;

        } else {
            alpha = allelicPON.getAlpha(site);
            beta = allelicPON.getBeta(site);
        }
        final double pi = state.outlierProbability();
        final double minorFraction = state.minorFractionInSegment(segment);
        final int a = count.getAltReadCount();
        final int r = count.getRefReadCount();

        if (indicator == AlleleFractionIndicator.OUTLIER) {
            return log(pi) + log10ToLog(log10Factorial(a) + log10Factorial(r) - log10Factorial(a + r + 1));
        } else {
            final double f = indicator == AlleleFractionIndicator.ALT_MINOR ? minorFraction : 1 - minorFraction;
            return log((1 - pi) / 2) + logPhi(alpha, beta, f, a, r);
        }
    }

    public static double logPhi(final double alpha, final double beta, final double f, final int a, final int r) {
        final double lambda0 = lambda0(alpha, beta, f, a, r);
        final int n = a + r;
        final double kappa = kappa(alpha, f, r, n, lambda0);
        final double rho = rho(lambda0, kappa);
        final double tau = tau(lambda0, kappa);
        final double logc = alpha*log(beta) - Gamma.logGamma(alpha) + a * log(f) + r * log(1 - f)
                + (r + alpha - rho) * log(lambda0) + (tau - beta) * lambda0 - n * log(f + (1 - f) * lambda0);
        return logc + Gamma.logGamma(rho) - rho * log(tau);
    }

    public static double lambda0(final double alpha, final double beta, final double f, final int a, final int r) {
        final double w = (1 - f) * (a - alpha + 1) + beta * f;
        return (sqrt(w * w + 4 * beta * f * (1 - f) * (r + alpha - 1)) - w) / (2 * beta * (1 - f));
    }

    public static double kappa(final double alpha, final double f, final int r, final int n, final double lambda0) {
        final double y = (1 - f)/(f + (1 - f) * lambda0);
        return n * y * y - (r + alpha - 1) / (lambda0 * lambda0);
    }

    public static double rho(final double lambda0, final double kappa) {
        return 1 - kappa * lambda0 * lambda0;
    }

    public static double tau(final double lambda0, final double kappa) {
        return -kappa * lambda0;
    }

    /**
     * the log likelihood summed (marginalized) over indicator states, which we use in the fully collapsed model
     * in which latent variables (bias and indicator) are marginalized out
     *
     * @param state allele fraction state
     * @param segment index of segment containing this het site
     * @param count AllelicCount of alt and ref reads
     * @return the log of the likelihood at this het site, marginalized over indicator states.
     */
    public static double collapsedHetLogLikelihood(final AlleleFractionState state, final int segment, final AllelicCount count, final AllelicPanelOfNormals allelicPON) {
        return logSumLog(hetLogLikelihood(state, segment, count, AlleleFractionIndicator.ALT_MINOR, allelicPON),
                hetLogLikelihood(state, segment, count, AlleleFractionIndicator.REF_MINOR, allelicPON),
                hetLogLikelihood(state, segment, count, AlleleFractionIndicator.OUTLIER, allelicPON));
    }

    /**
     * the log-likelihood of all the hets in a segment
     *
     * @param state allele fraction state
     * @param segment index of segment containijng this het site
     * @param counts AllelicCount of alt and ref reads in this segment
     * @return the sum of log-likelihoods over all het sites in a segment
     */
    public static double segmentLogLikelihood(final AlleleFractionState state, final int segment, final Collection<AllelicCount> counts, final AllelicPanelOfNormals allelicPON) {
        return counts.stream().mapToDouble(c -> collapsedHetLogLikelihood(state, segment, c, allelicPON)).sum();
    }

    /**
     * the total log likelihood of all segments
     * @param state current state
     * @param data data
     * @return sum of log likelihoods of all segments
     */
    public static double logLikelihood(final AlleleFractionState state, final AlleleFractionData data) {
        return IntStream.range(0, data.numSegments()).mapToDouble(s -> segmentLogLikelihood(state, s, data.countsInSegment(s), data.getPON())).sum();
    }

    protected static Function<Double, Double> logConditionalOnMinorFraction(final AlleleFractionState state,
                                                                            final AlleleFractionData data, final int segment) {
        return logConditionalOnMinorFraction(state, data, segment, AllelicPanelOfNormals.EMPTY_PON);
    }

    protected static Function<Double, Double> logConditionalOnMinorFraction(final AlleleFractionState state,
                                                                            final AlleleFractionData data, final int segment, final AllelicPanelOfNormals allelicPON) {
        return minorFraction -> {
            final AlleleFractionState proposal = new AlleleFractionState(state.meanBias(), state.biasVariance(), state.outlierProbability(), minorFraction);
            return AlleleFractionLikelihoods.segmentLogLikelihood(proposal, 0, data.countsInSegment(segment), allelicPON);
        };
    }

    //compute log(e^a + e^b + e^c) = log(e^M [e^(a-M) + e^(b-M) + e^(c-M)]) where M = max(a,b,c)
    private static double logSumLog(final double a, final double b, final double c) {
        final double maxValue = Doubles.max(a, b, c);
        return maxValue + Math.log(Math.exp(a-maxValue) + Math.exp(b-maxValue) + Math.exp(c-maxValue));
    }
}
