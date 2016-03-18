package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

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
    private final Map<SimpleInterval, HyperparameterValues> siteToHyperparameterPairMap = new HashMap<>();

    private final HyperparameterValues mleHyperparameterValues;

    public AllelicPanelOfNormals(final File inputFile) {
        Utils.nonNull(inputFile);
        Utils.regularReadableUserFile(inputFile);

        final AllelicCountCollection counts = new AllelicCountCollection(inputFile);
        mleHyperparameterValues = calculateMLEHyperparameterValues(counts);

        for (final AllelicCount count : counts.getCounts()) {
            final SimpleInterval site = count.getInterval();
            final HyperparameterValues hyperparameterValues = new HyperparameterValues(count.getRefReadCount(), count.getAltReadCount());
            if (siteToHyperparameterPairMap.containsKey(site)) {
                throw new UserException.BadInput("Input file for allelic panel of normals contains duplicate sites.");
            } else {
                siteToHyperparameterPairMap.put(site, hyperparameterValues);
            }
        }
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

        private HyperparameterValues(final int r, final int a) {
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

    private HyperparameterValues calculateMLEHyperparameterValues(final AllelicCountCollection counts) {
        //TODO MLE alpha beta

        final double alpha = 0.;
        final double beta = 0.;

        return new HyperparameterValues(alpha, beta);
    }
}
