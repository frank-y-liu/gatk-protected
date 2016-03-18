package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.utils.mcmc.AbstractParameterizedState;
import org.broadinstitute.hellbender.utils.mcmc.Parameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The state of the allele fraction model without hyperparameters for the allelic bias, containing: <p>
 *      1.  minor allele fractions for each segment <p>
 *      2.  a global outlier probability <p>
 * <p>
 * See docs/CNVs/CNV-methods.pdf for details.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionStateWithoutHyperparameters extends AbstractParameterizedState {
    @SuppressWarnings("serial")
    public static final class MinorFractions extends ArrayList<Double> {
        public MinorFractions() { super(); }
        public MinorFractions(final List<Double> other) {
            super(new ArrayList<>(other));
        }
    }

    public static final String P_OUTLIER_NAME = "pi";
    public static final String MINOR_FRACTIONS_NAME = "f";

    @Override
    protected <S extends AbstractParameterizedState> S copy(final Class<S> stateClass) {
        return stateClass.cast(new AlleleFractionStateWithoutHyperparameters(outlierProbability(), new MinorFractions(minorFractions())));
    }

    //copy the state with reference to the SAME minorFractions list data (to save copying) but different value of
    //one of the scalar parameters
    //This is dangerous and minorFractions should not be modified in the copy.
    //
    //The purpose of this is to make an MCMC proposal state for calculating a likelihood with one of the scalar parameters
    //modified (these are unboxed in the getters, so changing these in the copy is safe)
    public AlleleFractionStateWithoutHyperparameters shallowCopyWithProposedOutlierProbability(final double proposedOutlierProbability) {
        return new AlleleFractionStateWithoutHyperparameters(proposedOutlierProbability, minorFractions());
    }

    public AlleleFractionStateWithoutHyperparameters(final double outlierProbability,
                                                     final MinorFractions minorFractions) {
        super(Arrays.asList(
                new Parameter<>(P_OUTLIER_NAME, outlierProbability),
                new Parameter<>(MINOR_FRACTIONS_NAME, minorFractions)));
    }

    public double outlierProbability() {
        return get(P_OUTLIER_NAME, Double.class);
    }

    public MinorFractions minorFractions() {
        return get(MINOR_FRACTIONS_NAME, MinorFractions.class);
    }

    public double minorFractionInSegment(final int segment) {
        return get(MINOR_FRACTIONS_NAME, MinorFractions.class).get(segment);
    }
}
