package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.exome.CopyNumberTriStateFrequencies;
import org.broadinstitute.hellbender.tools.exome.CopyNumberTriStateSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Segment Evaluation detail output record.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class EvaluationSiteRecord implements Locatable {

    public final String sample;
    public final List<Pair<CopyNumberTriStateSegment, Set<EvaluationSegmentFilter>>> calls;
    public final List<Pair<CopyNumberTriStateSegment, Set<EvaluationSegmentFilter>>> truths;
    public final EvaluationClass evaluationClass;
    public final SimpleInterval interval;
    public final int targetCount;
    public final double deletionFrequency;
    public final double duplicationFrequency;

    private EvaluationSiteRecord(final String sample, final int targetCount,
                                 final EvaluationClass evaluationClass,
                                 final List<Pair<CopyNumberTriStateSegment, Set<EvaluationSegmentFilter>>> calls,
                                 final List<Pair<CopyNumberTriStateSegment, Set<EvaluationSegmentFilter>>> truths,
                                 final CopyNumberTriStateFrequencies frequencies) {
        this.sample = Utils.nonNull(sample);
        this.evaluationClass = Utils.nonNull(evaluationClass);
        this.calls = Collections.unmodifiableList(new ArrayList<>(Utils.nonNull(calls)));
        this.truths = Collections.unmodifiableList(new ArrayList<>(Utils.nonNull(truths)));
        if (truths.size() == 1) {
            interval = truths.get(0).getLeft().getInterval();
        } else if (calls.size() == 1) {
            interval = calls.get(0).getLeft().getInterval();
        } else {
            throw new IllegalArgumentException("either the truth or the call list must have exactly one element");
        }
        this.targetCount = targetCount;
        this.deletionFrequency = frequencies == null ? Double.NaN : frequencies.frequencyFor(CopyNumberTriState.DELETION, 1.0);
        this.duplicationFrequency = frequencies == null ? Double.NaN : frequencies.frequencyFor(CopyNumberTriState.DUPLICATION, 1.0);
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    private static EvaluationSiteRecord falseNegative(final String sample, final int targetCount, final Pair<CopyNumberTriStateSegment, Set<EvaluationSegmentFilter>> truth, final CopyNumberTriStateFrequencies frequencies) {
        return new EvaluationSiteRecord(sample, targetCount, EvaluationClass.FALSE_NEGATIVE,
                Collections.emptyList(), Collections.singletonList(Utils.nonNull(truth)), frequencies);
    }

    private static EvaluationSiteRecord unknownPositive(final String sample, final int targetCount, final Pair<CopyNumberTriStateSegment, Set<EvaluationSegmentFilter>> call, final CopyNumberTriStateFrequencies frequencies) {
        return new EvaluationSiteRecord(sample, targetCount, EvaluationClass.UNKNOWN_POSITIVE,
                Collections.singletonList(Utils.nonNull(call)), Collections.emptyList(), frequencies);
    }

    private static EvaluationSiteRecord otherPositive(final String sample, final int targetCount, final EvaluationClass positiveClass,
                                                      final Pair<CopyNumberTriStateSegment, Set<EvaluationSegmentFilter>> truth, final List<Pair<CopyNumberTriStateSegment, Set<EvaluationSegmentFilter>>> calls, final CopyNumberTriStateFrequencies frequencies) {
        return new EvaluationSiteRecord(sample, targetCount, positiveClass, Utils.nonNull(calls), Collections.singletonList(Utils.nonNull(truth)), frequencies);
    }

    /**
     * Returns the truth quality as the maximum quality amongst all the truth segments.
     * @return {@link Double#NaN} if there is no truth segment.
     */
    public double truthQuality() {
        return truths.stream()
                .map(Pair::getLeft)
                .mapToDouble(CopyNumberTriStateSegment::getSomeQuality).max().orElse(Double.NaN);
    }

    /**
     * Returns the called quality as the maximum quality amongst all called segments that
     * match the truth state call. If there is no truth segment, then we take the maximum quality amongst all
     * called segments listed.
     * @return {@link Double#NaN} if there is no called segments.
     */
    public double callQuality() {
        final Set<CopyNumberTriState> truthStates = truthCallStates();
        return calls.stream()
                .filter(s -> truthStates.contains(s.getLeft().getCall()))
                .map(Pair::getLeft)
                .mapToDouble(CopyNumberTriStateSegment::getSomeQuality).max().orElse(Double.NaN);
    }

    private Set<CopyNumberTriState> truthCallStates() {
        return truths.stream()
                    .map(Pair::getLeft)
                    .map(CopyNumberTriStateSegment::getCall)
                    .collect(Collectors.toSet());
    }

    /**
     * Returns the length in terms of target for the truth segment in this record.
     * <p>
     *     This is calculated as the maximum target count across all truth calls.
     * </p>
     * <p>
     *     When there is no truth, the length is 0.
     * </p>
     * @return 0 or greater.
     */
    public int truthTargetCount() {
        return (int) truths.stream()
                .map(Pair::getLeft)
                .mapToLong(CopyNumberTriStateSegment::getTargetCount).max().orElse(0);
    }

    /**
     * Returns the length in terms of target for the called segment in this record.
     * <p>
     *     This is calculated as the target count that best matches the truth.
     * <p>
     *     When there is no truth, the length is 0.
     * </p>
     * @return 0 or greater.
     */
    public int callTargetCount() {
        final Set<CopyNumberTriState> truthStates = truthCallStates();
        return (int) calls.stream()
                .map(Pair::getLeft)
                .filter(s -> truthStates.contains(s.getCall()))
                .mapToLong(CopyNumberTriStateSegment::getTargetCount).max().orElse(0);
    }

    @Override
    public String toString() {
        return interval.toString();
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }
}
