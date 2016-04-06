package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.exome.CopyNumberTriStateFrequencies;
import org.broadinstitute.hellbender.tools.exome.CopyNumberTriStateSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Segment Evaluation detail output record.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class EvaluationSiteRecord implements Locatable {

    public final String sample;
    public final EvaluationClass evaluationClass;
    public final VariantContext truth;
    public final Set<EvaluationSegmentFilter> filters;
    public final List<VariantContext> calls;
    public final SimpleInterval interval;
    public final int targetCount;
    public final VariantContext result;

    public EvaluationSiteRecord(final String sample,
                                final SimpleInterval interval,
                                final int targetCount,
                                final EvaluationClass evaluationClass,
                                final Set<EvaluationSegmentFilter> filters,
                                final VariantEvaluationContext result) {
        this.sample = Utils.nonNull(sample);
        this.result = new VariantContextBuilder(Utils.nonNull(result)).genotypes(Collections.singleton(result.getGenotype(sample))).make();
        this.evaluationClass = Utils.nonNull(evaluationClass);
        this.truth = result.getTruthVariantContext() == null ? null : new VariantContextBuilder(result.getTruthVariantContext()).genotypes(Collections.singleton(result.getTruthVariantContext().getGenotype(sample))).make();
        this.calls = result.getCallsVariantContexts() == null ? Collections.emptyList() : result.getCallsVariantContexts().stream().map(var -> new VariantContextBuilder(var).genotypes(var.getGenotypes(sample)).make()).collect(Collectors.toList());
        this.filters = Collections.unmodifiableSet(new LinkedHashSet<>(filters));
        this.interval = Utils.nonNull(interval);
        this.targetCount = targetCount;
        if (calls.contains(null)) {
            throw new IllegalArgumentException("calls cannot contain nulls");
        }
    }

    public String getFilterString() {
        if (filters.isEmpty()) {
            return VCFConstants.PASSES_FILTERS_v4;
        } else {
            return filters.stream().sorted().map(Enum<EvaluationSegmentFilter>::name).collect(Collectors.joining(VCFConstants.FILTER_CODE_SEPARATOR));
        }
    }

    public SimpleInterval getInterval() {
        return interval;
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
