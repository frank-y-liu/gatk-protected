package org.broadinstitute.hellbender.tools.exome.eval;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.EnumSet;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Filters applied to indicidual truth or called segments when taking them into account for evaluation.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public enum EvaluationSegmentFilter {
    /**
     * The segment has a low quality value.
     */
    LowQuality("LQ"),

    /**
     * The segment is to short to consider.
     */
    ShortEvent("SE");

    /**
     * Represent an empty filter set, i.e. the segment pass all filters.
     */
    public static final String PASS_ACRONYM = "PASS";

    public final String acronym;

    EvaluationSegmentFilter(final String acronym) {
        this.acronym = Utils.nonNull(acronym);
    }

    public static final String toString(final Set<EvaluationSegmentFilter> filters) {
        Utils.nonNull(filters);
        if (filters.isEmpty()) {
            return PASS_ACRONYM;
        } else {
            return filters.stream()
                    .map(f -> {
                        if (f == null) {
                            throw new IllegalArgumentException("the input filters set cannot contain a null");
                        } else {
                            return f.acronym;
                        }})
                    .collect(Collectors.joining(","));
        }
    }

    /**
     * Transform a single filter string representation into the corresponding filter enum element.
     * @param string the string
     * @return never code
     */
    public static EvaluationSegmentFilter toFilter(final String string) {
        final String trimmedString = Utils.nonNull(string, "the input string cannot be null").trim();
        if (trimmedString.equalsIgnoreCase(LowQuality.acronym)) {
            return LowQuality;
        } else if (trimmedString.equalsIgnoreCase(ShortEvent.acronym)) {
            return ShortEvent;
        } else {
            return EvaluationSegmentFilter.valueOf(trimmedString);
        }
    }

    /**
     * Returns an immutable set containing all the filters encoded in the input string.
     * @param string the input string.
     * @return never {@code null}, but perhaps the empty set.
     */
    public static Set<EvaluationSegmentFilter> toSet(final String string) {
        final String trimmedString = Utils.nonNull(string, "the input string cannot be null").trim();
        if (trimmedString.equalsIgnoreCase(PASS_ACRONYM)) {
            return EnumSet.noneOf(EvaluationSegmentFilter.class);
        } else {
            return EnumSet.copyOf(Stream.of(string.split(","))
                    .map(EvaluationSegmentFilter::toFilter)
                    .collect(Collectors.toSet()));
        }

    }

}
