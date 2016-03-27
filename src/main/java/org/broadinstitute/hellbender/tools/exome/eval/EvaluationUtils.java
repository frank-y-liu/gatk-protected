package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Comparator;
import java.util.Objects;

/**
 * Evaluation Utils.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class EvaluationUtils {

    /**
     * Declared to prevent instantiation of this helping class.
     */
    private EvaluationUtils() {}

    /**
     * Returns a comparator that sort locatable elements from the ones with better overlap with a query locatable to
     * the ones with a worse overlap.
     *
     * <p> The goodness of an overlap is quantified as the "ratio" between the overlap base pairs count divided by
     * the total union of the element and the query locatable. When the query and the target element actually overlap
     * this ratio is positive with a maximum value of 1.0 when both locatable have exactly the same coordinates.
     * <p>
     * In contrast, when both regions do not overlap, this ratios becomes negative and it is calculate as the distance
     * in base-pairs divided by the total size of the region that just enclosed both locatables. This can take values from
     * 0 to just under -1.0 for locatables located on the same contig.
     * </p>
     * When both locatables are based on different contigs then this ratio is effectively -Inf so sorted towards the end.
     * <p>
     * Null elements are accepted but systematically sorted at the very end.
     * </p>
     *
     * the ratio is negative when both regions don't overlap actually do
     * not overlap with larger negative values
     * </p>
     * @param query they query locatable.
     * @param <E> the locatable element type to compare.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code query == null}.
     */
    public static <E extends Locatable> Comparator<E> comparingOverlapToUnionRatio(
            final Locatable query) {
        Utils.nonNull(query);
        return new Comparator<E>() {
            @Override
            public int compare(final E a, final E b) {
                if (a == null) {
                    return b == null ? 0 : 1;
                } else if (b == null) {
                    return -1; // a != null at this point.
                } else if (Objects.equals(a.getContig(), b.getContig())) {
                    if (Objects.equals(a.getContig(), query.getContig())) {
                        final double aRatio = calculateOverlapToUnionRatio(a);
                        final double bRatio = calculateOverlapToUnionRatio(b);
                        return Double.compare(bRatio, aRatio);
                    } else {
                        return 0;
                    }
                } else if (Objects.equals(a.getContig(), query.getContig())) {
                    return -1;
                } else if (Objects.equals(b.getContig(), query.getContig())) {
                    return 1;
                } else {
                    // when neither a or b contig match the query then we consider both equal.
                    return 0;
                }
            }

            private double calculateOverlapToUnionRatio(final E target) {
                return (Math.min(target.getEnd(), query.getEnd()) - Math.max(target.getStart(), query.getStart()) + 1) /
                        (Math.max(target.getEnd(), query.getEnd()) - Math.max(target.getStart(), query.getStart()) + 1);
            }
        };
    }


}
