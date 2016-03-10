package org.broadinstitute.hellbender.tools.exome.eval;

/**
 * Segment evaluation classes: TP, FN, etc.
 */
public enum EvaluationClass {

    /**
     * When truth and calls overlap and are all compatible (all are deletion or all are duplication).
     */
    TRUE_POSITIVE("TP"),

    /**
     * When truth does not overlap with any call.
     */
    FALSE_NEGATIVE("FN"),

    /**
     * When truth overlaps with several calls that are discordant amongst them.
     */
    MIXED_POSITIVE("MP"),

    /**
     * When truth overlaps with several calls that are concordant amongst them but
     * the are discordant with the truth.
     */
    DISCORDANT_POSITIVE("DP"),

    /**
     * When a call does not overall with truth.
     */
    UNKNOWN_POSITIVE("UP");

    /**
     * Short acronym name for the class.
     */
    public final String acronym;

    /**
     * Creates a new evaluation class type given its acronym.
     * @param acronym
     */
    EvaluationClass(final String acronym) {
        this.acronym = acronym;
    }

    @Override
    public String toString() {
        return acronym;
    }
}
