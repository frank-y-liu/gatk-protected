package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Segment evaluation classes: TP, FN, etc.
 */
public enum EvaluationClass {

    /**
     * When truth and calls overlap and are all compatible (all are deletion or all are duplication).
     */
    TRUE_POSITIVE("TP", "When truth and calls overlap and are all compatible (all are deletion or all are duplication)"),

    /**
     * When truth does not overlap with any call.
     */
    FALSE_NEGATIVE("FN", "When truth does not overlap with any call"),

    /**
     * When truth overlaps with several calls that are discordant amongst them.
     */
    MIXED_POSITIVE("MP", "When truth overlaps with several calls that are discordant amongst them"),

    /**
     * When truth overlaps with several calls that are concordant amongst them but
     * the are discordant with the truth.
     */
    DISCORDANT_POSITIVE("DP", "When truth overlaps with several calls that are concordant amongst them but the are discordant with the truth"),

    /**
     * When a call does not overall with truth.
     */
    UNKNOWN_POSITIVE("UP", "When a call does not overall with truth");

    /**
     * The key for vcf header lines describing evaluation classes.
     */
    public static final String VCF_HEADER_KEY = "EVAL_CLASS";

    private static final Map<String, EvaluationClass> VALUE_BY_ACRONYM =
            Collections.unmodifiableMap(Stream.of(values())
                    .collect(Collectors.toMap(s -> s.acronym, s -> s)));

    /**
     * Short acronym name for the class.
     */
    public final String acronym;

    /**
     * Description of this class meaning.
     */
    public final String description;

    /**
     * Creates a new evaluation class type given its acronym.
     * @param acronym
     */
    EvaluationClass(final String acronym, final String description) {
        this.acronym = Utils.nonNull(acronym);
        this.description = Utils.nonNull(description);
    }

    @Override
    public String toString() {
        return acronym;
    }

    public VCFHeaderLine headerLine() {
        return new VCFSimpleHeaderLine(VCF_HEADER_KEY, acronym, description);
    }

    /**
     * Add all pertinent meta-data header lines to a vcf header.
     * @param header the target header.
     * @throws IllegalArgumentException if {@code header} is {@code null}.
     */
    public static void addHeaderLinesTo(final VCFHeader header) {
        for (final EvaluationClass value : values()) {
            header.addMetaDataLine(value.headerLine());
        }
    }

    public static EvaluationClass fromAcronym(final String acronym) {
        return VALUE_BY_ACRONYM.get(acronym);
    }
}
