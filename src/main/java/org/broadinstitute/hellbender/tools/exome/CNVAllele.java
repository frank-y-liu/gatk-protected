package org.broadinstitute.hellbender.tools.exome;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Enum of possible CNV alleles.
 */
public enum CNVAllele {
    REF(CopyNumberTriState.NEUTRAL, true),
    DEL(CopyNumberTriState.DELETION, false),
    DUP(CopyNumberTriState.DUPLICATION, false);

    public static final String ALT_KEY = "ALT";
    public static final String DEL_VCF_DESCRIPTION = "Represents a deletion with respect to reference copy number";
    public static final String DUP_VCF_DESCRIPTION = "Represents a duplication with respect to reference copy number";

    public final Allele allele;
    public final CopyNumberTriState state;

    public static final List<CNVAllele> ALL_CNV_ALLELES = Collections.unmodifiableList(Arrays.asList(values()));
    public static final List<CNVAllele> ALTERNATIVE_CNV_ALLELES = ALL_CNV_ALLELES.subList(1, 3);
    public static final List<Allele> ALL_ALLELES = Collections.unmodifiableList(ALL_CNV_ALLELES.stream()
                    .map(a -> a.allele).collect(Collectors.toList()));

    CNVAllele(final CopyNumberTriState state, final boolean isRef) {
        this.state = state;
        allele = isRef ? Allele.create("N", true) : Allele.create("<" + name() + ">");
    }

    /**
     * Adds the alternative allele ALT meta-data lines to a vcf-header.
     * @param header the header to add the lines to.
     * @throws IllegalArgumentException if {@code header} is {@code null}.
     */
    public static void addHeaderLinesTo(final VCFHeader header) {
        Utils.nonNull(header);
        header.addMetaDataLine(new VCFSimpleHeaderLine(ALT_KEY, DEL.allele.getBaseString(), Utils.nonNull(DEL_VCF_DESCRIPTION)));
        header.addMetaDataLine(new VCFSimpleHeaderLine(ALT_KEY, DUP.allele.getBaseString(), Utils.nonNull(DUP_VCF_DESCRIPTION)));
    }

    /**
     * Returns the value whose allele is the same as the one provided.
     * @param allele the query allele.
     * @return never {@code null}.
     * @throws IllegalArgumentException if the input allele is {@code null} or there is no
     * value with such an allele.
     */
    public static CNVAllele valueOf(final Allele allele) {
        Utils.nonNull(allele, "the input allele cannot be null");
        return Stream.of(values()).filter(v -> v.allele.equals(allele))
                .findFirst()
                .orElseThrow(() -> new IllegalArgumentException(String.format("unknown allele '%s'", allele)));
    }

    /**
     * Returns the value whose state is the same as the one provided.
     * @param state the query state.
     * @return never {@code null}.
     * @throws IllegalArgumentException if the input state is {@code null} or there is no
     * value with such a state.
     */
    public static CNVAllele valueOf(final CopyNumberTriState state) {
        Utils.nonNull(state, "the input state cannot be null");
        return Stream.of(values()).filter(v -> v.state.equals(state)).findFirst().orElseThrow(IllegalArgumentException::new);
    }
}
