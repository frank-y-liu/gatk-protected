package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.exome.CNVAllele;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by valentin on 3/30/16.
 */
public class VariantEvaluationContextBuilder extends VariantContextBuilder {

    private int[] truthAC = new int[CNVAllele.values().length];
    private int[] callsAC = new int[CNVAllele.values().length];

    public VariantEvaluationContextBuilder() {
        super();
    }

    public VariantEvaluationContextBuilder(final VariantContext init) {
        super(init);
        alleles(init.getAlleles());
        genotypes(init.getGenotypes());
    }

    @Override
    public VariantEvaluationContextBuilder alleles(final Collection<Allele> alleles) {
        ParamUtils.noNulls(alleles, "the input cannot contain null alleles");
        validateAlleleStrings(alleles.stream().map(Allele::getBaseString).collect(Collectors.toList()));
        if (!alleles.iterator().next().isReference()) {
            throw new IllegalArgumentException("the first allele must be a reference allele");
        }
        super.alleles(alleles);
        return this;
    }

    private void validateAlleleStrings(final List<String> allelesString) {
        ParamUtils.noNulls(allelesString, "the input cannot contain null strings");
        if (allelesString.isEmpty()) {
            throw new IllegalArgumentException("the allele list cannot be empty");
        } if (!allelesString.get(0).equals(CNVAllele.REF.allele.getBaseString())) {
            throw new IllegalArgumentException("the first allele must be the reference allele: " + CNVAllele.REF.allele);
        }
        final Set<CNVAllele> allelesFound = EnumSet.of(CNVAllele.REF);
        for (int i = 1; i < allelesString.size(); i++) {
            final String alleleString = allelesString.get(i);
            final CNVAllele allele;
            try {
                allele = CNVAllele.valueOf(Allele.create(alleleString, false));
            } catch(final IllegalArgumentException ex) {
                throw new IllegalArgumentException("unknown allele with string: " + alleleString);
            }
            if (!allelesFound.add(allele)) {
                throw new IllegalArgumentException("an allele cannot be listed twice.");
            }
        }
    }

    @Override
    public VariantEvaluationContextBuilder alleles(final List<String> alleleStrings) {
        validateAlleleStrings(alleleStrings);
        super.alleles(alleleStrings);
        return this;
    }

    @Override
    public VariantEvaluationContextBuilder alleles(final String... alleleStrings) {
        validateAlleleStrings(Arrays.asList(alleleStrings));
        super.alleles(alleleStrings);
        return this;
    }

    @Override
    public List<Allele> getAlleles() {
        final List<Allele> result = super.getAlleles();
        if (result == null) {
            throw new IllegalStateException("you must set the alleles before calling getAlleles");
        }
        return result;
    }

    @Override
    public VariantEvaluationContextBuilder genotypes(final GenotypesContext genotypes) {
        super.genotypes(genotypes);
        return this;
    }

    @Override
    public VariantEvaluationContextBuilder genotypesNoValidation(final GenotypesContext genotypes) {
        throw new UnsupportedOperationException("non-validated genotype assignation is not available");
    }

    @Override
    public VariantEvaluationContextBuilder genotypes(final Collection<Genotype> genotypes) {
        super.genotypes(genotypes);
        calculateAFs(genotypes);
        return this;
    }

    @Override
    public VariantEvaluationContextBuilder genotypes(final Genotype ... genotypes) {
        super.genotypes(genotypes);
        calculateAFs(Arrays.asList(genotypes));
        return this;
    }

    @Override
    public VariantEvaluationContextBuilder noGenotypes() {
        super.noGenotypes();
        Arrays.fill(truthAC, 0);
        Arrays.fill(callsAC, 0);
        return this;
    }

    private void calculateAFs(final Iterable<Genotype> genotypes) {
        final int[] truthAC = new int[CNVAllele.values().length];
        final int[] callsAC = new int[CNVAllele.values().length];
        for (final Genotype genotype : genotypes) {
            final List<Allele> alleles = genotype.getAlleles();
            if (alleles.size() > 1) {
                throw new GATKException("unexpected CNV genotype ploidy: " + alleles.size());
            } else if (!alleles.isEmpty()) {
                final int index = CNVAllele.valueOf(alleles.get(0)).ordinal();
                callsAC[index]++;
            }
            final String truthGT = String.valueOf(genotype.getExtendedAttribute(VariantEvaluationContext.TRUTH_GENOTYPE_KEY, VCFConstants.MISSING_VALUE_v4));
            final int truthAlleleIndex = truthGT.equals(VCFConstants.MISSING_VALUE_v4) ? -1 : Integer.parseInt(truthGT);
            if (truthAlleleIndex >= 0) {
                final List<Allele> contextAlleles = getAlleles();
                if (truthAlleleIndex >= contextAlleles.size()) {
                    throw new GATKException("unexpected CNV truth genotype makes reference to an unexistent allele: " + truthGT);
                }
                truthAC[CNVAllele.valueOf(contextAlleles.get(truthAlleleIndex)).ordinal()]++;
            }
            genotype.getAllele(0);
        }
        this.truthAC = truthAC;
        this.callsAC = callsAC;
    }

    @Override
    public VariantEvaluationContext make(final boolean leaveItModifiable) {
        final List<Allele> contextAlleles = getAlleles();
        final int[] alleleToCNVAlleleOrdinal = contextAlleles.stream()
                .mapToInt(a -> CNVAllele.valueOf(a).ordinal()).toArray();
        final int callsAN = (int) MathUtils.sum(this.callsAC);
        final int truthAN = (int) MathUtils.sum(this.truthAC);
        attribute(VariantEvaluationContext.CALLS_ALLELE_NUMBER_KEY, callsAN);
        attribute(VariantEvaluationContext.TRUTH_ALLELE_NUMBER_KEY, truthAN);
        if (callsAN > 0) {
            final double[] callsAF = IntStream.range(1, contextAlleles.size())
                    .mapToDouble(i -> callsAC[alleleToCNVAlleleOrdinal[i]] / (double) callsAN).toArray();
            attribute(VariantEvaluationContext.CALLS_ALLELE_FREQUENCY_KEY, callsAF);
        }
        if (truthAN > 0) {
            final double[] truthAF = IntStream.range(1, contextAlleles.size())
                    .mapToDouble(i -> truthAC[alleleToCNVAlleleOrdinal[i]] / (double) truthAN).toArray();
            attribute(VariantEvaluationContext.TRUTH_ALLELE_FREQUENCY_KEY, truthAF);
        }
        return new VariantEvaluationContext(super.make(leaveItModifiable));
    }

    @Override
    public VariantEvaluationContext make() {
        return make(false);
    }
}
