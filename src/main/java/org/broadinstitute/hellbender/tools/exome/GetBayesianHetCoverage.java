package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.File;

/**
 * Outputs base statistics and likelihoods of heterzygosity and homozygosity at probable heterozygous
 * SNP sites. The SNP sites must be provided.
 *
 * @author Mehrtash Babadi &lt;mehrtash@brgoadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Outputs base statistics and likelihoods of heterzygosity and homozygosity at probable heterozygous SNP sites.",
        oneLineSummary = "Outputs base statistics and likelihoods of heterzygosity and homozygosity at probable heterozygous SNP sites.",
        programGroup = CopyNumberProgramGroup.class
)
public final class GetBayesianHetCoverage extends CommandLineProgram {

    protected static final String READ_DEPTH_THRESHOLD_FULL_NAME = "readDepthThreshold";
    protected static final String READ_DEPTH_THRESHOLD_SHORT_NAME = "readThresh";

    protected static final String MINIMUM_MAPPING_QUALITY_FULL_NAME = "minimumMappingQuality";
    protected static final String MINIMUM_MAPPING_QUALITY_SHORT_NAME = "minMQ";

    protected static final String MINIMUM_BASE_QUALITY_FULL_NAME = "minimumBaseQuality";
    protected static final String MINIMUM_BASE_QUALITY_SHORT_NAME = "minBQ";

    protected static final String HET_CALLING_STRINGENCY_FULL_NAME = "hetCallingStringency";
    protected static final String HET_CALLING_STRINGENCY_SHORT_NAME = "hetS";

    protected static final String QUADRATURE_ORDER_FULL_NAME = "quadratureOrder";
    protected static final String QUADRATURE_ORDER_SHORT_NAME = "quad";

    protected static final String MINIMUM_ABNORMAL_FRACTION_FULL_NAME = "minimumAbnormalFraction";
    protected static final String MINIMUM_ABNORMAL_FRACTION_SHORT_NAME = "minAF";

    protected static final String MAXIMUM_ABNORMAL_FRACTION_FULL_NAME = "maximumAbnormalFraction";
    protected static final String MAXIMUM_ABNORMAL_FRACTION_SHORT_NAME = "maxAF";

    protected static final String MAXIMUM_COPY_NUMBER_FULL_NAME = "maximumCopyNumber";
    protected static final String MAXIMUM_COPY_NUMBER_SHORT_NAME = "maxCN";

    protected static final String ERROR_PROBABILITY_ADJUSTMENT_FACTOR_FULL_NAME = "errorAdjustmentFactor";
    protected static final String ERROR_PROBABILITY_ADJUSTMENT_FACTOR_SHORT_NAME = "errFact";

    @ArgumentCollection
    protected static final ReferenceInputArgumentCollection REFERENCE_ARGUMENTS =
            new RequiredReferenceInputArgumentCollection();

    @Argument(
            doc = "BAM file.",
            fullName = ExomeStandardArgumentDefinitions.GENERIC_BAM_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.GENERIC_BAM_FILE_SHORT_NAME,
            optional = false
    )
    protected File bamFile;

    @Argument(
            doc = "Interval-list file of common SNPs.",
            fullName = ExomeStandardArgumentDefinitions.SNP_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpFile;

    @Argument(
            doc = "Output file for base counts at probable heterozygous SNP sites.",
            fullName = ExomeStandardArgumentDefinitions.GENERIC_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.GENERIC_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File hetOutputFile;

    @Argument(
            doc = "Output file for base counts and statistics at probable heterozygous SNP sites.",
            fullName = ExomeStandardArgumentDefinitions.GENERIC_DETAILED_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.GENERIC_DETAILED_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = true
    )
    protected File hetDetailedOutputFile = null;

    @Argument(
            doc = "Minimum read depth; SNP sites with lower read depths will be ignored.",
            fullName = READ_DEPTH_THRESHOLD_FULL_NAME,
            shortName = READ_DEPTH_THRESHOLD_SHORT_NAME,
            optional = true
    )
    protected int readDepthThreshold = 15;

    @Argument(
            doc = "Minimum mapping quality; reads with lower quality will be filtered out of pileup.",
            shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME,
            fullName  = MINIMUM_MAPPING_QUALITY_FULL_NAME,
            optional = true
    )
    protected int minimumMappingQuality = 30;

    @Argument(
            doc = "Minimum base quality; base calls with lower quality will be filtered out of pileup.",
            shortName = MINIMUM_BASE_QUALITY_SHORT_NAME,
            fullName = MINIMUM_BASE_QUALITY_FULL_NAME,
            optional = true
    )
    protected int minimumBaseQuality = 20;

    @Argument(
            doc = "Validation stringency for all BAM files read by this program.  Setting stringency to SILENT " +
                    "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
                    "do not otherwise need to be decoded.",
            common=true)
    protected ValidationStringency VALIDATION_STRINGENCY = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;

    @Argument(
            doc = "Integration quadrature order; a good choice is the typical read depth",
            fullName = QUADRATURE_ORDER_FULL_NAME,
            shortName = QUADRATURE_ORDER_SHORT_NAME,
            optional = true
    )
    protected int quadratureOrder = 200;

    @Argument(
            doc = "Heterozygosous SNP calling stringency.",
            fullName = HET_CALLING_STRINGENCY_FULL_NAME,
            shortName = HET_CALLING_STRINGENCY_SHORT_NAME,
            optional = false
    )
    protected double hetCallingStringency = 5;

    @Argument(
            doc = "Estimated minimum abnormal cell fraction (for building the allele fraction prior)",
            fullName = MINIMUM_ABNORMAL_FRACTION_FULL_NAME,
            shortName = MINIMUM_ABNORMAL_FRACTION_SHORT_NAME,
            optional = true
    )
    protected double minimumAbnormalFraction = 0.5;

    @Argument(
            doc = "Estimated maximum abnormal cell fraction (for building the allele fraction prior)",
            fullName = MAXIMUM_ABNORMAL_FRACTION_FULL_NAME,
            shortName = MAXIMUM_ABNORMAL_FRACTION_SHORT_NAME,
            optional = true
    )
    protected double maximumAbnormalFraction = 0.8;

    @Argument(
            doc = "Estimated maximum copy number for either homologs in the abnormal portion of the sample (for " +
                    "building the allele fraction prior)",
            fullName = MAXIMUM_COPY_NUMBER_FULL_NAME,
            shortName = MAXIMUM_COPY_NUMBER_SHORT_NAME,
            optional = true
    )
    protected int maximumCopyNumber = 2;

    @Argument(
            doc = "(Experimental) error probability adjustment factor (could be any positive value)",
            fullName = ERROR_PROBABILITY_ADJUSTMENT_FACTOR_FULL_NAME,
            shortName = ERROR_PROBABILITY_ADJUSTMENT_FACTOR_SHORT_NAME,
            optional = true
    )
    protected double errorProbabilityAdjustmentFactor = 1.0;

    @Override
    protected Object doWork() {

        final BayesianHetPulldownCalculator calculator = new BayesianHetPulldownCalculator(REFERENCE_ARGUMENTS.getReferenceFile(),
                snpFile, minimumMappingQuality, minimumBaseQuality, readDepthThreshold, VALIDATION_STRINGENCY,
                quadratureOrder, minimumAbnormalFraction, maximumAbnormalFraction, maximumCopyNumber,
                errorProbabilityAdjustmentFactor);

        logger.info("Calculating Het pulldown...");
        final Pulldown hetPulldown = calculator.getHetPulldown(bamFile, hetCallingStringency);

        hetPulldown.writeWithBaseCounts(hetOutputFile);
        logger.info("Het pulldown with base counts written to " + hetOutputFile.toString());

        if (hetDetailedOutputFile != null) {
            hetPulldown.writeWithBaseCountsAndLikelihoods(hetDetailedOutputFile);
            logger.info("Het pulldown with base counts and statistics written to " + hetDetailedOutputFile.toString());
        }

        return "SUCCESS";
    }
}