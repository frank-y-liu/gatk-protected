package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.argumentcollections.ReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.RequiredReferenceInputArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.File;

@CommandLineProgramProperties(
        summary = "Gives read depth at provided intervals.",
        oneLineSummary = "Gives read depth at provided intervals.",
        programGroup = CopyNumberProgramGroup.class
)
public final class GetReadDepth extends CommandLineProgram {

    private final Logger logger = LogManager.getLogger(GetReadDepth.class);

    protected static final String READ_DEPTH_THRESHOLD_FULL_NAME = "readDepthThreshold";
    protected static final String READ_DEPTH_THRESHOLD_SHORT_NAME = "readThresh";

    protected static final String MINIMUM_MAPPING_QUALITY_FULL_NAME = "minimumMappingQuality";
    protected static final String MINIMUM_MAPPING_QUALITY_SHORT_NAME = "minMQ";

    protected static final String MINIMUM_BASE_QUALITY_FULL_NAME = "minimumBaseQuality";
    protected static final String MINIMUM_BASE_QUALITY_SHORT_NAME = "minBQ";

    @ArgumentCollection
    protected static final ReferenceInputArgumentCollection REFERENCE_ARGUMENTS =
            new RequiredReferenceInputArgumentCollection();

    @Argument(
            doc = "BAM file.",
            fullName = ExomeStandardArgumentDefinitions.BAM_FILE_SHORT_NAME,
            shortName = ExomeStandardArgumentDefinitions.BAM_FILE_LONG_NAME,
            optional = false
    )
    protected File bamFile;

    @Argument(
            doc = "IntervalList file.",
            fullName = ExomeStandardArgumentDefinitions.SNP_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SNP_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpFile;

    @Argument(
            doc = "Output file for read depths.",
            fullName = "readDepthFile",
            shortName = "RD",
            optional = false
    )
    protected File outputFile;

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

    @Override
    protected Object doWork() {

        BayesianHetPulldownCalculator calculator;

        calculator = new BayesianHetPulldownCalculator(REFERENCE_ARGUMENTS.getReferenceFile(),
                IntervalList.fromFile(snpFile), minimumMappingQuality, minimumBaseQuality, readDepthThreshold,
                VALIDATION_STRINGENCY, 1.0);

        Pulldown readDepthPulldown = calculator.getReadDepthPulldown(bamFile);

        logger.info("Writing read depth pulldown to " + outputFile.toString());
        readDepthPulldown.writeReadDepth(outputFile);

        return "SUCCESS";

    }
}