package org.broadinstitute.hellbender.tools.exome;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionData;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionInitializer;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionState;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Detects copy-number events using allelic-count data and GATK CNV output.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Detect copy-number events in a tumor sample using allelic-count data and GATK CNV output. " +
                "Allelic-count data (reference/alternate counts from the GetHetCoverage tool) is segmented using " +
                "circular binary segmentation; the result is combined with the target coverages " +
                "and called segments found by the GATK CNV tool. Bayesian parameter estimation of models for the " +
                "copy ratios and minor allele fractions in each segment is performed using Markov chain Monte Carlo.",
        oneLineSummary = "Detect copy-number events using allelic-count data and GATK CNV output.",
        programGroup = CopyNumberProgramGroup.class
)
public class AllelicCNV extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1l;

    private static final double INITIAL_ALLELIC_BIAS_GUESS = 1.;

    //filename tags for output
    protected static final String SNP_MAF_SEG_FILE_TAG = "MAF";
    protected static final String UNION_SEG_FILE_TAG = "union";
    protected static final String SMALL_MERGED_SEG_FILE_TAG = "no-small";
    protected static final String INITIAL_SEG_FILE_TAG = "sim-0";
    protected static final String INTERMEDIATE_SEG_FILE_TAG = "sim";
    protected static final String FINAL_SEG_FILE_TAG = "sim-final";
    protected static final String GATK_SEG_FILE_TAG = "cnv";
    protected static final String CGA_ACS_SEG_FILE_TAG = "acs";

    //CLI arguments
    protected static final String OUTPUT_PREFIX_LONG_NAME = "outputPrefix";
    protected static final String OUTPUT_PREFIX_SHORT_NAME = "pre";

    protected static final String SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD_LONG_NAME = "smallSegmentThreshold";
    protected static final String SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD_SHORT_NAME = "smallTh";

    protected static final String NUM_SAMPLES_COPY_RATIO_LONG_NAME = "numSamplesCopyRatio";
    protected static final String NUM_SAMPLES_COPY_RATIO_SHORT_NAME = "numSampCR";

    protected static final String NUM_BURN_IN_COPY_RATIO_LONG_NAME = "numBurnInCopyRatio";
    protected static final String NUM_BURN_IN_COPY_RATIO_SHORT_NAME = "numBurnCR";

    protected static final String NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME = "numSamplesAlleleFraction";
    protected static final String NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME = "numSampAF";

    protected static final String NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME = "numBurnInAlleleFraction";
    protected static final String NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME = "numBurnAF";

    protected static final String INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME = "intervalThresholdCopyRatio";
    protected static final String INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME = "simThCR";

    protected static final String INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME = "intervalThresholdAlleleFraction";
    protected static final String INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME = "simThAF";

    protected static final String MAX_NUM_SNP_SEGMENTATION_ITERATIONS_LONG_NAME = "maxNumIterationsSNPSeg";
    protected static final String MAX_NUM_SNP_SEGMENTATION_ITERATIONS_SHORT_NAME = "maxIterSNP";

    protected static final String MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME = "maxNumIterationsSimSeg";
    protected static final String MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_SHORT_NAME = "maxIterSim";

    protected static final String USE_ALL_COPY_RATIO_SEGMENTS_LONG_NAME = "useAllCopyRatioSegments";
    protected static final String USE_ALL_COPY_RATIO_SEGMENTS_SHORT_NAME = "useAllCRSeg";

    @Argument(
            doc = "Input file for tumor-sample ref/alt read counts at normal-sample heterozygous-SNP sites (output of GetHetCoverage tool).",
            fullName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TUMOR_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File snpCountsFile;

    @Argument(
            doc = "Input file for tumor-sample tangent-normalized target coverages (.tn.tsv output of GATK CNV tool).",
            fullName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.TANGENT_NORMALIZED_COUNTS_FILE_SHORT_NAME,
            optional = false
    )
    protected File tangentNormalizedCoverageFile;

    @Argument(
            doc = "Input file for tumor-sample target-coverage segments with calls (output of GATK CNV tool).",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            optional = false
    )
    protected File targetSegmentsFile;

    @Argument(
            doc = "Input file for allelic-bias panel of normals.",
            fullName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.ALLELIC_PON_FILE_SHORT_NAME,
            optional = true
    )
    protected File allelicPONFile;

    @Argument(
            doc = "Prefix for output files. Will also be used as the sample name if that is not provided." +
                    "(Note: if this is a file path or contains slashes (/), " +
                    "the string after the final slash will be used as the sample name if that is not provided.)",
            fullName = OUTPUT_PREFIX_LONG_NAME,
            shortName = OUTPUT_PREFIX_SHORT_NAME,
            optional = false
    )
    protected String outputPrefix;

    @Argument(
            doc = "Threshold for small-segment merging. If a segment has strictly less than this number of targets, " +
                    "it is considered small and will be merged with an adjacent segment.",
            fullName = SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD_LONG_NAME,
            shortName = SMALL_SEGMENT_TARGET_NUMBER_THRESHOLD_SHORT_NAME,
            optional = true
    )
    protected int smallSegmentTargetNumberThreshold = 3;

    @Argument(
            doc = "Total number of MCMC samples for copy-ratio model.",
            fullName = NUM_SAMPLES_COPY_RATIO_LONG_NAME,
            shortName = NUM_SAMPLES_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected int numSamplesCopyRatio = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for copy-ratio model.",
            fullName = NUM_BURN_IN_COPY_RATIO_LONG_NAME,
            shortName = NUM_BURN_IN_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected int numBurnInCopyRatio = 50;

    @Argument(
            doc = "Total number of MCMC samples for allele-fraction model.",
            fullName = NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected int numSamplesAlleleFraction = 200;

    @Argument(
            doc = "Number of burn-in samples to discard for allele-fraction model.",
            fullName = NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected int numBurnInAlleleFraction = 100;

    @Argument(
            doc = "Number of 95% credible-interval widths to use for copy-ratio similar-segment merging.",
            fullName = INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME,
            shortName = INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME,
            optional = true
    )
    protected double intervalThresholdCopyRatio = 2.;

    @Argument(
            doc = "Number of 95% credible-interval widths to use for allele-fraction similar-segment merging.",
            fullName = INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME,
            shortName = INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME,
            optional = true
    )
    protected double intervalThresholdAlleleFraction = 2.;

    @Argument(
            doc = "Maximum number of iterations allowed for SNP segmentation.",
            fullName = MAX_NUM_SNP_SEGMENTATION_ITERATIONS_LONG_NAME,
            shortName = MAX_NUM_SNP_SEGMENTATION_ITERATIONS_SHORT_NAME,
            optional = true
    )
    protected int maxNumSNPSegmentationIterations = 25;

    @Argument(
            doc = "Maximum number of iterations allowed for similar-segment merging.",
            fullName = MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_LONG_NAME,
            shortName = MAX_NUM_SIMILAR_SEGMENT_MERGING_ITERATIONS_SHORT_NAME,
            optional = true
    )
    protected int maxNumSimilarSegmentMergingIterations = 25;

    @Argument(
            doc = "Enable use of all copy-ratio--segment breakpoints. " +
                    "(Default behavior uses only breakpoints from segments not called copy neutral, " +
                    "if calls are available in output of GATK CNV provided, and none otherwise.)",
            fullName = USE_ALL_COPY_RATIO_SEGMENTS_LONG_NAME,
            shortName = USE_ALL_COPY_RATIO_SEGMENTS_SHORT_NAME,
            optional = true
    )
    protected boolean useAllCopyRatioSegments = false;

    @Override
    protected void runPipeline(final JavaSparkContext ctx) {
        final String originalLogLevel =
                (ctx.getLocalProperty("logLevel") != null) ? ctx.getLocalProperty("logLevel") : "INFO";
        ctx.setLogLevel("WARN");

        //the string after the final slash in the output prefix (which may be an absolute file path) will be used as the sample name
        final String sampleName = outputPrefix.substring(outputPrefix.lastIndexOf("/") + 1);

        logger.info("Starting workflow for " + sampleName + "...");

        //make Genome from input target coverages and SNP counts
        logger.info("Loading input files...");
        final Genome genome = new Genome(tangentNormalizedCoverageFile, snpCountsFile, sampleName);

        //load allelic-bias panel of normals if provided
        final AllelicPanelOfNormals allelicPON =
                allelicPONFile != null ? new AllelicPanelOfNormals(allelicPONFile) : AllelicPanelOfNormals.EMPTY_PON;

        //load target-coverage segments from input file
        final List<ModeledSegment> targetSegmentsWithCalls =
                SegmentUtils.readModeledSegmentsFromSegmentFile(targetSegmentsFile);
        logger.info("Number of target-coverage segments from CNV output: " + targetSegmentsWithCalls.size());

        //merge copy-neutral and uncalled segments (unless disabled) and fix up target-segment start breakpoints
        logger.info("Preparing target-coverage segments...");
        final List<SimpleInterval> targetSegments = prepareTargetSegments(genome, targetSegmentsWithCalls);

        //segment SNPs using CBS on per-SNP MLE minor allele fraction iteratively until convergence
        //(final segment file is output as a side effect)
        logger.info("Performing SNP segmentation...");
        final List<SimpleInterval> snpSegments = performSNPSegmentationStep(sampleName, genome);

        //combine SNP and target-coverage segments
        logger.info("Combining SNP and target-coverage segments...");
        final List<SimpleInterval> unionedSegments = SegmentUtils.unionSegments(targetSegments, snpSegments, genome);
        final File unionedSegmentsFile = new File(outputPrefix + "-" + UNION_SEG_FILE_TAG + ".seg");
        SegmentUtils.writeSegmentFileWithNumTargetsAndNumSNPs(unionedSegmentsFile, unionedSegments, genome);

        //small-segment merging (note that X and Y are always small segments and dropped, since GATK CNV drops them)
        logger.info("Merging small segments...");
        final SegmentedModel segmentedModelWithSmallSegments = new SegmentedModel(unionedSegments, genome);
        final SegmentedModel segmentedModel = segmentedModelWithSmallSegments.mergeSmallSegments(smallSegmentTargetNumberThreshold);
        final File segmentedModelFile = new File(outputPrefix + "-" + SMALL_MERGED_SEG_FILE_TAG + ".seg");
        segmentedModel.writeSegmentFileWithNumTargetsAndNumSNPs(segmentedModelFile);

        //initial MCMC model fitting performed by ACNVModeller constructor
        final ACNVModeller modeller = new ACNVModeller(segmentedModel, allelicPON,
                numSamplesCopyRatio, numBurnInCopyRatio, numSamplesAlleleFraction, numBurnInAlleleFraction, ctx);

        final File initialModeledSegmentsFile = new File(outputPrefix + "-" + INITIAL_SEG_FILE_TAG + ".seg");
        modeller.writeACNVModeledSegmentFile(initialModeledSegmentsFile);

        //similar-segment merging (segment files are output for each merge iteration)
        logger.info("Merging similar segments...");
        performSimilarSegmentMergingStep(modeller);

        //write final model fit to file
        final File finalModeledSegmentsFile = new File(outputPrefix + "-" + FINAL_SEG_FILE_TAG + ".seg");
        modeller.writeACNVModeledSegmentFile(finalModeledSegmentsFile);

        //write file for GATK CNV formatted seg file
        final File finalModeledSegmentsFileAsGatkCNV = new File(outputPrefix + "-" + FINAL_SEG_FILE_TAG + "." + GATK_SEG_FILE_TAG + ".seg");
        modeller.writeModeledSegmentFile(finalModeledSegmentsFileAsGatkCNV);

        //write file for ACS-compatible output to help Broad CGA
        final File finalACSModeledSegmentsFile = new File(outputPrefix + "-" + FINAL_SEG_FILE_TAG + "." + CGA_ACS_SEG_FILE_TAG + ".seg");
        modeller.writeAllelicCapSegFile(finalACSModeledSegmentsFile);

        ctx.setLogLevel(originalLogLevel);
        logger.info("SUCCESS: Allelic CNV run complete for sample " + sampleName + ".");
    }

    //merge copy-neutral and uncalled segments (unless disabled) and fix up target-segment start breakpoints
    private List<SimpleInterval> prepareTargetSegments(final Genome genome,
                                                       final List<ModeledSegment> targetSegmentsWithCalls) {
        final List<SimpleInterval> targetSegmentsWithUnfixedStarts;
        if (!useAllCopyRatioSegments) {
            logger.info("Merging copy-neutral and uncalled target-coverage segments...");
            targetSegmentsWithUnfixedStarts = SegmentMergeUtils.mergeNeutralSegments(targetSegmentsWithCalls);
            logger.info("Number of segments after copy-neutral merging: " + targetSegmentsWithUnfixedStarts.size());
        } else {
            targetSegmentsWithUnfixedStarts = targetSegmentsWithCalls.stream().map(ModeledSegment::getSimpleInterval)
                    .collect(Collectors.toList());
        }
        //fix up target-segment start breakpoints (convert from target-end--target-end to target-start--target-end)
        return SegmentUtils.fixTargetSegmentStarts(targetSegmentsWithUnfixedStarts, genome.getTargets());
    }

    //segment SNPs using CBS on per-SNP MLE minor allele fraction iteratively until convergence
    //(final segment file is output as a side effect)
    private List<SimpleInterval> performSNPSegmentationStep(final String sampleName, final Genome genome) {
        //initial segmentation of SNPs on per-SNP MLE minor allele fraction, assuming no allelic bias
        //(segments are written to temporary file)
        logger.info("Performing initial SNP segmentation (assuming no allelic bias)...");
        final File snpSegmentFile;
        try {
            snpSegmentFile = File.createTempFile(outputPrefix + "-" + SNP_MAF_SEG_FILE_TAG + "-tmp", ".seg");
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not create temporary SNP segmentation file.", e);
        }
        List<SimpleInterval> snpSegments =
                calculateAndWriteSNPSegmentation(sampleName, genome, snpSegmentFile, INITIAL_ALLELIC_BIAS_GUESS);
        logger.info("Initial number of SNP segments: " + snpSegments.size());
        final Set<List<SimpleInterval>> snpSegmentationsFound = new HashSet<>();
        //perform iterations of SNP segmentation until convergence
        for (int numIterations = 1; numIterations <= maxNumSNPSegmentationIterations; numIterations++) {
            logger.info("SNP-segmentation iteration: " + numIterations);
            snpSegmentationsFound.add(snpSegments);
            //use AlleleFractionInitializer to determine mean of global allelic-bias distribution given current SNP segmentation
            final double allelicBias = calculateAllelicBias(genome, snpSegmentFile);
            logger.info(String.format("Mean allelic bias: %.3f", allelicBias));
            //resegment SNPs on per-SNP MLE minor allele fraction assuming mean allelic bias found by AlleleFractionInitializer
            //(segments are written to temporary file)
            snpSegments = calculateAndWriteSNPSegmentation(sampleName, genome, snpSegmentFile, allelicBias);
            logger.info("Number of SNP segments: " + snpSegments.size());
            //stop if a previously found SNP segmentation is found again
            if (snpSegmentationsFound.contains(snpSegments)) {
                //take final SNP segmentation to be the same as the last one found
                logger.info("Final number of SNP segments: " + snpSegments.size());
                final File finalSNPSegmentFile = new File(outputPrefix + "-" + SNP_MAF_SEG_FILE_TAG + ".seg");
                snpSegments = calculateAndWriteSNPSegmentation(sampleName, genome, finalSNPSegmentFile, allelicBias);
                break;
            }
        }
        return snpSegments;
    }

    //return SNP segments from CBS (writes file as a side effect)
    private List<SimpleInterval> calculateAndWriteSNPSegmentation(final String sampleName, final Genome genome,
                                                                  final File snpSegmentFile, final double allelicBias) {
        SNPSegmenter.writeSegmentFile(genome.getSNPs(), sampleName, snpSegmentFile, allelicBias);
        return SegmentUtils.readIntervalsFromSegmentFile(snpSegmentFile);
    }

    //use AlleleFractionInitializer to determine mean of global allelic-bias distribution given SNP segmentation file
    private double calculateAllelicBias(final Genome genome, final File snpSegmentFile) {
        final SegmentedModel segmentedModelForSNPSegmentation = new SegmentedModel(snpSegmentFile, genome);
        final AlleleFractionData data = new AlleleFractionData(segmentedModelForSNPSegmentation);
        final AlleleFractionState initialState = new AlleleFractionInitializer(data).getInitializedState();
        return initialState.meanBias();
    }

    //similar-segment merging (segment files are output for each merge iteration)
    private void performSimilarSegmentMergingStep(final ACNVModeller modeller) {
        logger.info("Initial number of segments before similar-segment merging: " + modeller.getACNVModeledSegments().size());
        //perform iterations of similar-segment merging until all similar segments are merged
        for (int numIterations = 1; numIterations <= maxNumSimilarSegmentMergingIterations; numIterations++) {
            logger.info("Similar-segment merging iteration: " + numIterations);
            final int prevNumSegments = modeller.getACNVModeledSegments().size();
            modeller.performSimilarSegmentMergingIteration(intervalThresholdCopyRatio, intervalThresholdAlleleFraction);
            if (modeller.getACNVModeledSegments().size() == prevNumSegments) {
                break;
            }
            final File modeledSegmentsFile = new File(outputPrefix + "-" + INTERMEDIATE_SEG_FILE_TAG + "-" + numIterations + ".seg");
            modeller.writeACNVModeledSegmentFile(modeledSegmentsFile);
        }
        logger.info("Final number of segments after similar-segment merging: " + modeller.getACNVModeledSegments().size());
    }
}