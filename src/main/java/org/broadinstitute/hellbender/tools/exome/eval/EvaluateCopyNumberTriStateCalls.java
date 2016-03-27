package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.parquet.it.unimi.dsi.fastutil.objects.Object2IntLinkedOpenHashMap;
import org.apache.parquet.it.unimi.dsi.fastutil.objects.Object2IntMap;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Tool to evaluate the output of {@link DiscoverCopyNumberTriStateSegments}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Evaluate a set of call segments against the truth segments",
        oneLineSummary = "Evaluate a set of germline call segments",
        programGroup = CopyNumberProgramGroup.class
)
public final class EvaluateCopyNumberTriStateCalls extends CommandLineProgram {

    public static final String CALLS_FILE_SHORT_NAME = "calls";
    public static final String CALLS_FILE_FULL_NAME = "callsFile";
    public static final String TRUTH_FILE_SHORT_NAME = "truth";
    public static final String TRUTH_FILE_FULL_NAME = "truthFile";
    public static final String SAMPLES_SHORT_NAME = "sample";
    public static final String SAMPLES_FULL_NAME = "sample";
    public static final String SAMPLES_LIST_SHORT_NAME = "samples";
    public static final String SAMPLES_LIST_FULL_NAME = "sampleList";
    public static final String SAMPLE_SUMMARY_OUTPUT_SHORT_NAME = "summary";
    public static final String SAMPLE_SUMMARY_OUTPUT_FULL_NAME = "sampleSummaryOutput";
    public static final String DETAIL_CALL_OUTPUT_SHORT_NAME = "sites";
    public static final String DETAIL_CALL_OUTPUT_FULL_NAME = "siteDetailsOutput";
    public static final String OVERALL_SUMMARY_LINE_SHORT_NAME = "includeOverall";
    public static final String OVERALL_SUMMARY_LINE_FULL_NAME = "includeOverallSummaryOutputLine";
    public static final String OVERALL_SUMMARY_SAMPLE_SHORT_NAME = "overallSample";
    public static final String OVERALL_SUMMARY_SAMPLE_FULL_NAME = "overallSampleName";
    public static final String FREQUENCY_SMOOTHING_SHORT_NAME = "smooth";
    public static final String FREQUENCY_SMOOTHING_FULL_NAME = "frequencySmoothing";

    public static final String DEFAULT_OVERALL_SUMMARY_SAMPLE_NAME = "ALL";
    public static final double DEFAULT_FREQUENCY_SMOOTHING = 1.0;

    @ArgumentCollection
    protected TargetArgumentCollection targetArguments = new TargetArgumentCollection();

    @ArgumentCollection
    protected CallFiltersCollection filterArguments = new CallFiltersCollection();

    @Argument(
            doc = "File contained the called segments",
            shortName = CALLS_FILE_SHORT_NAME,
            fullName = CALLS_FILE_FULL_NAME)
    protected File callsFile;

    @Argument(
            doc = "File containing true events",
            shortName = TRUTH_FILE_SHORT_NAME,
            fullName = TRUTH_FILE_FULL_NAME,
            optional = true)
    protected File truthFile;

    @Argument(
            doc = "Output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME
    )
    protected File outputFile;

    @Argument(
            doc = "Samples to evaluate",
            shortName = SAMPLES_SHORT_NAME,
            fullName = SAMPLES_FULL_NAME,
            optional = true
    )
    protected Set<String> samples = new LinkedHashSet<>();

    @Argument(
            doc = "Samples list to evaluate",
            shortName = SAMPLES_LIST_SHORT_NAME,
            fullName  = SAMPLES_LIST_FULL_NAME,
            optional = true
    )
    protected File sampleListFile = null;

    @Argument(
            doc = "Sample summary output",
            shortName = SAMPLE_SUMMARY_OUTPUT_SHORT_NAME,
            fullName = SAMPLE_SUMMARY_OUTPUT_FULL_NAME,
            optional = true
    )
    protected File sampleSummaryOutputFile = null;

    @Argument(
            doc = "Site detail output",
            shortName = DETAIL_CALL_OUTPUT_SHORT_NAME,
            fullName = DETAIL_CALL_OUTPUT_FULL_NAME,
            optional = true
    )
    protected File segmentDetailOutputFile = null;

    @Argument(
            doc = "Whether to include the overall summary line in the sample summary output file",
            shortName = OVERALL_SUMMARY_LINE_SHORT_NAME,
            fullName = OVERALL_SUMMARY_LINE_FULL_NAME,
            optional = true
    )
    protected boolean includeOverallSummaryRecord = false;

    @Argument(
            doc = "Special sample name used for the overall summary output line",
            shortName = OVERALL_SUMMARY_SAMPLE_SHORT_NAME,
            fullName = OVERALL_SUMMARY_SAMPLE_FULL_NAME,
            optional = true
    )
    protected String overallSampleSummaryRecordSampleName = DEFAULT_OVERALL_SUMMARY_SAMPLE_NAME;

    @Argument(
            doc = "Reference copy number",
            shortName = ConvertGSVariantsToSegments.NEUTRAL_COPY_NUMBER_SHORT_NAME,
            fullName = ConvertGSVariantsToSegments.NEUTRAL_COPY_NUMBER_FULL_NAME,
            optional = true
    )
    protected int truthNeutralCopyNumber = ConvertGSVariantsToSegments.NEUTRAL_COPY_NUMBER_DEFAULT;

    @Argument(
            doc = "Frequency smoothing constant",
            shortName = FREQUENCY_SMOOTHING_SHORT_NAME,
            fullName = FREQUENCY_SMOOTHING_FULL_NAME,
            optional = true
    )
    protected double frequencySmoothing = DEFAULT_FREQUENCY_SMOOTHING;

    @Override
    protected Object doWork() {
        final TargetCollection<Target> targets = targetArguments.readTargetCollection(false);
        final VCFFileReader truthReader = openVCFReader(truthFile);
        final VCFFileReader callsReader = openVCFReader(callsFile);
        final Set<String> samples = this.samples.isEmpty() ? composeSetOfSamplesToEvaluate(callsReader) : this.samples;
        final VariantContextWriter outputWriter = openVCFWriter(outputFile, samples);
        final Map<String, EvaluationSampleSummaryRecord> sampleStats = samples.stream()
                .collect(Collectors.toMap(s -> s, EvaluationSampleSummaryRecord::new));
        final List<SimpleInterval> intervals = composeListOfProcessingIntervalsFromInputs(truthReader, callsReader);
        for (final SimpleInterval interval : intervals) {
            for (final VariantEvaluationContext vc : processInterval(truthReader, callsReader, interval, targets)) {
                outputWriter.add(vc);
                updateSampleStats(sampleStats, vc);
            }
        }
        truthReader.close();
        callsReader.close();
        outputWriter.close();
        writeSampleSummaryFile(sampleSummaryOutputFile, sampleStats);
        return "SUCCESS";
    }

    private void writeSampleSummaryFile(final File outputFile, final Map<String, EvaluationSampleSummaryRecord> sampleStats) {
        if (sampleSummaryOutputFile == null) {
            return;
        }
        try (final EvaluationSampleSummaryWriter writer = new EvaluationSampleSummaryWriter(outputFile,  overallSampleSummaryRecordSampleName)) {
            for (final EvaluationSampleSummaryRecord record : sampleStats.values()) {
                writer.writeRecord(record);
            }
            writeOverallSummaryRecordIfApplies(writer);
        } catch (final IOException ex) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, ex);
        }

    }

    private void updateSampleStats(final Map<String, EvaluationSampleSummaryRecord> sampleStats, final VariantEvaluationContext vc) {
        // ignore filtered out variants.
        if (vc.getFilters().stream().anyMatch(f -> !VCFConstants.PASSES_FILTERS_v4.equals(f))) {
            return;
        }
        for (final Genotype genotype : vc.getGenotypes()) {
            final String[] filters = GATKProtectedVariantContextUtils.getAttributeAsStringArray(genotype, VCFConstants.GENOTYPE_FILTER_KEY, () -> new String[] {VCFConstants.PASSES_FILTERS_v4}, "");
            // ignore filtered genotypes.
            if (Stream.of(filters)
                    .anyMatch(s -> !s.isEmpty() && !s.equals(VCFConstants.PASSES_FILTERS_v3))) {
                continue;
            }
            final String sample = genotype.getSampleName();
            final String evalClassString = GATKProtectedVariantContextUtils.getAttributeAsString(genotype, VariantEvaluationContext.EVALUATION_CLASS_KEY, null);
            if (evalClassString != null) {
                final EvaluationClass evalClass = EvaluationClass.fromAcronym(evalClassString);
                sampleStats.get(sample).increase(evalClass);
            }
        }
    }

    private VariantContextWriter openVCFWriter(final File outputFile, final Set<String> samples) {
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        builder.setOutputFile(outputFile);
        final VariantContextWriter result = builder.build();
        final VCFHeader header = new VCFHeader(Collections.emptySet(), samples);
        CNVAllele.addHeaderLinesTo(header);
        EvaluationClass.addHeaderLinesTo(header);

        // Format annotations.
        header.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.Character, "Called genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(VariantEvaluationContext.CALL_QUALITY_KEY, 1, VCFHeaderLineType.Float, "Quality of the call"));
        header.addMetaDataLine(new VCFFormatHeaderLine(VariantEvaluationContext.CALLED_SEGMENTS_COUNT_KEY, 1, VCFHeaderLineType.Integer, "Number of called segments that overlap with the truth"));
        header.addMetaDataLine(new VCFFormatHeaderLine(VariantEvaluationContext.CALLED_ALLELE_COUNTS_KEY, VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Called allele count for mixed calls"));
        header.addMetaDataLine(new VCFFormatHeaderLine(VariantEvaluationContext.TRUTH_COPY_FRACTION_KEY, 1, VCFHeaderLineType.Float, "Truth copy fraction estimated"));
        header.addMetaDataLine(new VCFFormatHeaderLine(VariantEvaluationContext.TRUTH_QUALITY_KEY, 1, VCFHeaderLineType.Float, "Truth call quality"));
        header.addMetaDataLine(new VCFFormatHeaderLine(VariantEvaluationContext.EVALUATION_CLASS_KEY, 1, VCFHeaderLineType.Character, "The evaluation class for the call or lack of call. It the values of the header key '" + EvaluationClass.VCF_HEADER_KEY + "'"));
        header.addMetaDataLine(new VCFFormatHeaderLine(VariantEvaluationContext.TRUTH_GENOTYPE_KEY, 1, VCFHeaderLineType.Character, "The truth genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine(VariantEvaluationContext.CALLED_SEGMENTS_LENGTH_KEY, 1, VCFHeaderLineType.Integer, "Number of targets covered by called segments"));
        header.addMetaDataLine(new VCFFormatHeaderLine(VariantEvaluationContext.CALL_QUALITY_KEY, 1, VCFHeaderLineType.Float, "The quality of the call (the maximum if ther are more than one segment"));

        // Info annotations.
        header.addMetaDataLine(new VCFInfoHeaderLine(VariantEvaluationContext.TRUTH_ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "The frequency of the alternative alleles in the truth callset"));
        header.addMetaDataLine(new VCFInfoHeaderLine(VariantEvaluationContext.TRUTH_ALLELE_NUMBER_KEY, 1, VCFHeaderLineType.Integer, "Total number of called alleles in the truth callset"));
        header.addMetaDataLine(new VCFInfoHeaderLine(VariantEvaluationContext.CALLS_ALLELE_FREQUENCY_KEY, VCFHeaderLineCount.A, VCFHeaderLineType.Float, "The frequency of the alternative alleles in the actual callset"));
        header.addMetaDataLine(new VCFInfoHeaderLine(VariantEvaluationContext.CALLS_ALLELE_NUMBER_KEY, 1, VCFHeaderLineType.Integer, "Total number of called alleles in the actual callset"));
        header.addMetaDataLine(new VCFInfoHeaderLine(VariantEvaluationContext.TRUTH_SEGMENT_LENGTH_KEY, 1, VCFHeaderLineType.Integer, "Number of targets overlapped by this variant"));

        // Filter annotations.
        for (final EvaluationSegmentFilter filter : EvaluationSegmentFilter.values()) {
            header.addMetaDataLine(new VCFFilterHeaderLine(filter.name(), filter.description));
            header.addMetaDataLine(new VCFFilterHeaderLine(filter.acronym, filter.description));
        }
        result.writeHeader(header);
        return result;
    }

    /**
     * Processes a cluster of truth and called variants that may overlap over a genome region.
     * @param truthReader reader to the truth variants.
     * @param callsReader reader to the called variants.
     * @param interval the interval to analyze.
     * @return never {@code null}, a least of evaluation variant records.
     */
    private List<VariantEvaluationContext> processInterval(final VCFFileReader truthReader, final VCFFileReader callsReader,
                                                 final SimpleInterval interval, final TargetCollection<Target> targets) {
        final List<VariantContext> truthVariants = variantQueryToList(truthReader, interval);
        final List<VariantContext> callsVariants = variantQueryToList(callsReader, interval);
        final List<VariantEvaluationContext> evaluatedVariants = new ArrayList<>(truthVariants.size() + callsVariants.size());
        for (final VariantContext truth : truthVariants) {
            final List<VariantContext> overlappingCalls = callsVariants.stream()
                    .filter(vc -> IntervalUtils.overlaps(truth, vc))
                    .collect(Collectors.toList());
            evaluatedVariants.add(composeTruePositive(truth, overlappingCalls, targets));
        }
        for (final VariantContext call : callsVariants) {
            final List<VariantContext> overlappingTruth = callsVariants.stream()
                    .filter(vc -> IntervalUtils.overlaps(call, vc))
                    .collect(Collectors.toList());
            if (overlappingTruth.isEmpty()) {
                evaluatedVariants.add(composeUnknownPositiveVariantContext(call)) ;
            }
        }
        return evaluatedVariants.stream()
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .map(this::applyVariantFilters)
                .collect(Collectors.toList());
    }

    private VariantEvaluationContext applyVariantFilters(final VariantEvaluationContext vec) {
        final Set<EvaluationSegmentFilter> filters = new LinkedHashSet<>();
        if (filterArguments.maximumTruthEventFrequency > 1.0 - vec.getTruthAlleleFrequency(CNVAllele.REF)) {
            filters.add(EvaluationSegmentFilter.CommonEvent);
        }
        if (filterArguments.minimumTruthSegmentLength > vec.getTargetCount()) {
            filters.add(EvaluationSegmentFilter.ShortEvent);
        }
        if (filterArguments.applyMultiAllelicTruthFilter && vec.getTruthAlleleFrequency(CNVAllele.DEL) > 0 &&
                vec.getTruthAlleleFrequency(CNVAllele.DUP) > 0) {
            filters.add(EvaluationSegmentFilter.MultiAllelicTruth);
        }
        final Set<String> filterStrings = filters.isEmpty() ? Collections.singleton(EvaluationSegmentFilter.PASS_ACRONYM) :
                filters.stream().map(Enum<EvaluationSegmentFilter>::name).collect(Collectors.toSet());

        final VariantEvaluationContextBuilder builder = new VariantEvaluationContextBuilder(vec);
        builder.filters(filterStrings);
        return builder.make();
    }

    private VariantEvaluationContext composeUnknownPositiveVariantContext(final VariantContext call) {
        final VariantEvaluationContextBuilder result = new VariantEvaluationContextBuilder();
        result.loc(call.getContig(), call.getStart(), call.getEnd());
        result.id(call.getID());
        result.alleles(call.getAlleles());
        result.genotypes(
                call.getGenotypes().stream()
                .map(this::markDiscoveredGenotypesAsUnknownPositive)
                .collect(Collectors.toList()));
        return result.make();
    }

    private Genotype markDiscoveredGenotypesAsUnknownPositive(final Genotype g) {
        if (!GenotypeCopyNumberTriStateSegments.DISCOVERY_TRUE.equals(g.getExtendedAttribute(GenotypeCopyNumberTriStateSegments.DISCOVERY_KEY))) {
                return g;
        } else {
                return new GenotypeBuilder(g).attribute(VariantEvaluationContext.EVALUATION_CLASS_KEY, EvaluationClass.UNKNOWN_POSITIVE).make();
        }
    }

    private VariantEvaluationContext composeTruePositive(final VariantContext truth, final List<VariantContext> calls, final TargetCollection<Target> targets) {
        final VariantEvaluationContextBuilder builder = new VariantEvaluationContextBuilder();
        builder.loc(truth.getContig(), truth.getStart(), truth.getEnd());
        builder.id(truth.getID());
        builder.alleles(CNVAllele.ALL_ALLELES);
        builder.genotypes(samples.stream()
                .map(s -> composeTruePositiveGenotype(s, truth, calls, targets))
                .collect(Collectors.toList()));
        return builder.make();
    }

    private Genotype composeTruePositiveGenotype(final String sample, final VariantContext truth, final List<VariantContext> calls,
                                                 final TargetCollection<Target> targets) {
        final Genotype truthGenotype = truth.getGenotype(sample);
        // if there is no truth genotype for that sample, we output the "empty" genotype.
        if (truthGenotype == null) {
            return GenotypeBuilder.create(sample, Collections.emptyList());
        }
        final int truthCopyNumber = GATKProtectedVariantContextUtils.getAttributeAsInt(truthGenotype,
                ConvertGSVariantsToSegments.GS_COPY_NUMBER_FORMAT, truthNeutralCopyNumber);
        final CNVAllele truthAllele = copyNumberToTrueAllele(truthCopyNumber);
        final List<Pair<VariantContext, Genotype>> sampleUnfilteredCalls = calls.stream()
                .map(vc -> new ImmutablePair<>(vc, vc.getGenotype(sample)))
                .filter(pair -> pair.getRight() != null)
                .filter(pair -> GATKProtectedVariantContextUtils.getAttributeAsString(pair.getRight(), GenotypeCopyNumberTriStateSegments.DISCOVERY_KEY,
                        GenotypeCopyNumberTriStateSegments.DISCOVERY_FALSE).equals(GenotypeCopyNumberTriStateSegments.DISCOVERY_TRUE))
                .collect(Collectors.toList());

        final List<Pair<VariantContext, Genotype>> sampleQualifyingCalls = sampleUnfilteredCalls.stream()
                // Filter vc/gt with no concrete call or that the call is the reference.
                .filter(pair -> !pair.getRight().getAlleles().isEmpty()
                        && !pair.getRight().getAlleles().get(0).equals(Allele.NO_CALL)
                        && !pair.getRight().getAlleles().get(0).equals(pair.getLeft().getReference()))
                // Filter vc/gt with low call quality (SQ).
                .filter(pair -> {
                        final double[] altSQs = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(pair.getRight(),
                                GenotypeCopyNumberTriStateSegments.SOME_QUALITY_KEY,
                                () -> new double[pair.getLeft().getAlleles().size() - 1], 0.0);
                        final Allele callAllele = pair.getLeft().getAlleles().get(0);
                        final int callSQIndex = pair.getLeft().getAlleles().indexOf(callAllele) - 1;
                        final double callSQ =  altSQs[callSQIndex];
                        return callSQ >= filterArguments.minimumCalledSegmentQuality;
                })
                // Filter vc/gt with small segments (small number of targets)
                .filter(pair -> targets.targetCount(pair.getLeft())
                        >= filterArguments.minimumCalledSegmentLength)
                // Filter vc/gt with that seem to be too common.
                // or if it is multi-allelic when applies.
                .filter(pair -> {
                        final double[] calledAF = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(pair.getLeft(),
                                VCFConstants.ALLELE_FREQUENCY_KEY, () -> new double[pair.getLeft().getAlleles().size() - 1],
                                0.0);
                        final double nonRefAF = MathUtils.sum(calledAF);
                        if (nonRefAF > filterArguments.maximumCalledEventFrequency) {
                            return false;
                        } else if (filterArguments.applyMultiAllelicCalledFilter && calledAF.length > 1 && calledAF[0] > 0 && calledAF[1] > 0) {
                            return false;
                        } else {
                            return true;
                        }
                })
                .collect(Collectors.toList());

        final Set<CNVAllele> calledAlleles = sampleQualifyingCalls.stream()
                .map(pair -> CNVAllele.valueOf(pair.getRight().getAllele(0)))
                .collect(Collectors.toSet());

        final Allele calledAllele = calledAlleles.size() == 1 ? calledAlleles.iterator().next().allele : Allele.NO_CALL;
        final GenotypeBuilder builder = new GenotypeBuilder();
        builder.alleles(Collections.singletonList(calledAllele));
        builder.attribute(VariantEvaluationContext.TRUTH_GENOTYPE_KEY, CNVAllele.ALL_ALLELES.indexOf(truthAllele.allele));
        builder.attribute(VariantEvaluationContext.CALLED_SEGMENTS_COUNT_KEY, sampleQualifyingCalls.size());
        // When there is more than one qualified type of event we indicate how many.
        if (calledAlleles.size() > 1) {
            builder.attribute(VariantEvaluationContext.CALLED_ALLELE_COUNTS_KEY,
                    CNVAllele.ALL_ALLELES.stream()
                        .mapToInt(allele -> (int) sampleQualifyingCalls.stream()
                                .filter(pair -> pair.getRight().getAllele(0).equals(allele, true))
                                .count())
                        .toArray());
        }
        builder.attribute(VariantEvaluationContext.CALLED_SEGMENTS_LENGTH_KEY,
                sampleQualifyingCalls.stream().mapToInt(pair -> targets.targetCount(pair.getLeft()))
                .sum());
        builder.attribute(VariantEvaluationContext.CALL_QUALITY_KEY,
                sampleQualifyingCalls.stream().mapToDouble(pair ->
                    GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(pair.getRight(), GenotypeCopyNumberTriStateSegments.SOME_QUALITY_KEY,
                            () -> new double[pair.getRight().getAlleles().size()], 0.0)[
                    pair.getLeft().getAlleles().indexOf(pair.getRight().getAllele(0))]).max());

        builder.attribute(VariantEvaluationContext.TRUTH_SEGMENT_LENGTH_KEY,
                targets.targetCount(truth));
        builder.attribute(VariantEvaluationContext.TRUTH_COPY_FRACTION_KEY,
                truthGenotype.getExtendedAttribute(ConvertGSVariantsToSegments.GS_COPY_NUMBER_FRACTION));

        final double[] truthPosteriors = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(truthGenotype, ConvertGSVariantsToSegments.GS_COPY_NUMBER_POSTERIOR, () -> new double[truthCopyNumber + 1], Double.NEGATIVE_INFINITY);
        final double truthPosteriorsSum = MathUtils.log10SumLog10(truthPosteriors);
        final double truthLog10Quality = MathUtils.log10SumLog10(truthPosteriors, 0, truthCopyNumber) +
                MathUtils.log10SumLog10(truthPosteriors, truthCopyNumber + 1, truthPosteriors.length) - truthPosteriorsSum;
        final double truthQuality = -10.0 * truthLog10Quality;
        builder.attribute(VariantEvaluationContext.TRUTH_QUALITY_KEY, truthQuality);
        if (truthQuality < filterArguments.minimumTruthSegmentQuality) {
            builder.filter(EvaluationSegmentFilter.LowQuality.acronym);
        } else {
            builder.filter(EvaluationSegmentFilter.PASS_ACRONYM);
        }

        final EvaluationClass evaluationClass;
        if (calledAllele.equals(Allele.NO_CALL)) {
            if (calledAlleles.size() > 1 && calledAlleles.contains(truthAllele)) {
                evaluationClass = EvaluationClass.MIXED_POSITIVE;
            } else if (calledAlleles.isEmpty()) {
                evaluationClass = EvaluationClass.FALSE_NEGATIVE;
            } else if (calledAlleles.contains(truthAllele)) {
                evaluationClass = EvaluationClass.TRUE_POSITIVE;
            } else  {
                evaluationClass = EvaluationClass.DISCORDANT_POSITIVE;
            }
            builder.attribute(VariantEvaluationContext.EVALUATION_CLASS_KEY, evaluationClass.acronym);
        }

        return builder.make();
    }

    private CNVAllele copyNumberToTrueAllele(final int cn) {
        if (cn == truthNeutralCopyNumber) {
            return CNVAllele.REF;
        } else if (cn < truthNeutralCopyNumber) {
            return CNVAllele.DEL;
        } else {
            return CNVAllele.DUP;
        }
    }

    /**
     * Collects all the variant in a coordinate in a list strictly sorted by the start position and the by the end
     *   position.
     *
     *
     * @param reader the variant source reader.
     * @param interval the query interval.
     * @return never {@code null}, potentially immutably
     */
    private List<VariantContext> variantQueryToList(final VCFFileReader reader, final Locatable interval) {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(reader.query(interval.getContig(),
                interval.getStart(), interval.getEnd()), Spliterator.NONNULL),false)
                .sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                .collect(Collectors.toList());
    }

    private List<SimpleInterval> composeListOfProcessingIntervalsFromInputs(final VCFFileReader truthReader, final VCFFileReader callsReader) {
        final Set<SimpleInterval> resultSet = new HashSet<>();
        for (final VariantContext vc : truthReader) {
            resultSet.add(new SimpleInterval(vc));
        }
        for (final VariantContext vc : callsReader) {
            resultSet.add(new SimpleInterval(vc));
        }
        if (resultSet.isEmpty()) {
            return Collections.emptyList();
        }
        final List<SimpleInterval> individualIntervals = new ArrayList<>(resultSet);
        Collections.sort(individualIntervals, IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR);
        final List<SimpleInterval> result = new ArrayList<>(individualIntervals.size());
        // The buffer will contain intervals that may still overlap with
        // further intervals in individual-intervals.
        final LinkedList<SimpleInterval> overlappingBuffer = new LinkedList<>();
        overlappingBuffer.add(individualIntervals.get(0));
        for (final SimpleInterval interval : individualIntervals) {
            overlappingBuffer.add(interval);
            final boolean overlaps = overlappingBuffer.getFirst().overlaps(overlappingBuffer.getLast());
            if (overlaps) {
                final SimpleInterval merged = IntervalUtils.getSpanningInterval(overlappingBuffer);
                overlappingBuffer.remove();
                overlappingBuffer.add(merged);
            } else {
                result.add(overlappingBuffer.pop());
                overlappingBuffer.add(interval);
            }
        }
        // Add the contents of the buffer.
        result.addAll(overlappingBuffer);
        return result;
    }

    private VCFFileReader openVCFReader(final File file) {
        try {
            return new VCFFileReader(file);
        } catch (final Exception ex) {
            throw new UserException.CouldNotReadInputFile(file, ex.getMessage(), ex);
        }
    }

    private void writeOverallSummaryRecordIfApplies(final EvaluationSampleSummaryWriter sampleSummaryOutputWriter) {
        if (includeOverallSummaryRecord) {
            if (sampleSummaryOutputWriter == null) {
                logger.warn("The overall sample summary record has been requested but the sample summary output file was not provided");
            } else {
                try {
                    sampleSummaryOutputWriter.writeOverallRecord();
                } catch (final IOException ex) {
                    throw new UserException.CouldNotCreateOutputFile(sampleSummaryOutputFile, ex);
                }
            }
        }
    }

    private EvaluationSiteRecordWriter createSegmentDetailOutputWriter() {
        if (segmentDetailOutputFile == null) {
            return null;
        } else {
            try {
                return new EvaluationSiteRecordWriter(segmentDetailOutputFile);
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(sampleSummaryOutputFile, ex);
            }
        }
    }


    private Set<String> composeSetOfSamplesToEvaluate(final VCFFileReader calls) {
        final List<String> sampleList = composeSampleFileList();
        if (samples.isEmpty() && sampleList.isEmpty()) {
            return composeSetOfSamplesFromInputs(calls);
        } else {
            final Set<String> result = new LinkedHashSet<>(sampleList.size() + samples.size());
            result.addAll(samples);
            result.addAll(sampleList);
            return result;
        }
    }

    private List<String> composeSampleFileList() {
        if (sampleListFile == null) {
            return Collections.emptyList();
        } else {
            try (final BufferedReader reader = new BufferedReader(new FileReader(sampleListFile))) {
                return reader.lines().collect(Collectors.toList());
            } catch (final IOException ex) {
                throw new UserException.CouldNotReadInputFile(sampleListFile, ex);
            }
        }
    }

    private Set<String> composeSetOfSamplesFromInputs(final VCFFileReader calls) {
        final LinkedHashSet<String> result = new LinkedHashSet<>();
        result.addAll(extractSampleNamesFromVCFFile(calls));
        return result;
    }

    private List<String> extractSampleNamesFromVCFFile(final VCFFileReader reader) {
        return reader.getFileHeader().getSampleNamesInOrder();
    }

    private static final class EvaluationSampleSummaryRecord {
        private final String sample;
        private Object2IntMap<EvaluationClass> countsByClass;
        private int total;

        private EvaluationSampleSummaryRecord(final String sampleName) {
            this.sample = sampleName;
            this.countsByClass = new Object2IntLinkedOpenHashMap<>(EvaluationClass.values().length);
            for (final EvaluationClass ec : EvaluationClass.values()) {
                countsByClass.put(ec, 0);
            }
            total = 0;
        }

        @Override
        public String toString() {

            final String countsString = countsByClass.entrySet().stream()
                    .map(entry -> String.format("%s = %d", entry.getKey(), entry.getValue()))
                    .collect(Collectors.joining(", "));
            return String.format("Overall stats: ALL = %d, %s", this.total, countsString);
        }

        public void increase(final EvaluationClass evalClass) {
            countsByClass.put(evalClass, countsByClass.getInt(evalClass) + 1);
        }

        public void add(final EvaluationSampleSummaryRecord record) {
            Utils.nonNull(record);
            for (final EvaluationClass evalClass : countsByClass.keySet()) {
                countsByClass.put(evalClass,
                        countsByClass.get(evalClass) + record.countsByClass.get(evalClass));
            }
        }
    }

    public static final class EvaluationSampleSummaryWriter extends TableWriter<EvaluationSampleSummaryRecord> {

        private static final String SAMPLE_NAME_COLUMN_NAME = "SAMPLE";
        private static final String TOTAL_COLUMN_NAME = "TOTAL";

        private final EvaluationSampleSummaryRecord overallRecord;

        private static final TableColumnCollection COLUMNS;

        static {
            final List<String> columnNames = new ArrayList<>(EvaluationClass.values().length + 2);
            columnNames.add(SAMPLE_NAME_COLUMN_NAME);
            for (final EvaluationClass evalClass : EvaluationClass.values()) {
                columnNames.add(evalClass.acronym);
            }
            columnNames.add(TOTAL_COLUMN_NAME);
            COLUMNS = new TableColumnCollection(columnNames);
        }

        public EvaluationSampleSummaryWriter(final File file, final String overallSample) throws IOException {
            super(file, COLUMNS);
            this.overallRecord = new EvaluationSampleSummaryRecord(overallSample);
        }

        @Override
        public void writeRecord(final EvaluationSampleSummaryRecord record) throws IOException {
            overallRecord.add(record);
            super.writeRecord(record);
        }

        @Override
        protected void composeLine(final EvaluationSampleSummaryRecord record, final DataLine dataLine) {
            dataLine.append(record.sample);
            for (final EvaluationClass evalClass : EvaluationClass.values()) {
                dataLine.append(record.countsByClass.getInt(evalClass));
            }
        }

        protected void writeOverallRecord() throws IOException {
            writeRecord(overallRecord);
        }
    }

}

