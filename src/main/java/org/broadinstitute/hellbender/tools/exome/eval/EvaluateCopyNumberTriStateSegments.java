package org.broadinstitute.hellbender.tools.exome.eval;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

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
public final class EvaluateCopyNumberTriStateSegments extends CommandLineProgram {

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
    public static final String FREQUENCY_FILE_SHORT_NAME = "frequencies";
    public static final String FREQUENCY_FILE_FULL_NAME = "frequenciesFile";
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
            fullName = CALLS_FILE_FULL_NAME,
            optional = false)
    protected File callsFile;

    @Argument(
            doc = "File containing true event frequency",
            shortName = FREQUENCY_FILE_SHORT_NAME,
            fullName = FREQUENCY_FILE_FULL_NAME,
            optional = false
    )
    protected File frequenciesFile;

    @Argument(
            doc = "File containing true events",
            shortName = TRUTH_FILE_SHORT_NAME,
            fullName = TRUTH_FILE_FULL_NAME,
            optional = true)
    protected File truthFile;

    @Argument(
            doc = "Samples to evaluate",
            shortName = SAMPLES_SHORT_NAME,
            fullName = SAMPLES_FULL_NAME,
            optional = true
    )
    protected List<String> samples = new ArrayList<>();

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
            doc = "Frequency smoothing constant",
            shortName = FREQUENCY_SMOOTHING_SHORT_NAME,
            fullName = FREQUENCY_SMOOTHING_FULL_NAME,
            optional = true
    )
    protected double frequencySmoothing = DEFAULT_FREQUENCY_SMOOTHING;

    @Override
    protected Object doWork() {
        final Set<String> samples = composeSetOfSamplesToEvaluate();
        final TargetCollection<Target> targets = targetArguments.readTargetCollection(false);
        final IntervalsSkipList<CopyNumberTriStateFrequencies> frequencies = readFrequencies();

        SummaryOutputWriter sampleSummaryOutputWriter = null;
        EvaluationSiteRecordWriter segmentDetailOutputWriter = null;
        try {
            sampleSummaryOutputWriter = createSampleSummaryOutputWriter();
            segmentDetailOutputWriter = createSegmentDetailOutputWriter();
            for (final String sample : samples) {
                doWorkPerSample(sample, targets, sampleSummaryOutputWriter, segmentDetailOutputWriter, frequencies);
            }
            writeOverallSummaryRecordIfApplies(sampleSummaryOutputWriter);
            closeOutputWriters(sampleSummaryOutputWriter, segmentDetailOutputWriter, true);
        } catch (final RuntimeException ex) {
            closeOutputWriters(sampleSummaryOutputWriter, segmentDetailOutputWriter, false);
            throw ex;
        }
        return "SUCCESS";
    }

    private IntervalsSkipList<CopyNumberTriStateFrequencies> readFrequencies() {
        if (frequenciesFile == null) {
            return new IntervalsSkipList<>(Collections.emptyList());
        } else {
            try (final CopyNumberTriStateFrequenciesReader reader = new CopyNumberTriStateFrequenciesReader(frequenciesFile)) {
                return new IntervalsSkipList<>(reader.stream().collect(Collectors.toList()));
            } catch (final IOException ex) {
                throw new UserException.CouldNotReadInputFile(frequenciesFile);
            }
        }
    }

    private void writeOverallSummaryRecordIfApplies(final SummaryOutputWriter sampleSummaryOutputWriter) {
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

    private void closeOutputWriters(final SummaryOutputWriter sampleSummaryOutputWriter,
                                    final EvaluationSiteRecordWriter segmentDetailOutputWriter, final boolean throwAnyException) {

        final List<UserException> exceptions = new ArrayList<>(2);
        closeOutputWriter(sampleSummaryOutputWriter, sampleSummaryOutputFile, exceptions);
        closeOutputWriter(segmentDetailOutputWriter, segmentDetailOutputFile, exceptions);
        if (!exceptions.isEmpty()) {
            exceptions.forEach(ex -> logger.error(ex.getMessage()));
            if (throwAnyException) {
                throw exceptions.get(0);
            }
        }
    }

    private void closeOutputWriter(final AutoCloseable writer, final File file, final List<UserException> exceptions) {
        try {
            if (writer != null) {
                writer.close();
            }
        } catch (final Exception ex) {
            exceptions.add(new UserException.CouldNotCreateOutputFile(file, ex));
        }
    }

    private SummaryOutputWriter createSampleSummaryOutputWriter() {
        if (sampleSummaryOutputFile == null) {
            return null;
        } else {
            try {
                return new SummaryOutputWriter(sampleSummaryOutputFile, overallSampleSummaryRecordSampleName);
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(sampleSummaryOutputFile, ex);
            }
        }
    }

    private void doWorkPerSample(final String sample,
                                 final TargetCollection<Target> targets,
                                 final SummaryOutputWriter summaryOutputWriter,
                                 final EvaluationSiteRecordWriter segmentDetailOutputWriter,
                                 final IntervalsSkipList<CopyNumberTriStateFrequencies> frequencies) {
        final List<CopyNumberTriStateSegment> truthSegments =
                readRelevantSegments(sample, truthFile, targets);
        final List<CopyNumberTriStateSegment> calledSegments =
                readRelevantSegments(sample, callsFile, targets);

        final IntervalsSkipList<CopyNumberTriStateSegment> indexedCalledSegments =
                new IntervalsSkipList<>(calledSegments);

        final List<EvaluationSiteRecord> records = new ArrayList<>(truthSegments.size());

        final SampleSummaryOutputRecord sampleSummaryOutputRecord = new SampleSummaryOutputRecord(sample);
        // Iterate over truth calls and do most of the counting.
        for (final CopyNumberTriStateSegment truthSegment : truthSegments) {
            final int targetCount = targets.targetCount(truthSegment);
            final List<CopyNumberTriStateSegment> overlappingCalls =
                    indexedCalledSegments.getOverlapping(truthSegment.getInterval());
            final CopyNumberTriStateFrequencies frequency = determineBestFrequencyRecordForInterval(truthSegment.getInterval(), frequencies);
            if (overlappingCalls.isEmpty()) {
                sampleSummaryOutputRecord.falseNegative++;
                records.add(EvaluationSiteRecord.falseNegative(sample, targetCount, truthSegment, frequency));
            } else if (isMixedCall(overlappingCalls)) {
                sampleSummaryOutputRecord.mixedPositive++;
                records.add(EvaluationSiteRecord.otherPositive(sample, targetCount, EvaluationClass.MIXED_POSITIVE, truthSegment, overlappingCalls, frequency));
            } else if (overlappingCalls.get(0).getCall() == truthSegment.getCall()) {
                sampleSummaryOutputRecord.truePositive++;
                records.add(EvaluationSiteRecord.otherPositive(sample, targetCount, EvaluationClass.TRUE_POSITIVE, truthSegment, overlappingCalls, frequency));
            } else {
                records.add(EvaluationSiteRecord.otherPositive(sample, targetCount, EvaluationClass.DISCORDANT_POSITIVE, truthSegment, overlappingCalls, frequency));
                sampleSummaryOutputRecord.discordantPositive++;
            }
        }

        // then we need to look into the calls not overlapping with any truth.
        final IntervalsSkipList<CopyNumberTriStateSegment> indexedTruthSegments =
                new IntervalsSkipList<>(truthSegments);

        for (final CopyNumberTriStateSegment callSegment : calledSegments) {
            final List<CopyNumberTriStateSegment> overlappingTruth =
                    indexedTruthSegments.getOverlapping(callSegment.getInterval());
            if (overlappingTruth.isEmpty()) {
                records.add(EvaluationSiteRecord.unknownPositive(sample, targets.targetCount(callSegment), callSegment, null));
                sampleSummaryOutputRecord.unknownPositive++;
            }
        }

        records.stream()
                .sorted(Comparator.comparing(EvaluationSiteRecord::getInterval,
                        IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR))
                .forEach(r -> writeSegmentDetailOutput(segmentDetailOutputWriter, filterArguments.applyFiltersOn(r)));
        writeSampleSummaryRecord(summaryOutputWriter, sampleSummaryOutputRecord);
    }

    /**
     * From all the overlapping frequencies object gets the one that best overlaps the given interval.
     *
     * <p>
     *     The "best" overlapping frequency object is the one where the ratio of the overlapping region length and
     *     the union region length is the largest. In case of a tie, is then the one with largest non-neutral copy number
     *     frequencies ({@code 1.0 - neutral copy frequency}).
     * </p>
     *
     * @param interval the target interval.
     * @param frequencies the potential overlapping frequency objects.
     * @return {@code null} if there is no overlapping frequencies objects.
     */
    private CopyNumberTriStateFrequencies determineBestFrequencyRecordForInterval(final SimpleInterval interval,
                                                                                        final IntervalsSkipList<CopyNumberTriStateFrequencies> frequencies) {
        final List<CopyNumberTriStateFrequencies> overlapping = frequencies.getOverlapping(interval);
        if (overlapping.isEmpty()) {
            return null;
        } else if (overlapping.size() == 1) {  // quick short cut for the trivial case of one overlapping frequency object.
            return overlapping.get(0);
        } else {

            final Comparator<CopyNumberTriStateFrequencies> bestOverlapping = EvaluationUtils.comparingOverlapToUnionRatio(interval);
            final Comparator<CopyNumberTriStateFrequencies> bestMatching = bestOverlapping
                .thenComparing(a -> a.frequencyFor(CopyNumberTriState.NEUTRAL, frequencySmoothing));
            return overlapping.stream().sorted(bestMatching).findFirst().get();
        }
    }

    private void writeSegmentDetailOutput(final EvaluationSiteRecordWriter segmentDetailOutputWriter, final EvaluationSiteRecord evaluationSiteRecord) {
        if (segmentDetailOutputWriter != null) {
            try {
                segmentDetailOutputWriter.writeRecord(evaluationSiteRecord);
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(segmentDetailOutputFile, ex);
            }
        }
    }

    private void writeSampleSummaryRecord(final SummaryOutputWriter summaryOutputWriter,
                                          final SampleSummaryOutputRecord sampleSummaryOutputRecord) {
        if (summaryOutputWriter != null) {
            try {
                summaryOutputWriter.writeRecord(sampleSummaryOutputRecord);
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(sampleSummaryOutputFile, ex);
            }
        }
    }

    private boolean isMixedCall(final List<CopyNumberTriStateSegment> overlappingCalls) {
        if (overlappingCalls.size() <= 1) {
            return false;
        } else {
            final CopyNumberTriState candidateCall = overlappingCalls.get(0).getCall();
            return overlappingCalls.stream().anyMatch(s -> s.getCall() != candidateCall);
        }
    }

    private List<CopyNumberTriStateSegment> readRelevantSegments(final String sample, final File file, final TargetCollection<Target> targets) {
        if (truthFile == null) {
            return Collections.emptyList();
        } else {
            try (final CopyNumberTriStateSegmentRecordReader reader = new CopyNumberTriStateSegmentRecordReader(file)) {
                return reader.stream()
                        .filter(record -> record.getSampleName().equals(sample))
                        .map(CopyNumberTriStateSegmentRecord::getSegment)
                        .filter(segment -> segment.getCall() != CopyNumberTriState.NEUTRAL)
                        .filter(segment -> targets.indexRange(segment).size() > 0)
                        .collect(Collectors.toList());
            } catch (final IOException ex) {
                throw new UserException.CouldNotReadInputFile(file, ex);
            }
        }
    }

    private Set<String> composeSetOfSamplesToEvaluate() {
        final List<String> sampleList = composeSampleFileList();
        if (samples.isEmpty() && sampleList.isEmpty()) {
            return composeSetOfSamplesFromInputs();
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

    private Set<String> composeSetOfSamplesFromInputs() {
        final LinkedHashSet<String> result = new LinkedHashSet<>();
        result.addAll(extractSampleNamesFromSegmentsFile(callsFile));
        result.addAll(extractSampleNamesFromSegmentsFile(truthFile));
        return result;
    }

    private List<String> extractSampleNamesFromSegmentsFile(final File file) {
        if (file == null) {
            return Collections.emptyList();
        } else {
            try (final CopyNumberTriStateSegmentRecordReader reader =
                         new CopyNumberTriStateSegmentRecordReader(file)) {
                return reader.stream()
                        .map(CopyNumberTriStateSegmentRecord::getSampleName)
                        .collect(Collectors.toList());
            } catch (final IOException ex) {
                throw new UserException.CouldNotReadInputFile(file, ex);
            }
        }
    }

    private static final class SampleSummaryOutputRecord {
        private final String sample;
        private int truePositive = 0;
        private int discordantPositive = 0;
        private int mixedPositive = 0;
        private int unknownPositive = 0;
        private int falseNegative = 0;

        private SampleSummaryOutputRecord(final String sampleName) {
            this.sample = sampleName;
        }

        public int calculateTotal() {
            return truePositive + discordantPositive + mixedPositive
                    + unknownPositive + falseNegative;
        }

        @Override
        public String toString() {
            return String.format("Overall stats: ALL = %d, TP = %d, FN = %d, UP = %d, MP = %d, DP = %d",
                    calculateTotal(), truePositive, falseNegative, unknownPositive, mixedPositive, discordantPositive);
        }
    }

    public static final class SummaryOutputWriter extends TableWriter<SampleSummaryOutputRecord> {

        private static final String SAMPLE_NAME_COLUMN_NAME = "SAMPLE";
        private static final String TRUE_POSITIVE_COLUMN_NAME = EvaluationClass.TRUE_POSITIVE.acronym;
        private static final String FALSE_NEGATIVE_COLUMN_NAME = EvaluationClass.FALSE_NEGATIVE.acronym;
        private static final String UNKNOWN_POSITIVE_COLUMN_NAME = EvaluationClass.UNKNOWN_POSITIVE.acronym;
        private static final String DISCORDANT_POSITIVE_COLUMN_NAME = EvaluationClass.DISCORDANT_POSITIVE.acronym;
        private static final String MIXED_POSITIVE_COLUMN_NAME = EvaluationClass.MIXED_POSITIVE.acronym;
        private static final String TOTAL_COLUMN_NAME = "TOTAL";

        private final SampleSummaryOutputRecord overallRecord;

        private static final TableColumnCollection COLUMNS = new TableColumnCollection(
                SAMPLE_NAME_COLUMN_NAME,
                TRUE_POSITIVE_COLUMN_NAME,
                FALSE_NEGATIVE_COLUMN_NAME,
                UNKNOWN_POSITIVE_COLUMN_NAME,
                DISCORDANT_POSITIVE_COLUMN_NAME,
                MIXED_POSITIVE_COLUMN_NAME,
                TOTAL_COLUMN_NAME
        );

        public SummaryOutputWriter(final File file, final String overallSample) throws IOException {
            super(file, COLUMNS);
            this.overallRecord = new SampleSummaryOutputRecord(overallSample);
        }

        @Override
        protected SampleSummaryOutputRecord beforeWriteRecord(final SampleSummaryOutputRecord record) {
            overallRecord.truePositive += record.truePositive;
            overallRecord.discordantPositive += record.discordantPositive;
            overallRecord.unknownPositive += record.unknownPositive;
            overallRecord.falseNegative += record.falseNegative;
            overallRecord.mixedPositive += record.mixedPositive;
            return record;
        }

        @Override
        protected void composeLine(final SampleSummaryOutputRecord record, final DataLine dataLine) {
            dataLine.append(record.sample)
                    .append(record.truePositive)
                    .append(record.falseNegative)
                    .append(record.unknownPositive)
                    .append(record.discordantPositive)
                    .append(record.mixedPositive)
                    .append(record.calculateTotal());
        }

        protected void writeOverallRecord() throws IOException {
            writeRecord(overallRecord);
        }
    }

}

