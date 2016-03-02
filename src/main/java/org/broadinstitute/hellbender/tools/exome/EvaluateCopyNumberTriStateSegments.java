package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
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

    public static final String DEFAULT_OVERALL_SUMMARY_SAMPLE_NAME = "ALL";

    @ArgumentCollection
    protected TargetArgumentCollection targetArguments = new TargetArgumentCollection();

    @Argument(
            doc = "File contained the called segments",
            shortName = CALLS_FILE_SHORT_NAME,
            fullName = CALLS_FILE_FULL_NAME,
            optional = false)
    protected File callsFile;

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

    @Override
    protected Object doWork() {
        final Set<String> samples = composeSetOfSamplesToEvaluate();
        final TargetCollection<Target> targets = targetArguments.readTargetCollection(false);

        SummaryOutputWriter sampleSummaryOutputWriter = null;
        SegmentDetailOutputWriter segmentDetailOutputWriter = null;
        try {
            final SampleSummaryOutputRecord overallSampleSummaryRecord =
                    new SampleSummaryOutputRecord(overallSampleSummaryRecordSampleName);
            sampleSummaryOutputWriter = createSampleSummaryOutputWriter();
            segmentDetailOutputWriter = createSegmentDetailOutputWriter();
            for (final String sample : samples) {
                doWorkPerSample(sample, targets, sampleSummaryOutputWriter, overallSampleSummaryRecord, segmentDetailOutputWriter);
            }
            writeOverallSummaryRecordIfApplies(sampleSummaryOutputWriter, overallSampleSummaryRecord);
            closeOutputWriters(sampleSummaryOutputWriter, segmentDetailOutputWriter, true);
        } catch (final RuntimeException ex) {
            closeOutputWriters(sampleSummaryOutputWriter, segmentDetailOutputWriter, false);
            throw ex;
        }
        return "SUCCESS";
    }

    private void writeOverallSummaryRecordIfApplies(final SummaryOutputWriter sampleSummaryOutputWriter, final SampleSummaryOutputRecord overallSampleSummaryRecord) {
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
        logger.info(overallSampleSummaryRecord.toString());
    }

    private SegmentDetailOutputWriter createSegmentDetailOutputWriter() {
        if (segmentDetailOutputFile == null) {
            return null;
        } else {
            try {
                return new SegmentDetailOutputWriter(segmentDetailOutputFile);
            } catch (final IOException ex) {
                throw new UserException.CouldNotCreateOutputFile(sampleSummaryOutputFile, ex);
            }
        }
    }

    private void closeOutputWriters(final SummaryOutputWriter sampleSummaryOutputWriter,
                                    final SegmentDetailOutputWriter segmentDetailOutputWriter, final boolean throwAnyException) {

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
                                 final SampleSummaryOutputRecord overallSampleSummaryRecord,
                                 final SegmentDetailOutputWriter segmentDetailOutputWriter) {
        final List<CopyNumberTriStateSegment> truthSegments =
                readRelevantSegments(sample, truthFile, targets);
        final List<CopyNumberTriStateSegment> calledSegments =
                readRelevantSegments(sample, callsFile, targets);

        final IntervalsSkipList<CopyNumberTriStateSegment> indexedCalledSegments =
                new IntervalsSkipList<>(calledSegments);

        // Iterate over truth calls and do most of the counting.
        for (final CopyNumberTriStateSegment truthSegment : truthSegments) {
            final List<CopyNumberTriStateSegment> overlappingCalls =
                    indexedCalledSegments.getOverlapping(truthSegment.getInterval());
            if (overlappingCalls.isEmpty()) {
                overallSampleSummaryRecord.falseNegative++;
                writeSegmentDetailOutput(segmentDetailOutputWriter,
                        SegmentDetailOutputRecord.falseNegative(sample, truthSegment));
            } else if (isMixedCall(overlappingCalls)) {
                overallSampleSummaryRecord.mixedPositive++;
                writeSegmentDetailOutput(segmentDetailOutputWriter,
                        SegmentDetailOutputRecord.otherPositive(sample, EvaluationClass.MIXED_POSITIVE, truthSegment, overlappingCalls));
            } else if (overlappingCalls.get(0).getCall() == truthSegment.getCall()) {
                overallSampleSummaryRecord.truePositive++;
                writeSegmentDetailOutput(segmentDetailOutputWriter,
                        SegmentDetailOutputRecord.otherPositive(sample, EvaluationClass.TRUE_POSITIVE, truthSegment, overlappingCalls));
            } else {
                writeSegmentDetailOutput(segmentDetailOutputWriter,
                        SegmentDetailOutputRecord.otherPositive(sample, EvaluationClass.DISCORDANT_POSITIVE, truthSegment, overlappingCalls));
                overallSampleSummaryRecord.discordantPositive++;
            }
        }

        // then we need to look into the calls not overlapping with any truth.
        final IntervalsSkipList<CopyNumberTriStateSegment> indexedTruthSegments =
                new IntervalsSkipList<>(truthSegments);

        for (final CopyNumberTriStateSegment callSegment : calledSegments) {
            final List<CopyNumberTriStateSegment> overlappingTruth =
                    indexedTruthSegments.getOverlapping(callSegment.getInterval());
            if (overlappingTruth.isEmpty()) {
                writeSegmentDetailOutput(segmentDetailOutputWriter,
                        SegmentDetailOutputRecord.unknownPositive(sample, callSegment));
                overallSampleSummaryRecord.unknownPositive++;
            }
        }
        writeSampleSummaryRecord(summaryOutputWriter, overallSampleSummaryRecord);
    }

    private void writeSegmentDetailOutput(final SegmentDetailOutputWriter segmentDetailOutputWriter, final SegmentDetailOutputRecord segmentDetailOutputRecord) {
        if (segmentDetailOutputWriter != null) {
            try {
                segmentDetailOutputWriter.writeRecord(segmentDetailOutputRecord);
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
            final CopyNumberTriState candidateCall = overlappingCalls.get(0).call;
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
            return super.beforeWriteRecord(record);
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

    private static class SegmentDetailOutputRecord {

        private final String sample;
        private final List<CopyNumberTriStateSegment> calls;
        private final List<CopyNumberTriStateSegment> truths;
        private final EvaluationClass evaluationClass;
        private final SimpleInterval interval;

        private SegmentDetailOutputRecord(final String sample, final EvaluationClass evaluationClass,
                                          final List<CopyNumberTriStateSegment> calls, final List<CopyNumberTriStateSegment> truths) {
            this.sample = Utils.nonNull(sample);
            this.evaluationClass = Utils.nonNull(evaluationClass);
            this.calls = Utils.nonNull(calls);
            this.truths = Utils.nonNull(truths);
            if (truths.size() == 1) {
                interval = truths.get(0).getInterval();
            } else if (calls.size() == 1) {
                interval = calls.get(0).getInterval();
            } else {
                throw new IllegalArgumentException("either the truth or the call list must have exactly one element");
            }
        }

        public SimpleInterval getInterval() {
            return interval;
        }

        private static SegmentDetailOutputRecord falseNegative(final String sample, final CopyNumberTriStateSegment truth) {
            return new SegmentDetailOutputRecord(sample, EvaluationClass.FALSE_NEGATIVE,
                    Collections.emptyList(), Collections.singletonList(Utils.nonNull(truth)));
        }

        private static SegmentDetailOutputRecord unknownPositive(final String sample, final CopyNumberTriStateSegment call) {
            return new SegmentDetailOutputRecord(sample, EvaluationClass.UNKNOWN_POSITIVE,
                    Collections.singletonList(Utils.nonNull(call)), Collections.emptyList());
        }

        private static SegmentDetailOutputRecord otherPositive(final String sample, final EvaluationClass positiveClass,
                                                        final CopyNumberTriStateSegment truth, final List<CopyNumberTriStateSegment> calls) {
            return new SegmentDetailOutputRecord(sample, positiveClass, Utils.nonNull(calls), Collections.singletonList(Utils.nonNull(truth)));
        }

        @Override
        public String toString() {
            return "";
        }
    }

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


        EvaluationClass(final String acronym) {
            this.acronym = acronym;
        }

        @Override
        public String toString() {
            return acronym;
        }
    }

    private static class SegmentDetailOutputWriter extends TableWriter<SegmentDetailOutputRecord> {

        private static final String SAMPLE_NAME_COLUMN = "SAMPLE";
        private static final String EVALUATION_CLASS_COLUMN = "CLASS";
        private static final String CONTIG_COLUMN = "CONTIG";
        private static final String START_COLUMN  = "START";
        private static final String END_COLUMN    = "END";
        private static final String TRUTH_COLUMN  = "TRUTH";
        private static final String CALL_COLUMN   = "CALL";
        private static final TableColumnCollection COLUMNS = new TableColumnCollection(
                SAMPLE_NAME_COLUMN,
                EVALUATION_CLASS_COLUMN,
                CONTIG_COLUMN,
                START_COLUMN,
                END_COLUMN,
                TRUTH_COLUMN,
                CALL_COLUMN
        );

        private static final String SEGMENT_SEPARATOR = ";";
        private static final String ATTRIBUTE_SEPARATOR = ":";
        private static final String NO_SEGMENT_STRING = ".";

        public SegmentDetailOutputWriter(final File file) throws IOException {
            super(file, COLUMNS);
        }

        @Override
        protected void composeLine(final SegmentDetailOutputRecord record,
                                   final DataLine dataLine) {
            dataLine.append(record.sample)
                    .append(record.evaluationClass.toString())
                    .append(composeSegmentsString(record.truths))
                    .append(composeSegmentsString(record.calls));

        }

        private String composeSegmentsString(final List<CopyNumberTriStateSegment> truths) {
            if (truths.isEmpty()) {
                return NO_SEGMENT_STRING;
            } else {
                return truths.stream()
                        .map(s ->  String.join(ATTRIBUTE_SEPARATOR, Arrays.asList("" + s.getStart(), "" + s.getEnd(),
                                "" + s.getCall(), "" + s.getSomeQuality(),
                                "" + s.getStartQuality(), "" + s.getEndQuality())))
                        .collect(Collectors.joining(SEGMENT_SEPARATOR));
            }
        }
    }

}

