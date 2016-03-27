package org.broadinstitute.hellbender.tools.exome.eval;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.tools.exome.CopyNumberTriStateSegment;
import org.broadinstitute.hellbender.utils.hmm.CopyNumberTriState;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * TODO document this.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
class EvaluationSiteRecordWriter extends TableWriter<EvaluationSiteRecord> {

    private static final String SAMPLE_NAME_COLUMN = "SAMPLE";
    private static final String EVALUATION_CLASS_COLUMN = "CLASS";
    private static final String CONTIG_COLUMN = "CONTIG";
    private static final String START_COLUMN = "START";
    private static final String END_COLUMN = "END";
    private static final String TARGET_COUNT_COLUMN = "NTARGETS";
    private static final String TRUTH_COLUMN = "TRUTH";
    private static final String CALL_COLUMN = "CALL";
    private static final String DELETION_FREQUENCY = "DELETION_FREQUENCY";
    private static final String DUPLICATION_FREQUENCY = "DUPLICATION_FREQUENCY";

    private static final TableColumnCollection COLUMNS = new TableColumnCollection(
            SAMPLE_NAME_COLUMN,
            EVALUATION_CLASS_COLUMN,
            CONTIG_COLUMN,
            START_COLUMN,
            END_COLUMN,
            TARGET_COUNT_COLUMN,
            TRUTH_COLUMN,
            CALL_COLUMN,
            DELETION_FREQUENCY,
            DUPLICATION_FREQUENCY
    );

    private static final String SEGMENT_SEPARATOR = ";";
    private static final String ATTRIBUTE_SEPARATOR = ":";
    private static final String NO_SEGMENT_STRING = ".";

    public EvaluationSiteRecordWriter(final File file) throws IOException {
        super(file, COLUMNS);
        writeComment("Possible classes in the CLASS column: ");
        for (final EvaluationClass clazz : EvaluationClass.values()) {
            writeComment(String.format("    %s : %s", clazz.acronym, clazz.name().replace("_", " ")));
        }
        writeComment(String.format("Content of %s and %s columns:", TRUTH_COLUMN, CALL_COLUMN));
        writeComment("   They may contain multiple segment information or none (represented with '.')");
        writeComment(String.format("   Each segment record is separated with a '%s' and each field in each segment is separated by '%s'", SEGMENT_SEPARATOR, ATTRIBUTE_SEPARATOR));
        writeComment("   The segments attributes are in this order:");
        writeComment("       " + String.join(ATTRIBUTE_SEPARATOR, "Start", "End", "Call", "TargetCount", "Mean Cov.", "Std Cov.", "Some Quality", "Start Quality", "End Quality"));
    }

    @Override
    protected void composeLine(final EvaluationSiteRecord record,
                               final DataLine dataLine) {
        dataLine.append(record.sample)
                .append(record.evaluationClass.toString())
                .append(record.interval.getContig())
                .append(record.interval.getStart())
                .append(record.interval.getEnd())
                .append(record.targetCount)
                .append(composeSegmentsString(record.truths))
                .append(composeSegmentsString(record.calls))
                .append(record.deletionFrequency)
                .append(record.duplicationFrequency);

    }

    private String composeSegmentsString(final List<Pair<CopyNumberTriStateSegment, Set<EvaluationSegmentFilter>>> segments) {
        if (segments.isEmpty()) {
            return NO_SEGMENT_STRING;
        } else {
            return segments.stream()
                    .map(Pair::getLeft)
                    .map(s -> Stream.of(s.getStart(), s.getEnd(),
                            s.getCall().toCallString(), s.getTargetCount(), s.getMean(), s.getStdev(), s.getSomeQuality(),
                            s.getStartQuality(), s.getEndQuality()).map(Object::toString)
                            .collect(Collectors.joining(ATTRIBUTE_SEPARATOR)))
                    .collect(Collectors.joining(SEGMENT_SEPARATOR));
        }
    }
}
