package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Interval;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Simple data structure to pass and read/write List of AllelicCounts.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicCountCollection {
    private final List<AllelicCount> counts;

    public AllelicCountCollection() {
        this.counts = new ArrayList<>();
    }

    /**
     * Constructor that reads (sequence, position, reference count, alternate count) from the specified file.
     * @param inputFile     file to read from
     */
    public AllelicCountCollection(final File inputFile) {
        try (final TableReader<AllelicCount> reader = TableUtils.reader(inputFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.containsAll(AllelicCountTableColumns.COLUMN_NAME_ARRAY))
                        throw formatExceptionFactory.apply("Missing required columns in AllelicCountCollection input file.");

                    // return the lambda to translate dataLines into AllelicCounts.
                    return (dataLine) -> {
                        final int position = dataLine.getInt(AllelicCountTableColumns.POSITION.toString());
                        final SimpleInterval interval = new SimpleInterval(dataLine.get(AllelicCountTableColumns.CONTIG.toString()), position, position);
                        final int refReadCount = dataLine.getInt(AllelicCountTableColumns.REF_COUNT.toString());
                        final int altReadCount = dataLine.getInt(AllelicCountTableColumns.ALT_COUNT.toString());
                        return new AllelicCount(interval, refReadCount, altReadCount);
                    };
                })) {
            this.counts = reader.stream().collect(Collectors.toList());
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    /**
     * Adds (interval, reference count, alternate count) to respective lists.
     * @param interval      site in 1-based format
     * @param refReadCount  number of reads at site matching the reference
     * @param altReadCount  number of reads at site different from the reference
     */
    public void add(final SimpleInterval interval, final int refReadCount, final int altReadCount) {
        counts.add(new AllelicCount(interval, refReadCount, altReadCount));
    }

    public void add(final SimpleInterval interval, final int refReadCount, final int altReadCount,
                    final Map<Character, Integer> baseCounts) {
        counts.add(new AllelicCount(interval, refReadCount, altReadCount, baseCounts));
    }

    /**
     * constructs {@link AllelicCount} from interval, ref/alt base counts, and het/hom likelihoods
     * @param interval site in 1-based format
     * @param refReadCount ref allele count
     * @param altReadCount alt allele count
     * @param baseCounts base counts map
     * @param hetLogLikelihood het log likelihood
     * @param homLogLikelihood hom log likelihood
     */
    public void add(final SimpleInterval interval, final int refReadCount, final int altReadCount,
                    final Map<Character, Integer> baseCounts, final double hetLogLikelihood, final double homLogLikelihood) {
        counts.add(new AllelicCount(interval, refReadCount, altReadCount, baseCounts, hetLogLikelihood, homLogLikelihood));
    }

    /**
     * constructs {@link AllelicCount} from interval, ref/alt base counts, and het/hom likelihoods
     * @param interval site in 1-based format
     * @param refReadCount ref allele count
     * @param altReadCount alt allele count
     * @param hetLogLikelihood het log likelihood
     * @param homLogLikelihood hom log likelihood
     */
    public void add(final SimpleInterval interval, final int refReadCount, final int altReadCount,
                    final double hetLogLikelihood, final double homLogLikelihood) {
        counts.add(new AllelicCount(interval, refReadCount, altReadCount, hetLogLikelihood, homLogLikelihood));
    }

    /** Returns an unmodifiable view of the list of AllelicCounts.   */
    public List<AllelicCount> getCounts() {
        return Collections.unmodifiableList(counts);
    }

    /**
     * Writes out (sequence, position, reference count, alternate count) to specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void write(final File outputFile) {
        try (final TableWriter<AllelicCount> writer = TableUtils.writer(outputFile,
                new TableColumnCollection(AllelicCountTableColumns.COLUMN_NAME_ARRAY),
                //lambda for filling an initially empty DataLine
                (count, dataLine) -> {
                    final Interval interval = count.getInterval();
                    final int refReadCount = count.getRefReadCount();
                    final int altReadCount = count.getAltReadCount();
                    dataLine.append(interval.getContig()).append(interval.getEnd(), refReadCount, altReadCount);
                })) {
            writer.writeAllRecords(counts);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    /**
     * Writes out (sequence, position, reference count, alternate count, A count, C count, G count, T_count)
     * to specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void writeWithBaseCounts(final File outputFile) {
        try (final TableWriter<AllelicCount> writer = TableUtils.writer(outputFile,
                new TableColumnCollection(AllelicCountTableColumns.COLUMN_NAME_ARRAY_WITH_BASE_COUNTS),
                //lambda for filling an initially empty DataLine
                (count, dataLine) -> {
                    final Interval interval = count.getInterval();
                    final int refReadCount = count.getRefReadCount();
                    final int altReadCount = count.getAltReadCount();
                    final int aCount = count.getBaseCount('A');
                    final int cCount = count.getBaseCount('C');
                    final int gCount = count.getBaseCount('G');
                    final int tCount = count.getBaseCount('T');
                    dataLine.append(interval.getContig())
                            .append(interval.getEnd(), refReadCount, altReadCount, aCount, cCount, gCount, tCount);
                })) {
            writer.writeAllRecords(counts);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    /**
     * Writes out (sequence, position, reference count, alternate count, A count, C count, G count, T_count,
     * het log likelihood, hom log likelihood) to specified file.
     * @param outputFile    file to write to (if it exists, it will be overwritten)
     */
    public void writeWithBaseCountsAndLikelihoods(final File outputFile) {
        try (final TableWriter<AllelicCount> writer = TableUtils.writer(outputFile,
                new TableColumnCollection(AllelicCountTableColumns.COLUMN_NAME_ARRAY_WITH_BASE_COUNTS_AND_LOG_LIKELIHOODS),
                //lambda for filling an initially empty DataLine
                (count, dataLine) -> {
                    final Interval interval = count.getInterval();
                    final int refReadCount = count.getRefReadCount();
                    final int altReadCount = count.getAltReadCount();
                    final int aCount = count.getBaseCount('A');
                    final int cCount = count.getBaseCount('C');
                    final int gCount = count.getBaseCount('G');
                    final int tCount = count.getBaseCount('T');
                    final double hetLogLikelihood = count.getHetLoglikelihood();
                    final double homLogLikelihood = count.getHomLoglikelihood();
                    dataLine.append(interval.getContig())
                            .append(interval.getEnd(), refReadCount, altReadCount, aCount, cCount, gCount, tCount)
                            .append(hetLogLikelihood, homLogLikelihood);
                })) {
            writer.writeAllRecords(counts);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof AllelicCountCollection)) {
            return false;
        }

        final AllelicCountCollection allelicCountCollection = (AllelicCountCollection) o;
        return counts.equals(allelicCountCollection.counts);
    }

    @Override
    public int hashCode() {
        return counts.hashCode();
    }
}
