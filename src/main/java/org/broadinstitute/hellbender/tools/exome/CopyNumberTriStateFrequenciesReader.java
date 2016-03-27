package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;

import java.io.File;
import java.io.IOException;

/**
 * Reads {@link CopyNumberTriStateFrequencies} table file.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class CopyNumberTriStateFrequenciesReader extends TableReader<CopyNumberTriStateFrequencies> {

    public CopyNumberTriStateFrequenciesReader(final File file) throws IOException {
        super(file);
    }

    @Override
    protected void processColumns(final TableColumnCollection columns) {
        columns.matchesExactly(CopyNumberTriStateFrequenciesWriter.COLUMNS.names().stream().toArray(String[]::new));
    }

    @Override
    protected CopyNumberTriStateFrequencies createRecord(final DataLine dataLine) {
        final SimpleInterval interval = new SimpleInterval(dataLine.get(0), dataLine.getInt(1), dataLine.getInt(2));
        final int deleteCount = dataLine.getInt(3);
        final int neutralCount = dataLine.getInt(4);
        final int duplicationCount = dataLine.getInt(5);
        return new CopyNumberTriStateFrequencies(interval, deleteCount, neutralCount, duplicationCount);
    }
}
