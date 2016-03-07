package org.broadinstitute.hellbender.tools.exome;

import java.util.Arrays;
import java.util.stream.Stream;

/**
 * Created by davidben on 11/30/15.
 */
public enum AllelicCountTableColumns {
    CONTIG("CONTIG"), POSITION("POS"), REF_COUNT("REF_COUNT"), ALT_COUNT("ALT_COUNT"),
    A_COUNT("A_COUNT"), C_COUNT("C_COUNT"), G_COUNT("G_COUNT"), T_COUNT("T_COUNT"),
    HET_LOG_LIKELIHOOD("HET_LOG_LIKELIHOOD"), HOM_LOG_LIKELIHOOD("HOM_LOG_LIKELIHOOD");

    private String columnName;

    AllelicCountTableColumns(String columnName) {this.columnName = columnName; }

    @Override
    public String toString() {
        return columnName;
    }

    public static final String[] COLUMN_NAME_ARRAY =
            Stream.of(Arrays.copyOfRange(values(), 0, 4)).map(AllelicCountTableColumns::toString).toArray(String[]::new);

    public static final String[] COLUMN_NAME_ARRAY_WITH_BASE_COUNTS =
            Stream.of(Arrays.copyOfRange(values(), 0, 8)).map(AllelicCountTableColumns::toString).toArray(String[]::new);

    public static final String[] COLUMN_NAME_ARRAY_WITH_BASE_COUNTS_AND_LOG_LIKELIHOODS =
            Stream.of(Arrays.copyOfRange(values(), 0, 10)).map(AllelicCountTableColumns::toString).toArray(String[]::new);

}