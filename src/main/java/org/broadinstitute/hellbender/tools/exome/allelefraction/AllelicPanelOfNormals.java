package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

/**
 * Represents the panel of normals used for allele-bias correction.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AllelicPanelOfNormals {
    private final Map<SimpleInterval, ReadCountPair> siteToCountPairMap = new HashMap<>();

    public AllelicPanelOfNormals(final File inputFile) {
        Utils.nonNull(inputFile);
        Utils.regularReadableUserFile(inputFile);

        final AllelicCountCollection counts = new AllelicCountCollection(inputFile);
        for (final AllelicCount count : counts.getCounts()) {
            final SimpleInterval site = count.getInterval();
            final ReadCountPair countPair = new ReadCountPair(count.getRefReadCount(), count.getAltReadCount());
            if (siteToCountPairMap.containsKey(site)) {
                throw new UserException.BadInput("Input file for allelic panel of normals contains duplicate sites.");
            } else {
                siteToCountPairMap.put(site, countPair);
            }
        }
    }

    public int getRefReadCount(final SimpleInterval site) {
        Utils.nonNull(site);

        if (siteToCountPairMap.containsKey(site)) {
            return siteToCountPairMap.get(site).getRefReadCount();
        } else {
            throw new IllegalArgumentException("Site is not present in the allelic panel of normals.");
        }
    }

    public int getAltReadCount(final SimpleInterval site) {
        Utils.nonNull(site);

        if (siteToCountPairMap.containsKey(site)) {
            return siteToCountPairMap.get(site).getAltReadCount();
        } else {
            throw new IllegalArgumentException("Site is not present in the allelic panel of normals.");
        }
    }

    private class ReadCountPair {
        private final Pair<Integer, Integer> readCountPair;

        private ReadCountPair(final int refReadCount, final int altReadCount) {
            ParamUtils.isPositiveOrZero(refReadCount, "Can't construct ReadCountPair with negative read counts.");
            ParamUtils.isPositiveOrZero(altReadCount, "Can't construct ReadCountPair with negative read counts.");
            readCountPair = Pair.of(refReadCount, altReadCount);
        }

        private int getRefReadCount() { return readCountPair.getLeft();  }
        private int getAltReadCount() { return readCountPair.getRight();  }
    }
}
