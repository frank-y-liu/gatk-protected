package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.engine.AssemblyRegion;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.clipping.ReadClipper;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.fragments.FragmentCollection;
import org.broadinstitute.hellbender.utils.fragments.FragmentUtils;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadCoordinateComparator;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by davidben on 9/8/16.
 */
public class AssemblyBasedCallerUtils {
    /**
     * Returns a map with the original read as a key and the realigned read as the value.
     * <p>
     *     Missing keys or equivalent key and value pairs mean that the read was not realigned.
     * </p>
     * @return never {@code null}
     */
    public static Map<GATKRead, GATKRead> realignReadsToTheirBestHaplotype(final ReadLikelihoods<Haplotype> originalReadLikelihoods, final Haplotype refHaplotype, final Locatable paddedReferenceLoc) {
        final Collection<ReadLikelihoods<Haplotype>.BestAllele> bestAlleles = originalReadLikelihoods.bestAlleles();
        final Map<GATKRead, GATKRead> result = new HashMap<>(bestAlleles.size());

        for (final ReadLikelihoods<Haplotype>.BestAllele bestAllele : bestAlleles) {
            final GATKRead originalRead = bestAllele.read;
            final Haplotype bestHaplotype = bestAllele.allele;
            final boolean isInformative = bestAllele.isInformative();
            final GATKRead realignedRead = AlignmentUtils.createReadAlignedToRef(originalRead, bestHaplotype, refHaplotype, paddedReferenceLoc.getStart(), isInformative);
            result.put(originalRead, realignedRead);
        }
        return result;
    }

    public static void finalizeRegion(final AssemblyRegion region,
                                      final boolean errorCorrectReads,
                                      final boolean dontUseSoftClippedBases,
                                      final byte minTailQuality,
                                      final SAMFileHeader readsHeader,
                                      final SampleList samplesList) {
        if ( region.isFinalized() ) {
            return;
        }

        // Loop through the reads hard clipping the adaptor and low quality tails
        final List<GATKRead> readsToUse = new ArrayList<>(region.getReads().size());
        for( final GATKRead myRead : region.getReads() ) {
            final byte minTailQualityToUse = errorCorrectReads ? HaplotypeCallerEngine.MIN_TAIL_QUALITY_WITH_ERROR_CORRECTION : minTailQuality;
            GATKRead clippedRead = ReadClipper.hardClipLowQualEnds(myRead, minTailQualityToUse);

            // remove soft clips if we cannot reliably clip off adapter sequence or if the user doesn't want to use soft clips at all
            // otherwie revert soft clips so that we see the alignment start and end assuming the soft clips are all matches
            // TODO -- WARNING -- still possibility that unclipping the soft clips will introduce bases that aren't
            // TODO -- truly in the extended region, as the unclipped bases might actually include a deletion
            // TODO -- w.r.t. the reference.  What really needs to happen is that kmers that occur before the
            // TODO -- reference haplotype start must be removed
            clippedRead = dontUseSoftClippedBases || ! ReadUtils.hasWellDefinedFragmentSize(clippedRead) ?
                    ReadClipper.hardClipSoftClippedBases(clippedRead) : ReadClipper.revertSoftClippedBases(clippedRead);

            clippedRead = clippedRead.isUnmapped() ? clippedRead : ReadClipper.hardClipAdaptorSequence(clippedRead);
            if ( ! clippedRead.isEmpty() && clippedRead.getCigar().getReadLength() > 0 ) {
                clippedRead = ReadClipper.hardClipToRegion( clippedRead, region.getExtendedSpan().getStart(), region.getExtendedSpan().getEnd() );
                if ( region.readOverlapsRegion(clippedRead) && clippedRead.getLength() > 0 ) {
                    readsToUse.add(clippedRead);
                }
            }
        }

        // TODO -- Performance optimization: we partition the reads by sample 4 times right now; let's unify that code.
        // final List<GATKRead> downsampledReads = DownsamplingUtils.levelCoverageByPosition(ReadUtils.sortReadsByCoordinate(readsToUse), maxReadsInRegionPerSample, minReadsPerAlignmentStart);
        Collections.sort(readsToUse, new ReadCoordinateComparator(readsHeader)); // TODO: sort may be unnecessary here

        // handle overlapping read pairs from the same fragment
        cleanOverlappingReadPairs(readsToUse, samplesList, readsHeader);

        region.clearReads();
        region.addAll(readsToUse);
        region.setFinalized(true);
    }

    /**
     * Clean up reads/bases that overlap within read pairs
     *
     * @param reads the list of reads to consider
     */
    private static void cleanOverlappingReadPairs(final List<GATKRead> reads, final SampleList samplesList, final SAMFileHeader readsHeader) {
        for ( final List<GATKRead> perSampleReadList : splitReadsBySample(samplesList, readsHeader, reads).values() ) {
            final FragmentCollection<GATKRead> fragmentCollection = FragmentCollection.create(perSampleReadList);
            for ( final List<GATKRead> overlappingPair : fragmentCollection.getOverlappingPairs() ) {
                FragmentUtils.adjustQualsOfOverlappingPairedFragments(overlappingPair);
            }
        }
    }

    public static Map<String, List<GATKRead>> splitReadsBySample( final SampleList samplesList, final SAMFileHeader header, final Collection<GATKRead> reads ) {
        final Map<String, List<GATKRead>> returnMap = new HashMap<>();
        for (final String sample : samplesList.asListOfSamples()) {
            returnMap.put(sample, new ArrayList<>());
        }

        for ( final GATKRead read : reads ) {
            returnMap.get(ReadUtils.getSampleName(read, header)).add(read);
        }

        return returnMap;
    }

    /**
     * Helper function to create the reference haplotype out of the active region and a padded loc
     * @param region the active region from which to generate the reference haplotype
     * @param paddedReferenceLoc the interval which includes padding and shows how big the reference haplotype should be
     * @return a non-null haplotype
     */
    public static Haplotype createReferenceHaplotype(final AssemblyRegion region, final SimpleInterval paddedReferenceLoc, final IndexedFastaSequenceFile referenceReader) {
        return ReferenceConfidenceModel.createReferenceHaplotype(region, region.getAssemblyRegionReference(referenceReader), paddedReferenceLoc);
    }

    public static SimpleInterval getPaddedReferenceLoc(final AssemblyRegion region, final int referencePadding, final IndexedFastaSequenceFile referenceReader) {
        final int padLeft = Math.max(region.getExtendedSpan().getStart() - referencePadding, 1);
        final int padRight = Math.min(region.getExtendedSpan().getEnd() + referencePadding, referenceReader.getSequenceDictionary().getSequence(region.getExtendedSpan().getContig()).getSequenceLength());
        return new SimpleInterval(region.getExtendedSpan().getContig(), padLeft, padRight);
    }

    public static CachingIndexedFastaSequenceFile createReferenceReader(final String reference) {
        try {
            // fasta reference reader to supplement the edges of the reference sequence
            return new CachingIndexedFastaSequenceFile(new File(reference));
        } catch( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(new File(reference), e);
        }
    }

    /**
     * Instantiates the appropriate likelihood calculation engine.
     *
     * @return never {@code null}.
     */
    public static ReadLikelihoodCalculationEngine createLikelihoodCalculationEngine(final LikelihoodEngineArgumentCollection likelihoodArgs,
                                                                              final double log10GlobalReadMismappingRate) {
        switch ( likelihoodArgs.likelihoodEngineImplementation) {
            case PairHMM:
                return new PairHMMLikelihoodCalculationEngine((byte) likelihoodArgs.gcpHMM, likelihoodArgs.pairHMM, log10GlobalReadMismappingRate, likelihoodArgs.pcrErrorModel, likelihoodArgs.BASE_QUALITY_SCORE_THRESHOLD);
            case Random:
                return new RandomLikelihoodCalculationEngine();
            default:
                throw new UserException("Unsupported likelihood calculation engine.");
        }
    }

    public static double calculateLog10GlobalReadMismappingRate(final int phredScaledGlobalReadMismappingRate) {
        return phredScaledGlobalReadMismappingRate < 0 ? - Double.MAX_VALUE : QualityUtils.qualToErrorProbLog10(phredScaledGlobalReadMismappingRate);
    }

    // create the assembly using just high quality reads (eg Q20 or higher).  We may want to use lower
    // quality reads in the PairHMM downstream, so we can't use a ReadFilter
    public static AssemblyRegion assemblyRegionWithWellMappedReads(final AssemblyRegion originalAssemblyRegion,
                                                                   final int minMappingQuality,
                                                                   final SAMFileHeader readsHeader) {
        final AssemblyRegion result = new AssemblyRegion(originalAssemblyRegion.getSpan(), originalAssemblyRegion.getSupportingStates(), originalAssemblyRegion.isActive(), originalAssemblyRegion.getExtension(), readsHeader);
        originalAssemblyRegion.getReads().stream()
                .filter(rec -> rec.getMappingQuality() >= minMappingQuality)
                .forEach(result::add);
        return result;
    }
}
