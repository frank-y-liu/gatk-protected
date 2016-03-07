package org.broadinstitute.hellbender.tools.exome;

import com.google.common.annotations.VisibleForTesting;
import com.google.protobuf.Internal;
import htsjdk.samtools.*;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegrator;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegratorFactory;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import scala.Char;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A Bayesian heterozygous SNP pulldown calculator. Base qualities are taken into account
 * to increase precision (see CNV-methods.pdf for details).
 *
 * TODO
 *
 * - The quadrature order can be adaptively chosen based on the pileup size and the accuracy
 *   required for likelihood estimation. In theory, a Gaussian quadrature of order N yields
 *   the exact result for a pileup of size N.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */

public final class BayesianHetPulldownCalculator {

    private final Logger logger = LogManager.getLogger(HetPulldownCalculator.class);

    private final File refFile;
    private final IntervalList snpIntervals;

    private final int readDepthThreshold;
    private final int minMappingQuality;
    private final int minBaseQuality;
    private final ValidationStringency validationStringency;

    private final double errorProbabilityAdjustmentFactor;

    /* these are for building a prior for allele fraction */
    private final double minAbnormalFraction;
    private final double maxAbnormalFraction;
    private final double maxCopyNumber;

    /* these are derived */
    private final double minHetAlleleFraction;
    private final double breakpointHetAlleleFraction;

    /* integration quadrature */
    @VisibleForTesting
    final ArrayList<Double> gaussIntegrationWeights = new ArrayList<>();
    @VisibleForTesting
    final ArrayList<Double> gaussIntegrationLogWeights = new ArrayList<>();
    @VisibleForTesting
    final ArrayList<Double> gaussIntegrationAbscissas = new ArrayList<>();

    /* allele fraction prior for Het sites */
    @VisibleForTesting
    final ArrayList<Double> alleleFractionPriors = new ArrayList<>();
    @VisibleForTesting
    final ArrayList<Double> alleleFractionLogPriors = new ArrayList<>();

    /* minimum order of the integration quadrature */
    private static final int MIN_QUADRATURE_ORDER = 50;

    /* interval threshold for indexing for SamLocusIterator */
    private static final int MAX_INTERVALS_FOR_INDEX = 25000;

    /* default priors */
    private static final double DEFAULT_PRIOR_REF_HOM = 0.5; /* a homozygous site being the ref allele */
    private static final double DEFAULT_PRIOR_HET = 0.5; /* a site being heterozygous */

    /* a third of the minimum sequencing error (for safeguarding log likelihood calculations) */
    private static final double DEFAULT_MIN_BASE_ERROR_THIRD = 1e-6;

    /* approximate number of status updates printed to log */
    private static final int NUMBER_OF_LOG_UPDATES = 20;

    /**
     * a simple class to handle base read and mapping error probabilitites
     */
    @VisibleForTesting
    public static final class BaseQuality {

        final double readErrorProb, mappingErrorProb;

        public BaseQuality(final double readErrorProb, final double mappingErrorProb) {
            this.readErrorProb = readErrorProb;
            this.mappingErrorProb = mappingErrorProb;
        }

        public double getReadErrorProb() { return readErrorProb; }
        public double getMappingErrorProb() { return mappingErrorProb; }

    }

    public BayesianHetPulldownCalculator(final File refFile, final File snpFile,
                                         final int minMappingQuality, final int minBaseQuality,
                                         final int readDepthThreshold, final ValidationStringency validationStringency,
                                         final int quadratureOrder,
                                         final double minAbnormalFraction, final double maxAbnormalFraction,
                                         final int maxCopyNumber, final double errorProbabilityAdjustmentFactor) {

        ParamUtils.isPositiveOrZero(minMappingQuality, "Minimum mapping quality must be nonnegative.");
        ParamUtils.isPositiveOrZero(minBaseQuality, "Minimum base quality must be nonnegative.");

        /* read related members */
        this.refFile = refFile;
        this.snpIntervals = IntervalList.fromFile(snpFile);
        this.minMappingQuality = ParamUtils.isPositive(minMappingQuality, "Minimum mapping quality must be a positive integer");
        this.minBaseQuality = ParamUtils.isPositive(minBaseQuality, "Minimum base quality must be a positive integer");
        this.readDepthThreshold = ParamUtils.isPositive(readDepthThreshold, "Read depth threshold must be a positive integer");
        this.validationStringency = validationStringency;

        /* parameters for building the allele ratio prior at Het sites */
        this.minAbnormalFraction = ParamUtils.inRange(minAbnormalFraction, 0.0, 1.0, "Minimum fraction of abnormal" +
                " cells must be between 0 and 1.");
        this.maxAbnormalFraction = ParamUtils.inRange(maxAbnormalFraction, this.minAbnormalFraction, 1.0, "Maximum fraction of abnormal" +
                " cells must be greater than the provided minimum and less than 1.");
        this.maxCopyNumber = ParamUtils.isPositive(maxCopyNumber, "Maximum copy number must be positive");
        this.errorProbabilityAdjustmentFactor = ParamUtils.isPositive(errorProbabilityAdjustmentFactor,
                "Error adjustment factor must be positive.");

        /* auxiliary member functions */
        this.minHetAlleleFraction = (1 - this.maxAbnormalFraction) / (this.maxCopyNumber * this.maxAbnormalFraction +
                2 * (1 - this.maxAbnormalFraction));
        this.breakpointHetAlleleFraction = (1 - this.minAbnormalFraction) / (this.maxCopyNumber * this.minAbnormalFraction +
                2 * (1 - this.minAbnormalFraction));

        /* initialize the integration quadrature and calculate the allele fraction prior on the abscissas */
        initializeIntegrationQuadrature(ParamUtils.isPositive(quadratureOrder - MIN_QUADRATURE_ORDER,
                "Quadrature order must be greater than " + MIN_QUADRATURE_ORDER) + MIN_QUADRATURE_ORDER);
        initializeHetAlleleFractionPrior();
    }

    /**
     * core computational methods
     */

    /**
     * Initilizes the quadrature for calculating allele ratio integrals in getHetLogLikelihood
     * @param numIntegPoints  number of points in the quadrature
     */
    private void initializeIntegrationQuadrature(final int numIntegPoints) {

        /* get Gauss-Legendre quadrature factory of order @numIntegPoints */
        GaussIntegratorFactory integratorFactory = new GaussIntegratorFactory();
        GaussIntegrator gaussIntegrator = integratorFactory.legendre(numIntegPoints,
                minHetAlleleFraction, 1.0 - minHetAlleleFraction);

        /* abscissas */
        gaussIntegrationAbscissas.addAll(IntStream.range(0, numIntegPoints).
                mapToDouble(gaussIntegrator::getPoint).boxed().collect(Collectors.toList()));

        /* weights */
        gaussIntegrationWeights.addAll(IntStream.range(0, numIntegPoints).
                mapToDouble(gaussIntegrator::getWeight).boxed().collect(Collectors.toList()));

        /* log of weights */
        gaussIntegrationLogWeights.addAll(gaussIntegrationWeights.stream().
                mapToDouble(FastMath::log).boxed().collect(Collectors.toList()));
    }

    /**
     * (advanced) calculate a simple prior (and log prior) for the allele fraction based on
     * (1) minimum purity of the sample, and (2) maximum copy number. See CNV-methods.pdf for details.
     */
    private void initializeHetAlleleFractionPrior() {

        alleleFractionPriors.addAll(gaussIntegrationAbscissas.stream()
                .mapToDouble(alleleFraction -> {
                    double minorAlleleFraction = (alleleFraction < 0.5) ? alleleFraction : 1 - alleleFraction;
                    if (minorAlleleFraction < minHetAlleleFraction) {
                        return 0;
                    } else if (minorAlleleFraction < breakpointHetAlleleFraction) {
                        double denom = 2 * FastMath.pow((1 - minorAlleleFraction) * minorAlleleFraction * maxCopyNumber, 2) *
                                maxAbnormalFraction * (maxAbnormalFraction - minAbnormalFraction);
                        double num = (-1 + (-1 + minorAlleleFraction * maxCopyNumber) * maxAbnormalFraction) *
                                (-1 + maxAbnormalFraction + minorAlleleFraction * (2 + (-2 + maxCopyNumber) * maxAbnormalFraction)) +
                                2 * (1 + minorAlleleFraction * (-2 + minorAlleleFraction * maxCopyNumber)) * maxAbnormalFraction *
                                        FastMath.log(FastMath.abs(((1 + minorAlleleFraction * (-2 + maxCopyNumber)) * maxAbnormalFraction)) /
                                                (1 - 2 * minorAlleleFraction));
                        return num / denom;
                    } else { /* breakpointHetAlleleFraction < minorAlleleFraction < 1/2 */
                        double denom = 2 * FastMath.pow((1 - minorAlleleFraction) * minorAlleleFraction * maxCopyNumber, 2) *
                                maxAbnormalFraction * minAbnormalFraction * (maxAbnormalFraction - minAbnormalFraction);
                        double num = (maxAbnormalFraction - minAbnormalFraction) * (-1 + 2 * minorAlleleFraction +
                                (1 + minorAlleleFraction * (-2 + maxCopyNumber)) * (-1 + minorAlleleFraction * maxCopyNumber) *
                                maxAbnormalFraction * minAbnormalFraction) + 2 * (1 + minorAlleleFraction * (-2 + minorAlleleFraction *
                                maxCopyNumber)) * maxAbnormalFraction * minAbnormalFraction * FastMath.log(maxAbnormalFraction / minAbnormalFraction);
                        return num / denom;
                    }
                })
                .boxed().collect(Collectors.toList()));

        /* calculate the log prior */
        alleleFractionLogPriors.addAll(alleleFractionPriors.stream()
                .mapToDouble(FastMath::log).boxed().collect(Collectors.toList()));
    }

    /**
     * calculate the log probability, safeguarded with a minimum probability DEFAULT_MIN_BASE_ERROR_THIRD
     * @param prob probability
     * @return safeguarded-log of probability
     */
    private double getSafeguardedLogProbability(final double prob) {
        return FastMath.log(FastMath.max(DEFAULT_MIN_BASE_ERROR_THIRD, prob));
    }

    /**
     * calculate the log likelihood of a SNP site being homozygous for a given read pileup
     * (see CNV-method.pdf for details)
     * @param baseQualities map of bases to list of their calling error probabilities at the SNP site
     * @param alleleRef the ref allele base
     * @param alleleAlt the alt allele base
     * @param homRefPrior the prior probability of the ref allele given that the site is homozygous
     * @return the log likelihood
     */
    @VisibleForTesting
    public double getHomLogLikelihood(final Map<Character, ArrayList<BaseQuality>> baseQualities,
                                      final Character alleleRef, final Character alleleAlt,
                                      final double homRefPrior) {

        /* initilize the log likehoods of ref and alt with the priors ... */
        double homRefLogLikelihood = FastMath.log(Utils.nonNull(homRefPrior));
        double homAltLogLikelihood = FastMath.log(1.0 - Utils.nonNull(homRefPrior));

        /* ... and add the log likelihood of the reads */
        for (Character base : baseQualities.keySet()) {
            for (BaseQuality currentBaseQuality : baseQualities.get(base)) {
                homRefLogLikelihood += getSafeguardedLogProbability(currentBaseQuality.getReadErrorProb() / 3 +
                        currentBaseQuality.getMappingErrorProb() / 4 +
                        ((base == alleleRef) ? (1 - 4 * currentBaseQuality.getReadErrorProb() / 3
                                - currentBaseQuality.getMappingErrorProb()) : 0));
                homAltLogLikelihood += getSafeguardedLogProbability(currentBaseQuality.getReadErrorProb() / 3 +
                                currentBaseQuality.getMappingErrorProb() / 4 +
                                ((base == alleleAlt) ? (1 - 4 * currentBaseQuality.getReadErrorProb() / 3
                                        - currentBaseQuality.getMappingErrorProb()) : 0));
            }
        }

        /* return the sum of |hom,ref) and |hom,alt) likelihoods */
        return GATKProtectedMathUtils.naturalLogSumExp(homRefLogLikelihood, homAltLogLikelihood);
    }

    /**
     * [internal helper function]
     * calculate the log likelihood of hetrozygosity from just alt and ref pileup for a given @alleleFraction
     * (see CNV-method.pdf for details)
     * @param alleleFraction ref-to-alt allele fraction
     * @param alphaList [internal to getHetLogLikelihood]
     * @param betaList [internal to getHetLogLikelihood]
     * @return log likelihood
     */
    private double getRefAltHetLogLikelihoodFixedAlleleFraction(final double alleleFraction,
                                                                final ArrayList<Double> alphaList,
                                                                final ArrayList<Double> betaList) {
        double logLikelihood = 0;
        for (int i=0; i<alphaList.size(); i++) {
            logLikelihood += getSafeguardedLogProbability(alphaList.get(i) + alleleFraction * betaList.get(i));
        }

        return logLikelihood;
    }

    /**
     * calculate the log likelihood of a SNP site being heterozygous for a given read pileup
     * (see CNV-method.pdf for details)
     * @param baseQualities map of bases to list of their calling error probabilities at the SNP site
     * @param alleleRef the ref allele base
     * @param alleleAlt the alt allele base
     * @return the log likelihood
     */
    @VisibleForTesting
    public double getHetLogLikelihood(final Map<Character, ArrayList<BaseQuality>> baseQualities,
                                       final Character alleleRef, final Character alleleAlt) {

        double errorLogLikelihood = 0.0;
        double refAltLogLikelihood;
        ArrayList<Double> alphaList = new ArrayList<>();
        ArrayList<Double> betaList = new ArrayList<>();

        for (Character base : baseQualities.keySet()) {
            if (base == alleleRef) {
                for (BaseQuality currentBaseQuality : baseQualities.get(base)) {
                    alphaList.add(currentBaseQuality.getReadErrorProb() / 3 +
                        currentBaseQuality.getMappingErrorProb() / 4);
                    betaList.add(1 - 4 * currentBaseQuality.getReadErrorProb() / 3 -
                        currentBaseQuality.getMappingErrorProb() / 4);
                }
            }
            else if (base == alleleAlt) {
                for (BaseQuality currentBaseQuality : baseQualities.get(base)) {
                    alphaList.add(1 - currentBaseQuality.getReadErrorProb() -
                        3 * currentBaseQuality.getMappingErrorProb() / 4);
                    betaList.add(-1 + 4 * currentBaseQuality.getReadErrorProb() / 3 +
                        currentBaseQuality.getMappingErrorProb());
                }
            }
            else {
                for (BaseQuality currentBaseQuality : baseQualities.get(base)) {
                    errorLogLikelihood += getSafeguardedLogProbability(currentBaseQuality.getReadErrorProb() / 3 +
                        currentBaseQuality.getMappingErrorProb() / 4);
                }
            }
        }

        /** marginalize allele fraction by integrating {@link getRefAltHetLogLikelihoodFixedAlleleFraction} over
         * the precomputed allele fraction prior */
        ArrayList<Double> refAltLogLikelihoodList = new ArrayList<>(gaussIntegrationAbscissas.size());
        refAltLogLikelihoodList.addAll(
                gaussIntegrationAbscissas.stream()
                .mapToDouble(Double::doubleValue)
                .map(f -> getRefAltHetLogLikelihoodFixedAlleleFraction(f, alphaList, betaList))
                .boxed().collect(Collectors.toList()));

        ArrayList<Double> logLikelihoodIntegrandWithPriorAndWeights = new ArrayList<>(gaussIntegrationAbscissas.size());
        logLikelihoodIntegrandWithPriorAndWeights.addAll(
                IntStream.range(0, gaussIntegrationAbscissas.size())
                .mapToDouble(i -> refAltLogLikelihoodList.get(i) + gaussIntegrationLogWeights.get(i) +
                        alleleFractionLogPriors.get(i))
                .boxed().collect(Collectors.toList()));

        refAltLogLikelihood = GATKProtectedMathUtils.naturalLogSumExp(logLikelihoodIntegrandWithPriorAndWeights);

        return errorLogLikelihood + refAltLogLikelihood;
    }

    /**
     * wrapper methods
     */

    /**
     * Returns map of base-pair to error probabilities at a given locus. All reads are considered (not just ACTG)
     * @param locus locus
     * @return map of base-pair to error probabilities
     */
    private Map<Character, ArrayList<BaseQuality>> getPileupBaseQualities(final SamLocusIterator.LocusInfo locus) {

        Map<Character, ArrayList<BaseQuality>> baseQualities = new HashMap<>();

        /* make sure that we have keys at least for A, C, T, and G */
        for (Character base : new Character[]{'A', 'C', 'T', 'G'}) {
            baseQualities.put(base, new ArrayList<>());
        }

        for (final SamLocusIterator.RecordAndOffset rec : locus.getRecordAndPositions()) {
            Character currentBase = (char)(Character.toUpperCase(rec.getReadBase() & 0xFF));
            ArrayList<BaseQuality> currentBaseQualityList = baseQualities.get(currentBase);
            if (currentBaseQualityList == null) { /* this can occur if a non-ACTG read is found */
                currentBaseQualityList = new ArrayList<>();
                baseQualities.put(currentBase, currentBaseQualityList);
            }
            currentBaseQualityList.add(new BaseQuality(
                    errorProbabilityAdjustmentFactor * QualityUtils.qualToErrorProb(rec.getBaseQuality()),
                    errorProbabilityAdjustmentFactor * QualityUtils.qualToErrorProb(rec.getRecord().getMappingQuality()))
            );
        }

        return baseQualities;
    }

    /**
     * get base counts map from base quality map
     * @param baseQualities map from bases to list of {@link BaseQuality}
     * @return base counts map
     */
    private Map<Character, Integer> getBaseCountsFromBaseQualities(final Map<Character, ArrayList<BaseQuality>> baseQualities) {
        Map<Character, Integer> baseCounts = new HashMap<>();
        for (Character base : new Character[]{'A', 'C', 'T', 'G'}) {
            if (baseQualities.containsKey(base)) {
                baseCounts.put(base, baseQualities.get(base).size());
            } else {
                baseCounts.put(base, 0);
            }
        }
        return baseCounts;
    }

    /**
     * guess the alt allele base from the pileup. we simply choose the base with the highest frequency
     *
     * * Remark: while this is definitely not the right way to do it, it has no serious pitfalls:
     *
     * - if the input is balanced between betweet ref and alt, it correctly chooses alt.
     * - if the input is all ref, then alt is chosen randomly, but it's OK since there are no alt counts
     *   in the pileup anyway.
     * - if the input is all alt, we are good.
     * - if the input in all ref + a few errors, then alt will be choosen as the most frequent erroneous
     *   base; again, this is OK because Hom likelihood is far greater than Het likelihood in this case
     * - if the input is all alt + a few errors, we are good.

     * @param baseQualities map from bases to error probabilities
     * @param refBase the ref allele base
     * @return the likely alt allele
     */
    @VisibleForTesting
    public static char guessAltFromPileup(final Map<Character, ArrayList<BaseQuality>> baseQualities,
                                          final char refBase) {

        /* sort the bases in the descending order by their frequency */
        Character[] bases = new Character[]{'A', 'C', 'T', 'G'};
        Arrays.sort(bases, (L, R) -> Integer.compare(baseQualities.get(R).size(), baseQualities.get(L).size()));
        /* pick the base with highest frequency, skip over ref */
        for (char base : bases) {
            if (base != refBase) {
                return base;
            }
        }

        /* we shouldn't be here unless the baseErrorProbabilities is malformed */
        return 'N';
    }

    /**
     * For a normal or tumor sample, returns a data structure giving (intervals, reference counts, alternate counts),
     * where intervals give positions of likely heterozygous SNP sites.
     *
     * The analysis for normal and tumor is identical. The IntervalList gives common SNP sites in 1-based
     * format.
     *
     * The @hetCallingStrigency parameters sets the threshold Het posterior for calling:
     *
     *      hetPosteriorThreshold = 1 - 10^{-hetCallingStringency}
     *      hetThresholdLogOdds = log(hetPosteriorThreshold/(1-hetPosteriorThreshold))
     *                          = log(10^{hetCallingStringency} - 1)
     *
     * @param bamFile sorted BAM file for sample
     * @param hetCallingStringency strigency for calling a Het site
     * @return Pulldown of heterozygous SNP sites in 1-based format
     */
    @VisibleForTesting
    public Pulldown getHetPulldown(final File bamFile, final double hetCallingStringency) {

        final double hetThresholdLogOdds = FastMath.log(FastMath.pow(10, hetCallingStringency) - 1);

        try (final SamReader bamReader = SamReaderFactory.makeDefault().validationStringency(validationStringency)
                .referenceSequence(refFile).open(bamFile);
             final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(this.refFile)) {
            if (bamReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
                throw new UserException.BadInput("BAM file " + bamFile.toString() + " must be coordinate sorted.");
            }

            final Pulldown hetPulldown = new Pulldown(bamReader.getFileHeader());

            final int totalNumberOfSNPs = snpIntervals.size();
            final SamLocusIterator locusIterator = new SamLocusIterator(bamReader, snpIntervals,
                    totalNumberOfSNPs < MAX_INTERVALS_FOR_INDEX);

            /* set read and locus filters */
            final List<SamRecordFilter> samFilters = Arrays.asList(new NotPrimaryAlignmentFilter(),
                    new DuplicateReadFilter());
            locusIterator.setSamFilters(samFilters);
            locusIterator.setEmitUncoveredLoci(false);
            locusIterator.setIncludeNonPfReads(false);
            locusIterator.setMappingQualityScoreCutoff(minMappingQuality);
            locusIterator.setQualityScoreCutoff(minBaseQuality);

            logger.info("Examining " + totalNumberOfSNPs + " sites...");

            final int iterationsPerLogUpdate =
                    Math.max((int) Math.floor((float) totalNumberOfSNPs / NUMBER_OF_LOG_UPDATES), 1);

            int locusCount = 1;
            for (final SamLocusIterator.LocusInfo locus : locusIterator) {

                if (locusCount % iterationsPerLogUpdate == 0) {
                    logger.info("Examined " + locusCount + " out of " + totalNumberOfSNPs + " sites.");
                }
                locusCount++;

                final int totalReadCount = locus.getRecordAndPositions().size();
                if (totalReadCount <= readDepthThreshold) {
                    continue;
                }

                final Map<Character, ArrayList<BaseQuality>> baseQualities = getPileupBaseQualities(locus);

                final char refBase = (char) refWalker.get(locus.getSequenceIndex()).getBases()[locus.getPosition() - 1];

                /**
                 * At the moment, we don't have the alt allele information (Sam has promised to create a new common SNP
                 * list with the minor allele). To keep the ball rolling, I will take the next common base in the read
                 * as the alt allele base.
                 *
                 * (please refer to {@link guessAltFromPileup} for details)
                 */
                final char altBase = guessAltFromPileup(baseQualities, refBase);

                /* calculate Het log odds */
                final double hetLogLikelihood = getHetLogLikelihood(baseQualities, refBase, altBase);
                final double homLogLikelihood = getHomLogLikelihood(baseQualities, refBase, altBase,
                        DEFAULT_PRIOR_REF_HOM);
                final double hetLogOdds = (hetLogLikelihood + FastMath.log(DEFAULT_PRIOR_HET)) -
                        (homLogLikelihood + FastMath.log(1 - DEFAULT_PRIOR_HET));

                if (hetLogOdds > hetThresholdLogOdds) {
                    hetPulldown.add(new SimpleInterval(locus.getSequenceName(), locus.getPosition(), locus.getPosition()),
                            baseQualities.get(refBase).size(), baseQualities.get(altBase).size(),
                            getBaseCountsFromBaseQualities(baseQualities), hetLogLikelihood, homLogLikelihood);
                }
            }

            logger.info("Examined " + totalNumberOfSNPs + " sites.");

            return hetPulldown;

        } catch (final IOException | SAMFormatException e) {
            throw new UserException(e.getMessage());
        }
    }

    /* this is for debugging purposes -- will be removed */
    public static void printLocusInfo(final Map<Character, ArrayList<Double>> baseErrorProbabilities,
                                      final double hetLogLikelihood, final double homLogLikelihood) {

        System.out.println('A');
        baseErrorProbabilities.get('A').stream().forEach(System.out::println);
        System.out.println('C');
        baseErrorProbabilities.get('C').stream().forEach(System.out::println);
        System.out.println('T');
        baseErrorProbabilities.get('T').stream().forEach(System.out::println);
        System.out.println('G');
        baseErrorProbabilities.get('G').stream().forEach(System.out::println);

        System.out.println("Het log likelihood:");
        System.out.println(hetLogLikelihood);

        System.out.println("Hom log likelihood:");
        System.out.println(homLogLikelihood);

    }
}