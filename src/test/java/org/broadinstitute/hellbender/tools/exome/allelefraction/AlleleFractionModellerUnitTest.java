package org.broadinstitute.hellbender.tools.exome.allelefraction;

import htsjdk.samtools.util.Log;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.Decile;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Tests the MCMC inference of the {@link AlleleFractionModeller}.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionModellerUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final File ALLELIC_PON_NORMAL_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-pon-normal.tsv");
    private static final File ALLELIC_PON_EVENT_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-pon-event.tsv");

    private static final File SAMPLE_NORMAL_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-sample-normal.tsv");
    private static final File SAMPLE_EVENT_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-sample-event.tsv");

    private static final File SEGMENTS_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-segments.seg");

    @Test
    public void testMCMCWithoutAllelicPON() {
        final double meanBias = 1.1;
        final double biasVariance = 0.01;

        testMCMC(meanBias, biasVariance, AllelicPanelOfNormals.EMPTY_PON);
    }

    /**
     * Note that these MCMC tests were written to use simulated hets before the allelic PON was introduced.
     * Rather than generate a simulated PON on the fly, we simply use a fixed PON loaded from a file
     * and check that its hyperparameters are recovered correctly---i.e., the PON is not actually used to
     * to correct reference bias in the simulated data in any way. This latter functionality is tested
     * below instead.
     * TODO change so that simulated data is made with different alpha and beta than PON and that PON values are recovered
     */
    @Test
    public void testMCMCWithAllelicPON() {
        final double meanBias = 1.083;
        final double biasVariance = 0.0181;
        final AllelicPanelOfNormals allelicPON = new AllelicPanelOfNormals(ALLELIC_PON_NORMAL_FILE);
        testMCMC(meanBias, biasVariance, allelicPON);
    }

    private void testMCMC(final double meanBias, final double biasVariance, final AllelicPanelOfNormals allelicPON) {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);

        final int numSamples = 300;
        final int numBurnIn = 100;

        final double averageHetsPerSegment = 50;
        final int numSegments = 100;
        final int averageDepth = 50;


        final double outlierProbability = 0.02;

        // note: the following tolerances could actually be made much smaller if we used more segments and/or
        // more hets -- most of the error is the sampling error of a finite simulated data set, not numerical error of MCMC
        final double minorFractionTolerance = 0.02;
        final double meanBiasTolerance = 0.02;
        final double biasVarianceTolerance = 0.01;
        final double outlierProbabilityTolerance = 0.02;
        final AlleleFractionSimulatedData simulatedData = new AlleleFractionSimulatedData(averageHetsPerSegment, numSegments,
                averageDepth, meanBias, biasVariance, outlierProbability);

        final AlleleFractionModeller model = new AlleleFractionModeller(simulatedData.getSegmentedModel(), allelicPON);
        model.fitMCMC(numSamples, numBurnIn);

        final List<Double> meanBiasSamples = model.getmeanBiasSamples();
        Assert.assertEquals(meanBiasSamples.size(), numSamples - numBurnIn);

        final List<Double> biasVarianceSamples = model.getBiasVarianceSamples();
        Assert.assertEquals(biasVarianceSamples.size(), numSamples - numBurnIn);

        final List<Double> outlierProbabilitySamples = model.getOutlierProbabilitySamples();
        Assert.assertEquals(outlierProbabilitySamples.size(), numSamples - numBurnIn);

        final List<AlleleFractionState.MinorFractions> minorFractionsSamples = model.getMinorFractionsSamples();
        Assert.assertEquals(minorFractionsSamples.size(), numSamples - numBurnIn);
        for (final AlleleFractionState.MinorFractions sample : minorFractionsSamples) {
            Assert.assertEquals(sample.size(), numSegments);
        }

        final List<List<Double>> minorFractionsSamplesBySegment = model.getMinorFractionSamplesBySegment();

        final double mcmcMeanBias = meanBiasSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final double mcmcBiasVariance = biasVarianceSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final double mcmcOutlierProbabilityr = outlierProbabilitySamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final List<Double> mcmcMinorFractions = minorFractionsSamplesBySegment
                .stream().map(list -> list.stream().mapToDouble(x -> x).average().getAsDouble())
                .collect(Collectors.toList());

        double totalSegmentError = 0.0;
        for (int segment = 0; segment < numSegments; segment++) {
            totalSegmentError += Math.abs(mcmcMinorFractions.get(segment) - simulatedData.getTrueState().minorFractionInSegment(segment));
        }

        Assert.assertEquals(mcmcMeanBias, meanBias, meanBiasTolerance);
        Assert.assertEquals(mcmcBiasVariance, biasVariance, biasVarianceTolerance);
        Assert.assertEquals(mcmcOutlierProbabilityr, outlierProbability, outlierProbabilityTolerance);
        Assert.assertEquals(totalSegmentError / numSegments, 0.0, minorFractionTolerance);
    }

    @Test
    private void testEvent() {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final int numSamples = 100;
        final int numBurnIn = 25;

        final AllelicPanelOfNormals allelicPON = new AllelicPanelOfNormals(ALLELIC_PON_EVENT_FILE);
        final AllelicCountCollection sample = new AllelicCountCollection(SAMPLE_EVENT_FILE);
        final List<TargetCoverage> emptyTargets = new ArrayList<>();
        final Genome genome = new Genome(emptyTargets, sample.getCounts(), "test");
        final List<SimpleInterval> segments = SegmentUtils.readIntervalsFromSegmentFile(SEGMENTS_FILE);
        final SegmentedModel segmentedModel = new SegmentedModel(segments, genome);

        final AlleleFractionModeller modeller = new AlleleFractionModeller(segmentedModel, allelicPON);
        modeller.fitMCMC(numSamples, numBurnIn);

        final List<PosteriorSummary> minorAlleleFractionPosteriorSummaries =
                modeller.getMinorAlleleFractionsPosteriorSummaries(0.05, ctx);
        System.out.println(minorAlleleFractionPosteriorSummaries.stream().map(s -> s.getDeciles().getAll()).collect(Collectors.toList()));
    }
}