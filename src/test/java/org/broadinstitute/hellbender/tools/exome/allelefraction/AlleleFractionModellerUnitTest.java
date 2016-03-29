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
import org.testng.annotations.DataProvider;
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

    private static final File SAMPLE_EVENT_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-sample-event.tsv");

    private static final File SEGMENTS_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-segments.seg");

    /**
     * Test MCMC inference on simulated data.  Note that hyperparameter values used to generate the data should be recovered
     * along with outlier probability and minor fractions.
     */
    @Test
    public void testMCMCWithoutAllelicPON() {
        final double meanBiasSimulated = 1.2;
        final double biasVarianceSimulated = 0.04;
        testMCMC(meanBiasSimulated, biasVarianceSimulated, meanBiasSimulated, biasVarianceSimulated, AllelicPanelOfNormals.EMPTY_PON);
    }

    /**
     * Test MCMC inference on simulated data using an allelic PON.  Note that these MCMC tests were written to use
     * simulated hets before the allelic PON was introduced.  Rather than generate a simulated PON on the fly,
     * we simply use a fixed PON loaded from a file and check that its MLE hyperparameters are "sampled" correctly
     * by simply taking the MLE PON values---i.e., the PON does not actually cover the simulated sites and
     * hence is not used to correct reference bias in the simulated data in any way.
     * This latter functionality is tested on fixed data loaded from files in
     * {@link AlleleFractionModellerUnitTest#testBiasCorrection} instead.
     */
    @Test
    public void testMCMCWithAllelicPON() {
        final double meanBiasSimulated = 1.2;
        final double biasVarianceSimulated = 0.04;
        final double meanBiasOfPON = 1.083;         // alpha = 65
        final double biasVarianceOfPON = 0.0181;    // beta = 60
        final AllelicPanelOfNormals allelicPON = new AllelicPanelOfNormals(ALLELIC_PON_NORMAL_FILE);
        testMCMC(meanBiasSimulated, biasVarianceSimulated, meanBiasOfPON, biasVarianceOfPON, allelicPON);
    }

    private void testMCMC(final double meanBiasSimulated, final double biasVarianceSimulated,
                          final double meanBiasExpected, final double biasVarianceExpected,
                          final AllelicPanelOfNormals allelicPON) {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);

        final int numSamples = 100;
        final int numBurnIn = 25;

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
                averageDepth, meanBiasSimulated, biasVarianceSimulated, outlierProbability);

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

        Assert.assertEquals(mcmcMeanBias, meanBiasExpected, meanBiasTolerance);
        Assert.assertEquals(mcmcBiasVariance, biasVarianceExpected, biasVarianceTolerance);
        Assert.assertEquals(mcmcOutlierProbabilityr, outlierProbability, outlierProbabilityTolerance);
        Assert.assertEquals(totalSegmentError / numSegments, 0.0, minorFractionTolerance);
    }

    @DataProvider(name = "biasCorrection")
    public Object[][] dataBiasCorrection() {
        return new Object[][]{
                {new AllelicPanelOfNormals(ALLELIC_PON_EVENT_FILE), 0.5},
                {new AllelicPanelOfNormals(ALLELIC_PON_NORMAL_FILE), 0.4}
        };
    }

    /**
     * Tests that the allelic PoN is appropriately used to correct bias.  The basic set up for the test data is
     * simulated hets at 1000 sites (1:1-1000) across 3 segments in a normal.  In the middle segment consisting of
     * 100 sites (1:451-550), all of the sites have relatively high reference bias
     * (alpha = 9, beta = 6 -> mean bias = 1.5, as compared to alpha = 65, beta = 60 -> mean bias = 1.083 for the
     * rest of the sites).  In this segment, without using a PON that knows about the high reference bias of these sites
     * we would infer a minor-allele fraction of 6 / (6 + 9) = 0.40; however, with the PON we correctly infer that all
     * of the segments are balanced.
     *
     * <p>
     *     Note that alpha and beta are not actually correctly recovered in this PON via MLE because the biases
     *     drawn from a mixture of gamma distributions (as opposed to a single gamma distribution as assumed in the model).
     *     TODO https://github.com/broadinstitute/gatk-protected/issues/421
     * </p>
     */
    @Test(dataProvider = "biasCorrection")
    public void testBiasCorrection(final AllelicPanelOfNormals allelicPON, final double minorFractionExpectedInMiddleSegment) {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final double minorFractionTolerance = 0.025;

        // set up test data
        final AllelicCountCollection sample = new AllelicCountCollection(SAMPLE_EVENT_FILE);
        final List<TargetCoverage> emptyTargets = new ArrayList<>();    // no targets in test data
        final Genome genome = new Genome(emptyTargets, sample.getCounts(), "test");
        final List<SimpleInterval> segments = SegmentUtils.readIntervalsFromSegmentFile(SEGMENTS_FILE);
        final SegmentedModel segmentedModel = new SegmentedModel(segments, genome);

        final int numSamples = 100;
        final int numBurnIn = 25;
        final AlleleFractionModeller modeller = new AlleleFractionModeller(segmentedModel, allelicPON);
        modeller.fitMCMC(numSamples, numBurnIn);

        final double credibleAlpha = 0.05;
        final List<PosteriorSummary> minorAlleleFractionPosteriorSummaries =
                modeller.getMinorAlleleFractionsPosteriorSummaries(credibleAlpha, ctx);
        final List<Double> minorFractionsResult = minorAlleleFractionPosteriorSummaries.stream().map(PosteriorSummary::getCenter).collect(Collectors.toList());

        final double minorFractionBalanced = 0.5;
        final List<Double> minorFractionsExpected = Arrays.asList(minorFractionBalanced, minorFractionExpectedInMiddleSegment, minorFractionBalanced);
        for (int segment = 0; segment < 3; segment++) {
            Assert.assertEquals(minorFractionsResult.get(segment), minorFractionsExpected.get(segment), minorFractionTolerance);
        }
    }
}