package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Test the MCMC inference of the {@link AlleleFractionModeller}.
 *
 *  @author David Benjamin
 */
public final class AlleleFractionModellerUnitTest extends BaseTest {
    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final File ALLELIC_PON_FILE = new File(TEST_SUB_DIR, "allelic-pon-for-acnv-modeller.tsv");

    @Test
    public void testMCMC() {
        final int numSamples = 300;
        final int numBurnIn = 100;

        final double averageHetsPerSegment = 50;
        final int numSegments = 100;
        final int averageDepth = 50;

        final double meanBias = 1.1;
        final double biasVariance = 0.01;
        final double outlierProbability = 0.02;

        // note: the following tolerances could actually be made much smaller if we used more segments and/or
        // more hets -- most of the error is the sampling error of a finite simulated data set, not numerical error of MCMC
        final double minorFractionTolerance = 0.02;
        final double meanBiasTolerance = 0.02;
        final double biasVarianceTolerance = 0.01;
        final double outlierProbabilityTolerance = 0.02;
        final AlleleFractionSimulatedData simulatedData = new AlleleFractionSimulatedData(averageHetsPerSegment, numSegments,
                averageDepth, meanBias, biasVariance, outlierProbability);

        final AlleleFractionModeller model = new AlleleFractionModeller(simulatedData.getSegmentedModel());
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
    public void testMCMCWithAllelicPON() {
        final int numSamples = 300;
        final int numBurnIn = 100;

        final double averageHetsPerSegment = 50;
        final int numSegments = 100;
        final int averageDepth = 50;

        final double meanBias = 1.1;
        final double biasVariance = 0.01;
        final double outlierProbability = 0.02;

        // note: the following tolerances could actually be made much smaller if we used more segments and/or
        // more hets -- most of the error is the sampling error of a finite simulated data set, not numerical error of MCMC
        final double minorFractionTolerance = 0.02;
        final double meanBiasTolerance = 0.02;
        final double biasVarianceTolerance = 0.01;
        final double outlierProbabilityTolerance = 0.02;
        final AlleleFractionSimulatedData simulatedData = new AlleleFractionSimulatedData(averageHetsPerSegment, numSegments,
                averageDepth, meanBias, biasVariance, outlierProbability);

        final AllelicPanelOfNormals allelicPON = new AllelicPanelOfNormals(ALLELIC_PON_FILE);

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
}