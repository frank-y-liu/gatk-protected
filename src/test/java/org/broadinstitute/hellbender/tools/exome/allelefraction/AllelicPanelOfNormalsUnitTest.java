package org.broadinstitute.hellbender.tools.exome.allelefraction;

import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Tests for the {@link AllelicPanelOfNormals}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class AllelicPanelOfNormalsUnitTest extends BaseTest {

    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome/";

    private static final File ALLELIC_PON_NORMAL_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-pon-normal.tsv");
    private static final File ALLELIC_PON_EVENT_FILE = new File(TEST_SUB_DIR, "allelic-pon-test-pon-event.tsv");

    private static final double ALPHA_EXPECTED = 65;
    private static final double BETA_EXPECTED = 60;
    private static final double MEAN_BIAS_EXPECTED = meanBias(ALPHA_EXPECTED, BETA_EXPECTED);
    private static final double BIAS_VARIANCE_EXPECTED = biasVariance(ALPHA_EXPECTED, BETA_EXPECTED);
    private static final double DELTA = 0.5;

    @Test
    public void testNormalPoNHyperparameterInitialization() {
        testPoNHyperparameterInialization(ALLELIC_PON_NORMAL_FILE);
    }

    @Test
    public void testEventPoNHyperparameterInitialization() {
        testPoNHyperparameterInialization(ALLELIC_PON_EVENT_FILE);

    }

    private void testPoNHyperparameterInialization(final File allelicPONFile) {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final AllelicPanelOfNormals allelicPON = new AllelicPanelOfNormals(allelicPONFile);

        final SimpleInterval siteNotInPON = new SimpleInterval("2", 1, 1);  //all sites in PON are from chr1
        final double alphaNotInPON = allelicPON.getAlpha(siteNotInPON);
        final double betaNotInPON = allelicPON.getBeta(siteNotInPON);
        final double meanBias = allelicPON.getMLEMeanBias();
        final double biasVariance = allelicPON.getMLEBiasVariance();

        Assert.assertEquals(alphaNotInPON, ALPHA_EXPECTED, DELTA);
        Assert.assertEquals(betaNotInPON, BETA_EXPECTED, DELTA);
        Assert.assertEquals(meanBias, MEAN_BIAS_EXPECTED, DELTA);
        Assert.assertEquals(biasVariance, BIAS_VARIANCE_EXPECTED, DELTA);
    }

    private static double meanBias(final double alpha, final double beta) {
        return alpha / beta;
    }

    private static double biasVariance(final double alpha, final double beta) {
        return alpha / (beta * beta);
    }
}