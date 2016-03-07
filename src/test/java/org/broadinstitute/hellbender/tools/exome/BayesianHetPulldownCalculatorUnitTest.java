package org.broadinstitute.hellbender.tools.exome;

import com.sun.tools.javac.util.List;
import htsjdk.samtools.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Unit tests for {@link BayesianHetPulldownCalculator}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class BayesianHetPulldownCalculatorUnitTest extends BaseTest {

    private static final String TEST_SUB_DIR = publicTestDir + "org/broadinstitute/hellbender/tools/exome";

    /*****************************************************************************************************************/

    private static final File NORMAL_BAM_FILE = new File(TEST_SUB_DIR, "normal.sorted.bam");
    private static final File SNP_FILE = new File(TEST_SUB_DIR, "common_SNP.interval_list");
    private static final File TUMOR_BAM_FILE = new File(TEST_SUB_DIR, "tumor.sorted.bam");
    private static final File REF_FILE = new File(hg19MiniReference);

    private static final int READ_DEPTH_THRESHOLD = 10;
    private static final int MINIMUM_MAPPING_QUALITY = 30;
    private static final int MINIMUM_BASE_QUALITY = 20;
    private static final int QUADRATURE_ORDER = 200;
    private static final double MIN_ABNORMAL_FRACTION = 0.5;
    private static final double MAX_ABNORMAL_FRACTION = 0.8;
    private static final int MAX_COPY_NUMBER = 4;
    private static final double ERROR_PROBABILITY_ADJUSTMENT_FACTOR = 1.0;

    private static SAMFileHeader normalHeader;
    private static SAMFileHeader tumorHeader;

    private static final BayesianHetPulldownCalculator calculator = new BayesianHetPulldownCalculator(
            REF_FILE, SNP_FILE, MINIMUM_MAPPING_QUALITY, MINIMUM_BASE_QUALITY, READ_DEPTH_THRESHOLD,
            ValidationStringency.STRICT, QUADRATURE_ORDER, MIN_ABNORMAL_FRACTION, MAX_ABNORMAL_FRACTION,
            MAX_COPY_NUMBER, ERROR_PROBABILITY_ADJUSTMENT_FACTOR);

    private static int numPileupEntries;
    private static ArrayList<Map<Character, ArrayList<BayesianHetPulldownCalculator.BaseQuality>>>
            fakePileupBaseQualities = new ArrayList<>();
    private static ArrayList<Double> fakePileupHetLogLikelihoodArray = new ArrayList<>();
    private static ArrayList<Double> fakePileupHomLogLikelihoodArray = new ArrayList<>();

    /*****************************************************************************************************************/

    private static final File FAKE_PILEUP_FILE = new File(TEST_SUB_DIR, "test_fake_generated_pileup.txt");

    public interface StringMapper <T> {
        T of(String input);
    }

    @BeforeClass
    public void initHeaders() throws IOException {
        try (final SamReader normalBamReader = SamReaderFactory.makeDefault().open(NORMAL_BAM_FILE);
             final SamReader tumorBamReader = SamReaderFactory.makeDefault().open(TUMOR_BAM_FILE)) {
            normalHeader = normalBamReader.getFileHeader();
            tumorHeader = tumorBamReader.getFileHeader();
        }
    }

    /**
     * load the fake fileup from file
     * TODO add mapping qualities to the fake pileup; all is set to 60.0 at the moment
     *
     * @throws IOException
     */
    @BeforeClass
    public void loadPileup() throws IOException {

        /* load fake pileup and likelihoods from disk */
        StringMapper<String> stringStripper = s -> s.substring(s.indexOf(':') + 1).replaceAll("\\s+", "");
        StringMapper<Double> parseToDouble = s -> Double.parseDouble(stringStripper.of(s));
        StringMapper<ArrayList<Double>> parseToDoubleArray = s -> {
            String[] tokenizedLine = stringStripper.of(s).split(",");
            ArrayList<Double> errorList = new ArrayList<>();
            if (tokenizedLine.length >= 1 && !tokenizedLine[0].equals("")) {
                errorList.addAll(List.from(tokenizedLine).stream()
                        .map(Double::parseDouble)
                        .collect(Collectors.toList()));
            }
            return errorList;
        };

        Scanner reader = new Scanner(new FileInputStream(FAKE_PILEUP_FILE));
        while (reader.hasNextLine()) {

            Map<Character, ArrayList<BayesianHetPulldownCalculator.BaseQuality>> baseQualities = new HashMap<>();

            for (Character base : new Character[]{'A', 'C', 'T', 'G'}) {
                /* load base read error list */
                ArrayList<Double> readErrorList = parseToDoubleArray.of(reader.nextLine());

                /* TODO load mapping error list */
                ArrayList<Double> mappingErrorList = new ArrayList<>();
                mappingErrorList.addAll(IntStream.range(0, readErrorList.size())
                        .mapToDouble(i -> 60.0).boxed().collect(Collectors.toList()));

                /* contruct the BaseQuality list*/
                ArrayList<BayesianHetPulldownCalculator.BaseQuality> baseQualityList = new ArrayList<>();
                baseQualityList.addAll(IntStream.range(0, readErrorList.size()).mapToObj(i -> new BayesianHetPulldownCalculator.BaseQuality(
                        readErrorList.get(i), mappingErrorList.get(i))).collect(Collectors.toList()));

                baseQualities.put(base, baseQualityList);
            }

            fakePileupBaseQualities.add(baseQualities);
            fakePileupHetLogLikelihoodArray.add(parseToDouble.of(reader.nextLine()));
            fakePileupHomLogLikelihoodArray.add(parseToDouble.of(reader.nextLine()));
        }

        numPileupEntries = fakePileupBaseQualities.size();

    }

    @Test
    public void testQuadrature() {
        double sumWeights = calculator.gaussIntegrationWeights.stream().mapToDouble(Double::doubleValue).sum();
        double minHetAlleleFraction = (1 - MAX_ABNORMAL_FRACTION) / (MAX_COPY_NUMBER * MAX_ABNORMAL_FRACTION +
                2 * (1 - MAX_ABNORMAL_FRACTION));
        Assert.assertEquals(sumWeights, 1 - 2 * minHetAlleleFraction, 1e-8);
        Assert.assertEquals(calculator.gaussIntegrationAbscissas.size(), QUADRATURE_ORDER);
    }

    @Test
    public void testAlelleRatioPrior() {
        /* the prior should integrate to 1 */
        Assert.assertEquals(IntStream.range(0, calculator.gaussIntegrationAbscissas.size())
                .mapToDouble(i -> calculator.gaussIntegrationWeights.get(i) * calculator.alleleFractionPriors.get(i))
                .sum(), 1.0, 1e-3);
        calculator.alleleFractionPriors.stream().forEach(x -> Assert.assertTrue(x > 0));
        // calculator.gaussIntegrationAbscissas.stream().forEach(System.out::println);
        // calculator.alleleFractionPriors.stream().forEach(System.out::println);
    }

    @DataProvider(name = "inputTestGetHomLogLikelihood")
    public Object[][] inputTestGetHomLogLikelihood() {

        Character alleleRef = 'A';
        Character alleleAlt = 'T';
        double homRefPrior = 0.5;

        Object[][] out = new Object[numPileupEntries][];
        for (int i=0; i<numPileupEntries; i++) {
            out[i] = new Object[]{
                    fakePileupBaseQualities.get(i),
                    alleleRef,
                    alleleAlt,
                    homRefPrior,
                    fakePileupHomLogLikelihoodArray.get(i)
            };
        }

        return out;
    }

//    @Test(dataProvider = "inputTestGetHomLogLikelihood")
//    public void testGetHomLogLikelihood(final Map<Character, ArrayList<Double>> baseErrorProbabilities,
//                                        final Character alleleRef, final Character alleleAlt,
//                                        final double homRefPrior, final double expectedHomLogLikelihood) {
//        for (int i=0; i<numPileupEntries; i++) {
//            Assert.assertEquals(calculator.getHomLogLikelihood(baseErrorProbabilities, alleleRef, alleleAlt,
//                    homRefPrior), expectedHomLogLikelihood, 1e-4);
//        }
//    }
//
//    @DataProvider(name = "inputTestGetHetLogLikelihood")
//    public Object[][] inputTestGetHetLogLikelihood() {
//
//        Character alleleRef = 'A';
//        Character alleleAlt = 'T';
//
//        Object[][] out = new Object[numPileupEntries][];
//        for (int i=0; i<numPileupEntries; i++) {
//            out[i] = new Object[]{
//                    fakePileupBaseErrorProbabilitiesList.get(i),
//                    alleleRef,
//                    alleleAlt,
//                    fakePileupHetLogLikelihoodArray.get(i)
//            };
//        }
//
//        return out;
//    }
//
//    @Test(dataProvider = "inputTestGetHetLogLikelihood")
//    public void testGetHetLogLikelihood(final Map<Character, ArrayList<Double>> baseErrorProbabilities,
//                                        final Character alleleRef, final Character alleleAlt,
//                                        final double expectedHetLogLikelihood) {
//        for (int i=0; i<numPileupEntries; i++) {
//            Assert.assertEquals(calculator.getHetLogLikelihood(baseErrorProbabilities, alleleRef, alleleAlt),
//                    expectedHetLogLikelihood, 1e-2);
//        }
//    }

//    @Test
//    public void testGuessAlleleAltFromPileup() {
//
//        /* test 1 */
//        Map<Character, ArrayList<Double>> baseErrorProbabilities = new HashMap<>();
//        baseErrorProbabilities.put('A', new ArrayList<>(Arrays.asList(new Double[]{
//                0.0, 0.5, 0.4, 0.0
//        })));
//        baseErrorProbabilities.put('C', new ArrayList<>(Arrays.asList(new Double[]{
//                0.0, 0.1
//        })));
//        baseErrorProbabilities.put('T', new ArrayList<>(Arrays.asList(new Double[]{
//                0.2
//        })));
//        baseErrorProbabilities.put('G', new ArrayList<>(Arrays.asList(new Double[]{
//                0.0
//        })));
//
//        Assert.assertEquals(BayesianHetPulldownCalculator.guessAltFromPileup(baseErrorProbabilities, 'A'), 'C');
//        Assert.assertEquals(BayesianHetPulldownCalculator.guessAltFromPileup(baseErrorProbabilities, 'C'), 'A');
//        Assert.assertEquals(BayesianHetPulldownCalculator.guessAltFromPileup(baseErrorProbabilities, 'T'), 'A');
//        Assert.assertEquals(BayesianHetPulldownCalculator.guessAltFromPileup(baseErrorProbabilities, 'G'), 'A');
//
//        /* test 2 */
//        baseErrorProbabilities.put('A', new ArrayList<>(Arrays.asList(new Double[]{
//        })));
//        baseErrorProbabilities.put('C', new ArrayList<>(Arrays.asList(new Double[]{
//                0.0, 0.1, 0.0
//        })));
//        baseErrorProbabilities.put('T', new ArrayList<>(Arrays.asList(new Double[]{
//        })));
//        baseErrorProbabilities.put('G', new ArrayList<>(Arrays.asList(new Double[]{
//        })));
//
//        Assert.assertEquals(BayesianHetPulldownCalculator.guessAltFromPileup(baseErrorProbabilities, 'A'), 'C');
//
//        /* in this case, A is chosen simply because sorting places A right after C; in theory, T and G are equally
//           valid "guesses"
//         */
//        Assert.assertEquals(BayesianHetPulldownCalculator.guessAltFromPileup(baseErrorProbabilities, 'C'), 'A');
//
//        Assert.assertEquals(BayesianHetPulldownCalculator.guessAltFromPileup(baseErrorProbabilities, 'T'), 'C');
//        Assert.assertEquals(BayesianHetPulldownCalculator.guessAltFromPileup(baseErrorProbabilities, 'G'), 'C');
//
//    }

    @Test
    public void testGetHetPulldown() {

        Pulldown testPulldown, expectedPulldown;
        File tempFile;

        /* write Pulldown to file while checking the results */
        try {

            /* test 1: normal, loose threshold */
            testPulldown = calculator.getHetPulldown(NORMAL_BAM_FILE, 2);
            tempFile = File.createTempFile("testPulldownNormalLoose", ".txt");
            testPulldown.writeWithBaseCountsAndLikelihoods(tempFile);

            expectedPulldown = new Pulldown(normalHeader);
            expectedPulldown.add(new SimpleInterval("1", 11522, 11522), 7, 4);
            expectedPulldown.add(new SimpleInterval("1", 12098, 12098), 8, 6);
            expectedPulldown.add(new SimpleInterval("1", 14630, 14630), 9, 8);
            expectedPulldown.add(new SimpleInterval("2", 14689, 14689), 6, 9);
            expectedPulldown.add(new SimpleInterval("2", 14982, 14982), 6, 5);

            Assert.assertEquals(testPulldown, expectedPulldown);

            /* test 2: normal, tight threshold */
            testPulldown = calculator.getHetPulldown(NORMAL_BAM_FILE, 12);
            tempFile = File.createTempFile("testPulldownNormalTight", ".txt");
            testPulldown.writeWithBaseCountsAndLikelihoods(tempFile);

            expectedPulldown = new Pulldown(normalHeader);
            expectedPulldown.add(new SimpleInterval("1", 12098, 12098), 8, 6);
            expectedPulldown.add(new SimpleInterval("1", 14630, 14630), 9, 8);
            expectedPulldown.add(new SimpleInterval("2", 14689, 14689), 6, 9);

            Assert.assertEquals(testPulldown, expectedPulldown);

            /* test 3: tumor, loose threshold */
            testPulldown = calculator.getHetPulldown(TUMOR_BAM_FILE, 2);
            tempFile = File.createTempFile("testPulldownTumorLoose", ".txt");
            testPulldown.writeWithBaseCountsAndLikelihoods(tempFile);

            expectedPulldown = new Pulldown(tumorHeader);
            expectedPulldown.add(new SimpleInterval("1", 11522, 11522), 7, 4);
            expectedPulldown.add(new SimpleInterval("1", 12098, 12098), 8, 6);
            expectedPulldown.add(new SimpleInterval("1", 14630, 14630), 9, 8);
            expectedPulldown.add(new SimpleInterval("2", 14689, 14689), 6, 9);
            expectedPulldown.add(new SimpleInterval("2", 14982, 14982), 6, 5);

            Assert.assertEquals(testPulldown, expectedPulldown);

            /* test 4: tumor, tight threshold */
            testPulldown = calculator.getHetPulldown(TUMOR_BAM_FILE, 12);
            tempFile = File.createTempFile("testPulldownTumorTight", ".txt");
            testPulldown.writeWithBaseCountsAndLikelihoods(tempFile);

            expectedPulldown = new Pulldown(tumorHeader);
            expectedPulldown.add(new SimpleInterval("1", 14630, 14630), 9, 8);

            Assert.assertEquals(testPulldown, expectedPulldown);

        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not write pulldown to to file.", e);
        }

    }
}
