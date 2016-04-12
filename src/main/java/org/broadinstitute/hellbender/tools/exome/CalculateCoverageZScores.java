package org.broadinstitute.hellbender.tools.exome;

import org.apache.zookeeper.Op;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * Calculates Z scores of coverage.
 * <p>
 *   This tool takes a read counts file and outputs a read counts file in which the coverage unit is the Z score
 *   relative to the distribution of input coverage at the corresponding target.
 * </p>
 * <p>
 *   The input and output format for the coverage file is described in {@link ReadCountCollectionUtils}.
 * </p>
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Expresses coverage of a read counts file in terms of Z scores.",
        oneLineSummary = "Expresses coverage of a read counts file in terms of Z scores.",
        programGroup = CopyNumberProgramGroup.class
)
public final class CalculateCoverageZScores extends CommandLineProgram {
    @Argument(
            doc = "Input coverage file",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            optional = false
    )
    protected File inputFile;

    @Argument(
            doc = "Output coverage file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            optional = false
    )
    protected File outputFile;

    @Argument(
            doc = "File containing post-tangent normalization weights (inverse variance) of each target in the panel of normals." +
                    "  By default, this simple text file is written next to the pon file with the extension " + CreatePanelOfNormals.TARGET_WEIGHTS_FILE_APPEND + "." +
            "If no file given, per-sample z scores will be calculated instead of per-target z-scores.",
            shortName = CreatePanelOfNormals.TARGET_WEIGHTS_SHORT_NAME,
            fullName = CreatePanelOfNormals.TARGET_WEIGHTS_FULL_NAME,
            optional = true
    )
    protected File targetWeightsFile;

    @Override
    public Object doWork() {
        final ReadCountCollection inputCounts = getInputCounts();
        final ReadCountCollection zScoreCounts = calculateZScores(inputCounts);
        writeZScoreCounts(zScoreCounts);
        return "SUCCESS";
    }

    private ReadCountCollection getInputCounts() {
        final ReadCountCollection inputCounts;
        try {
            inputCounts = ReadCountCollectionUtils.parse(inputFile);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, "cannot reach or read the input target coverage file");
        }
        return inputCounts;
    }

    private ReadCountCollection calculateZScores(final ReadCountCollection inputCounts) {
        final ReadCountCollection zScoreCounts;
        try {
            if (targetWeightsFile == null) {
                zScoreCounts = inputCounts.zScoreCountsBySample();
            } else {
                final double[] targetInverseVariances = ParamUtils.readValuesFromFile(targetWeightsFile);
                final double[] targetStandardDeviations = Arrays.stream(targetInverseVariances)
                        .map(x -> 1.0 / Math.sqrt(x)).toArray();
                final double[] zeroMeans = new double[targetInverseVariances.length];   // zero mean for tangent-normalized counts
                zScoreCounts = inputCounts.zScoreCounts(zeroMeans, targetStandardDeviations);
            }
        } catch (final IllegalArgumentException e) {
            throw new UserException.BadInput("Could not read target weights file.");
        }
        return zScoreCounts;
    }

    private void writeZScoreCounts(final ReadCountCollection zScoreCounts) {
        try {
            ReadCountCollectionUtils.write(outputFile, zScoreCounts);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, "cannot  write to the given output file");
        }
    }

}
