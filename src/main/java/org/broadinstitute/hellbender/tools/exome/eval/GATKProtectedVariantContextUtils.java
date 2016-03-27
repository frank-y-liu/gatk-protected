package org.broadinstitute.hellbender.tools.exome.eval;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.function.Supplier;
import java.util.function.ToDoubleFunction;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Miscellaneous utilities to work with Variant data structures.
 */
public class GATKProtectedVariantContextUtils {

    /**
     * Composes the double array from a genotype annotation.
     *
     * @param variantContext the target variant-context.
     * @param key the name of the attribute containing the double array.
     * @param defaultValue the double array to return in case there is no such an annotation.
     * @param missingValue value to use to fill up positions with a missing value (e.g. '.').
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code variantContext} is {@code null} or {@code key} is {@code null}.
     */
    public static double[] getAttributeAsDoubleArray(final VariantContext variantContext, final String key,
                                                     final Supplier<double[]> defaultValue, final double missingValue) {
        Utils.nonNull(variantContext);
        return attributeValueToDoubleArray(variantContext.getAttribute(key), key, defaultValue, missingValue);
    }

    /**
     * Composes the double array from a genotype annotation.
     *
     * @param genotype the target variant-context.
     * @param key the name of the attribute containing the double array.
     * @param defaultValue the double array to return in case there is no such an annotation.
     * @param missingValue value to use to fill up positions with a missing value (e.g. '.').
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code variantContext} is {@code null} or {@code key} is {@code null}.
     */
    public static double[] getAttributeAsDoubleArray(final Genotype genotype, final String key,
                                                     final Supplier<double[]> defaultValue, final double missingValue) {
        Utils.nonNull(genotype);
        return attributeValueToDoubleArray(genotype.getExtendedAttribute(key), key, defaultValue, missingValue);
    }

    private static double[] attributeValueToDoubleArray(final Object value, final String key, final Supplier<double[]> defaultResult, final double missingValue) {
        Utils.nonNull(key);
        final ToDoubleFunction<Object> doubleConverter = o -> {
            if (o == null) {
                return missingValue;
            } else {
                final String s = String.valueOf(o);
                if (s.equals(VCFConstants.MISSING_VALUE_v4)) {
                    return missingValue;
                } else {
                    try {
                        return Double.parseDouble(s);
                    } catch (final NumberFormatException ex) {
                        throw new GATKException(String.format("INFO annotation '%s' contains a non-double value '%s'", key, s), ex);
                    }
                }
            }
        };

        if (value == null) {
            return defaultResult.get();
        } else if (value.getClass().isArray()) {
            final double[] result = new double[Array.getLength(value)];
            for (int i = 0; i < result.length; i++) {
                result[i] = doubleConverter.applyAsDouble(String.valueOf(Array.get(value, i)));
            }
            return result;
        } else if (value.getClass().isAssignableFrom(Iterable.class)) {
            return StreamSupport.stream(((Iterable<?>)value).spliterator(), false)
                    .mapToDouble(doubleConverter).toArray();
        } else { // as a last resort with transform it into an String and try to parse an array out of it.
            return Stream.of(String.valueOf(value).trim().replaceAll("\\[\\]", "")
                    .split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR))
                    .mapToDouble(doubleConverter).toArray();
        }
    }

    public static String getAttributeAsString(final Genotype genotype, final String key, final String defaultValue) {
        Utils.nonNull(genotype);
        Utils.nonNull(key);
        final Object value = genotype.getExtendedAttribute(key);
        if (value == null) {
            return defaultValue;
        } else {
            return String.valueOf(value);
        }
    }

    public static int getAttributeAsInt(final Genotype genotype, final String key, final int defaultValue) {
        Utils.nonNull(genotype);
        Utils.nonNull(key);
        final Object value = genotype.getExtendedAttribute(key);
        if (value == null) {
            return defaultValue;
        } else {
            try {
                return Integer.parseInt(String.valueOf(value));
            } catch (final NumberFormatException ex) {
                throw new IllegalArgumentException(
                        String.format("attribute '%s' does not have a valid integer value: '%s'", key, String.valueOf(value)));
            }
        }
    }

    public static String[] getAttributeAsStringArray(final Genotype genotype, final String key, final Supplier<String[]> defaultValue, final String missingValue) {
        Utils.nonNull(key);
        Utils.nonNull(genotype);
        final Object value = genotype.getExtendedAttribute(key);
        if (value == null) {
            return defaultValue.get();
        } else if (value.getClass().isArray()) {
            final String[] result = new String[Array.getLength(value)];
            for (int i = 0; i < result.length; i++) {
                result[i] = String.valueOf(Array.get(value, i));
                if (result[i].equals(VCFConstants.MISSING_VALUE_v4)) {
                    result[i] =  missingValue;
                }
            }
            return result;
        } else if (value.getClass().isAssignableFrom(Iterable.class)) {
            return StreamSupport.stream(((Iterable<?>)value).spliterator(), false)
                    .map(String::valueOf)
                    .map(s -> s.equals(VCFConstants.MISSING_VALUE_v4) ? missingValue : s).toArray(String[]::new);
        } else { // as a last resort with transform it into an String and try to parse an array out of it.
            return Stream.of(String.valueOf(value).trim().replaceAll("\\[\\]", "")
                    .split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR))
                    .map(String::valueOf)
                    .map(s -> s.equals(VCFConstants.MISSING_VALUE_v4) ? missingValue : s).toArray(String[]::new);
        }
    }
}
