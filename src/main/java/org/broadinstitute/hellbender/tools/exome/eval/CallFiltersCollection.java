/*
*  By downloading the PROGRAM you agree to the following terms of use:
*
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
*
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/
package org.broadinstitute.hellbender.tools.exome.eval;

import org.broadinstitute.hellbender.cmdline.Argument;

import java.util.function.Predicate;

/**
 * TODO document this.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class CallFiltersCollection {

    public static final String MINIMUM_TRUTH_SEGMENT_LENGTH_SHORT_NAME = "minTruthLen";
    public static final String MINIMUM_TRUTH_SEGMENT_LENGTH_FULL_NAME = "minimumTruthSegmentLength";
    public static final String MINIMUM_CALLED_SEGMENT_LENGTH_SHORT_NAME = "minCallLen";
    public static final String MINIMUM_CALLED_SEGMENT_LENGTH_FULL_NAME = "minimumCalledSegmentLength";
    public static final String MINIMUM_TRUTH_SEGMENT_QUALITY_SHORT_NAME = "minTruthQual";
    public static final String MINIMUM_TRUTH_SEGMENT_QUALITY_FULL_NAME = "minimumTruthSegmentQuality";
    public static final String MINIMUM_CALLED_SEGMENT_QUALITY_SHORT_NAME = "minCallQual";
    public static final String MINIMUM_CALLED_SEGMENT_QUALITY_FULL_NAME = "minimumCalledSegmentQuality";
    public static final String APPLY_MULTI_ALLELIC_TRUTH_FILTER_SHORT_NAME = "applyMATFilter";
    public static final String APPLY_MULTI_ALLELIC_TRUTH_FILTER_FULL_NAME = "applyMultiAllelicTruthFilter";
    public static final String APPLY_MULTI_ALLELIC_CALLED_FILTER_SHORT_NAME = "applyMACFilter";
    public static final String APPLY_MULTI_ALLELIC_CALLED_FILTER_FULL_NAME = "applyMultiAllelicCalledFilter";
    public static final String MAXIMUM_TRUTH_EVENT_FREQUENCY_SHORT_NAME = "maxTruthFreq";
    public static final String MAXIMUM_TRUTH_EVENT_FREQUENCY_FULL_NAME = "maximumTruthEventFrequency";
    public static final String MAXIMUM_CALLED_EVENT_FREQUENCY_SHORT_NAME = "maxCallFreq";
    public static final String MAXIMUM_CALLED_EVENT_FREQUENCY_FULL_NAME = "maximumCalledEventFrequency";

    public static final int DEFAULT_MINIMUM_TRUTH_SEGMENT_LENGTH = 1;
    public static final int DEFAULT_MINIMUM_CALLED_SEGMENT_LENGTH = 1;
    public static final double DEFAULT_MINIMUM_TRUTH_SEGMENT_QUALITY = 0.0;
    public static final double DEFAULT_MINIMUM_CALLED_SEGMENT_QUALITY = 0.0;
    public static final double DEFAULT_MAXIMUM_TRUTH_EVENT_FREQUENCY = 1.0;
    public static final double DEFAULT_MAXIMUM_CALLED_EVENT_FREQUENCY = 1.0;

    @Argument(
            doc = "Minimum Truth segment length to consider the segment callable.",
            shortName = MINIMUM_TRUTH_SEGMENT_LENGTH_SHORT_NAME,
            fullName = MINIMUM_TRUTH_SEGMENT_LENGTH_FULL_NAME,
            optional = true
    )
    protected int minimumTruthSegmentLength = DEFAULT_MINIMUM_TRUTH_SEGMENT_LENGTH;

    @Argument(
            doc = "Minimum called segment length to consider the segment trust-worthy.",
            shortName = MINIMUM_CALLED_SEGMENT_LENGTH_SHORT_NAME,
            fullName = MINIMUM_CALLED_SEGMENT_LENGTH_FULL_NAME,
            optional = true
    )
    protected int minimumCalledSegmentLength = DEFAULT_MINIMUM_CALLED_SEGMENT_LENGTH;

    @Argument(
            doc = "Minimum Truth segment quality.",
            shortName = MINIMUM_TRUTH_SEGMENT_QUALITY_SHORT_NAME,
            fullName = MINIMUM_TRUTH_SEGMENT_QUALITY_FULL_NAME,
            optional = true
    )
    protected double minimumTruthSegmentQuality = DEFAULT_MINIMUM_TRUTH_SEGMENT_QUALITY;

    @Argument(
            doc = "Minimum called segment quality.",
            shortName = MINIMUM_CALLED_SEGMENT_QUALITY_SHORT_NAME,
            fullName = MINIMUM_CALLED_SEGMENT_QUALITY_FULL_NAME,
            optional = true
    )
    protected double minimumCalledSegmentQuality = DEFAULT_MINIMUM_CALLED_SEGMENT_QUALITY;

    @Argument(
            doc = "Apply the multi-allelic truth segment filter",
            shortName = APPLY_MULTI_ALLELIC_TRUTH_FILTER_SHORT_NAME,
            fullName = APPLY_MULTI_ALLELIC_TRUTH_FILTER_FULL_NAME,
            optional = true
    )
    protected boolean applyMultiAllelicTruthFilter = false;

    @Argument(
            doc = "Apply the multi-allelic called segment filter",
            shortName = APPLY_MULTI_ALLELIC_CALLED_FILTER_SHORT_NAME,
            fullName = APPLY_MULTI_ALLELIC_CALLED_FILTER_FULL_NAME,
            optional = true
    )
    protected boolean applyMultiAllelicCalledFilter = false;

    @Argument(
            doc = "Maximum frequency for truth events at a position",
            shortName = MAXIMUM_TRUTH_EVENT_FREQUENCY_SHORT_NAME,
            fullName = MAXIMUM_TRUTH_EVENT_FREQUENCY_FULL_NAME,
            optional = true
    )
    protected double maximumTruthEventFrequency = DEFAULT_MAXIMUM_TRUTH_EVENT_FREQUENCY;

    @Argument(
            doc = "Maximum frequency for truth events at a position",
            shortName = MAXIMUM_CALLED_EVENT_FREQUENCY_SHORT_NAME,
            fullName = MAXIMUM_CALLED_EVENT_FREQUENCY_FULL_NAME,
            optional = true
    )
    protected double maximumCalledEventFrequency = DEFAULT_MAXIMUM_CALLED_EVENT_FREQUENCY;



    public EvaluationSiteRecord applyFiltersOn(
            final EvaluationSiteRecord record) {

    }

    /**
     * Evaluation filters.
     */
    public enum Filter {
        /**
         * The truth event segment is too small (target count) to be easy to call; so won't count towards FN.
         * <p>Controlled by {@link #minimumTruthSegmentLength}.</p>
         */
        ShortTrueSegment,

        /**
         * The called event segment is too short (target count wise) to be reliable; so won't count towards positives
         * evaluation classes (TP, UP, DP... etc.)
         * <p>Controlled by {@link #minimumCalledSegmentLength}.</p>
         */
        ShortCalledSegment,

        /**
         * The truth event segment has a low quality so no easy to call or even not a real event; so won't count towards anything.
         * <p>Controlled by {@link #minimumTruthSegmentQuality}.</p>
         */
        LowTrueSegmentQuality,

        /**
         * The called segment has a low quality so it won't count towards positive evaluation classes.
         * <p>Controlled by {@link #minimumCalledSegmentQuality}.</p>
         */
        LowCallQuality,

        /**
         * The are a few called segments across sample with divergent copy number changes (mixture of del. and dups.).
         * So it won't count towards positive evaluation classes.
         * <p>Controlled by {@link #applyMultiAllelicCalledFilter}.</p>
         */
        MultiAllelicCall,

        /**
         * The are a few called segments across sample with divergent copy number changes (mixture of del. and dups.).
         * So it won't count towards false negatives.
         * <p>Controlled by {@link #applyMultiAllelicCalledFilter}.</p>
         */
        MultiAllelicTruth,

        /**
         * Too many samples have overlapping calls.
         * <p>Controlled by {@link #maximumCalledEventFrequency}.</p>
         */
        CommonCalledVariant,

        /**
         * Too many sample have overlapping true segments.
         * <p>Controlled by {@link #maximumTruthEventFrequency}.</p>
         */
        CommonTruthVariant,
    }

    public interface CallFilterPredicate extends Predicate<EvaluationSiteRecord> {

        Filter getFilter();

        @Override
        boolean test(final EvaluationSiteRecord evaluationSiteRecord);
    }

    private class CommonTruthVariant implements CallFilterPredicate {

        @Override
        public Filter getFilter() {
            return Filter.CommonCalledVariant;
        }

        @Override
        public boolean test(final EvaluationSiteRecord evaluationSiteRecord) {
            return maximumCalledEventFrequency >= evaluationSiteRecord.deletionFrequency +
                    evaluationSiteRecord.duplicationFrequency;
        }
    }

    private class MultiAllelicTruth implements CallFilterPredicate {

        @Override
        public Filter getFilter() {
            return Filter.MultiAllelicTruth;
        }

        @Override
        public boolean test(final EvaluationSiteRecord evaluationSiteRecord) {
            return evaluationSiteRecord.deletionFrequency <= 0 || evaluationSiteRecord.duplicationFrequency <= 0;
        }
    }

    private class LowTrueSegmentQuality implements  CallFilterPredicate {

        @Override
        public Filter getFilter() {
            return Filter.LowTrueSegmentQuality;
        }

        @Override
        public boolean test(final EvaluationSiteRecord evaluationSiteRecord) {
            return minimumCalledSegmentQuality <= evaluationSiteRecord.truthQuality();
        }
    }



}
