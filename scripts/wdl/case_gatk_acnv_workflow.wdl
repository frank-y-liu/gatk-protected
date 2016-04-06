#
# Workflow a single pair of case-control samples.
#
# Notes:
#
# - file names will use the entity ID specified, but inside the file, the bam SM tag will typically be used.
#
# - THIS SCRIPT SHOULD BE CONSIDERED OF "BETA" QUALITY
#
###########

workflow case_gatk_acnv_workflow {
    String wf_entity_id_tumor
    String wf_entity_id_normal
    File target_bed
    File ref_fasta
    File ref_fasta_dict
    File ref_fasta_fai
    File common_snp_list
    File tumor_bam
    File tumor_bam_idx
    File normal_bam
    File normal_bam_idx
    File jar_file
    File PoN
    String is_disable_reference_validation

  call PadTargets {
    input:
        target_bed=target_bed,
        jar_file=jar_file
  }

  call CalculateTargetCoverage {
    input:
        entity_id=wf_entity_id_tumor,
        padded_target_bed=PadTargets.padded_target_bed,
        input_bam=tumor_bam,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_fai,
        ref_fasta_dict=ref_fasta_dict,
        jar_file=jar_file,
        disable_reference_validation=is_disable_reference_validation
  }

  call NormalizeSomaticReadCounts {
    input:
        entity_id=wf_entity_id_tumor,
        coverage_file=CalculateTargetCoverage.gatk_cnv_coverage_file,
        padded_target_bed=PadTargets.padded_target_bed,
        pon=PoN,
        jar_file=jar_file
  }

  call PerformSegmentation {
    input:
        entity_id=wf_entity_id_tumor,
        jar_file=jar_file,
        tn_file=NormalizeSomaticReadCounts.tn_file,
        mem=2
  }

  call Caller {
    input:
        entity_id=wf_entity_id_tumor,
        jar_file=jar_file,
        tn_file=NormalizeSomaticReadCounts.tn_file,
        seg_file=PerformSegmentation.seg_file,
        mem=2
  }

  call HetPulldown {
    input:
        entity_id_tumor=wf_entity_id_tumor,
        entity_id_normal=wf_entity_id_normal,
        jar_file=jar_file,
        mem = 4,
        ref_fasta=ref_fasta,
        ref_fasta_fai=ref_fasta_fai,
        ref_fasta_dict=ref_fasta_dict,
        tumor_bam=tumor_bam,
        tumor_bam_idx=tumor_bam_idx,
        normal_bam=normal_bam,
        normal_bam_idx=normal_bam_idx,
        common_snp_list=common_snp_list
  }

  call AllelicCNV {
    input:
        entity_id=wf_entity_id_tumor,
        jar_file=jar_file,
        mem = 4,
        tumor_hets=HetPulldown.tumor_hets,
        called_file=Caller.called_file,
        tn_file=NormalizeSomaticReadCounts.tn_file
  }
}

# Pad the target file.  This was found to help sensitivity and specificity.  This step should only be altered
#  by advanced users.  Note that by changing this, you need to have a PoN that also reflects the change.
task PadTargets {
    File target_bed
    Int pd = 250
    Int mem = 1
    File jar_file
    command {
        java -Xmx${mem}g -Djava.library.path=/usr/lib/jni/ -jar ${jar_file} PadTargets  --targets ${target_bed} --output targets.padded.bed --padding ${pd}  --help false --version false --verbosity INFO --QUIET false
    }
    output {
        File padded_target_bed = "targets.padded.bed"
    }
    #runtime {
    #    docker: "gatk-protected/a1"
    #}
}

# Calculate the target coverage
task CalculateTargetCoverage {
    String entity_id
    File padded_target_bed
    String transform = "PCOV"
    String grouping = "SAMPLE"
    Boolean keepduplicatereads = true
    Boolean disable_all_read_filters = false
    File input_bam
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    Int mem = 4
    File jar_file
    String disable_reference_validation = "false"

    command {
        java -Xmx${mem}g -Djava.library.path=/usr/lib/jni/ -jar ${jar_file} CalculateTargetCoverage --output ${entity_id}.coverage.tsv --groupBy ${grouping} --transform ${transform} --targets ${padded_target_bed} --targetInformationColumns FULL --keepduplicatereads ${keepduplicatereads} --input ${input_bam} --reference ${ref_fasta}  --disable_all_read_filters ${disable_all_read_filters} --interval_set_rule UNION --interval_padding 0 --secondsBetweenProgressUpdates 10.0 --disableSequenceDictionaryValidation ${disable_reference_validation} --createOutputBamIndex true --help false --version false --verbosity INFO --QUIET false
    }
    output {
        File gatk_cnv_coverage_file = "${entity_id}.coverage.tsv"
    }

    #runtime {
    #    docker: "gatk-protected/a1"
    #}
}

task NormalizeSomaticReadCounts {
    String entity_id
    File coverage_file
    File padded_target_bed
    File pon
    Int mem = 2
    File jar_file
    command {
        java -Xmx${mem}g -Djava.library.path=/usr/lib/jni/ -jar ${jar_file} NormalizeSomaticReadCounts  --input ${coverage_file} \
        --targets ${padded_target_bed} --panelOfNormals ${pon} --factorNormalizedOutput ${entity_id}.fnt.tsv --tangentNormalized ${entity_id}.tn.tsv \
        --betaHatsOutput ${entity_id}.betaHats.tsv --preTangentNormalized  ${entity_id}.preTN.tsv  --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File tn_file = "${entity_id}.tn.tsv"
        File pre_tn_file = "${entity_id}.preTN.tsv"
        File betahats_file = "${entity_id}.betaHats.tsv"
    }
    #runtime {
    #    docker: "gatk-protected/a1"
    #}
}

task PerformSegmentation {
    Int mem = 2
    String entity_id
    File jar_file
    File tn_file

    command {
        java -Xmx${mem}g -Djava.library.path=/usr/lib/jni/ -jar ${jar_file} PerformSegmentation  --targets ${tn_file} \
        --output ${entity_id}.seg --log2Input true  --alpha 0.01 --nperm 10000 --pmethod HYBRID --minWidth 2 --kmax 25 \
        --nmin 200 --eta 0.05 --trim 0.025 --undoSplits NONE --undoPrune 0.05 --undoSD 3 --help false --version false \
        --verbosity INFO --QUIET false
    }

    output {
        File seg_file = "${entity_id}.seg"
    }
}

task Caller {
    Int mem = 2
    String entity_id
    File jar_file
    File tn_file
    File seg_file

    command {
        java -Xmx${mem}g -Djava.library.path=/usr/lib/jni/ -jar ${jar_file} CallSegments  --targets ${tn_file} \
         --segments ${seg_file} --output ${entity_id}.called --threshold 2.0  --legacy false --experimental false \
          --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File called_file="${entity_id}.called"
    }
}

# entity IDs can be the same value
task HetPulldown {
    String entity_id_tumor
    String entity_id_normal
    File jar_file
    Int mem = 2
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File tumor_bam
    File tumor_bam_idx
    File normal_bam
    File normal_bam_idx
    File common_snp_list

    command {
        java -Xmx${mem}g -Djava.library.path=/usr/lib/jni/ -jar ${jar_file} GetHetCoverage  --reference ${ref_fasta} \
        --normal ${normal_bam} --tumor ${tumor_bam} --snpIntervals ${common_snp_list}  \
        --normalHets ${entity_id_normal}.normal.hets.tsv --tumorHets ${entity_id_tumor}.tumor.hets.tsv --pvalueThreshold 0.05 \
        --help false --version false --verbosity INFO --QUIET false --VALIDATION_STRINGENCY LENIENT
    }

    output {
        File normal_hets="${entity_id}.normal.hets.tsv"
        File tumor_hets="${entity_id}.tumor.hets.tsv"
    }
}

task AllelicCNV {
    String entity_id
    File jar_file
    Int mem = 2
    File tumor_hets
    File called_file
    File tn_file

    command {
        java -Xmx${mem}g -Djava.library.path=/usr/lib/jni/ -jar ${jar_file} AllelicCNV  --tumorHets ${tumor_hets} \
         --tangentNormalized ${tn_file} --segments ${called_file} --outputPrefix ${entity_id} \
         --smallSegmentThreshold 3 --numSamplesCopyRatio 100 --numBurnInCopyRatio 50 --numSamplesAlleleFraction 100 \
         --numBurnInAlleleFraction 50 --intervalThresholdCopyRatio 5.0 --intervalThresholdAlleleFraction 2.0  \
         --help false --version false --verbosity INFO --QUIET false
    }

    output {
        File acnv_final_segs="${entity_id}-sim-final.seg"
        File acnv_final_segs_gatk_compatible="${entity_id}-sim-final.cnv.seg"
        File acnv_final_segs_acs_compatible="${entity_id}-sim-final.acs.seg"
    }
}