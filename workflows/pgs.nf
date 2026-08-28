nextflow.enable.dsl=2

def chromosomeList(value) {
    def text = value.toString()
    if (text ==~ /\d+\.\.\d+/) {
        def bounds = text.split(/\.\./)*.toInteger()
        return (bounds[0]..bounds[1]).toList()
    }
    return text.split(',')*.trim().findAll { it }*.toInteger()
}

def flagEnabled(value) {
    return value instanceof Boolean ? value : value?.toString()?.toBoolean()
}

process PREPROCESS_CHROMOSOME {
    tag "chr${chr}"
    label 'small'
    publishDir "${params.outdir}/01_chromosomes", mode: 'copy', pattern: 'chr*.*'

    input:
    tuple val(chr), path(vcf), path(rsid_map)

    output:
    tuple val(chr), path("chr${chr}.pgen"), path("chr${chr}.pvar"), path("chr${chr}.psam"), emit: pfiles

    script:
    def infoFilter = params.r2 != null ? "--extract-if-info R2 >= ${params.r2}" : (params.aq != null ? "--extract-if-info AQ >= ${params.aq}" : '')
    def memMb = Math.max(1000, task.memory.toMega() - 1000)
    """
    plink2 \
      --vcf ${vcf} \
      --vcf-half-call missing \
      --snps-only just-acgt \
      --max-alleles 2 \
      --var-filter \
      --mac ${params.mac} \
      --geno ${params.geno} \
      --set-all-var-ids '@:#:${'$'}r:${'$'}a' \
      ${infoFilter} \
      --make-pgen \
      --out filtered \
      --threads ${task.cpus} \
      --memory ${memMb}

    awk '{print ${'$'}2}' ${rsid_map} > mapped_ids.txt
    plink2 \
      --pfile filtered \
      --update-name ${rsid_map} 2 1 \
      --extract mapped_ids.txt \
      --maf ${params.maf} \
      --make-pgen \
      --out chr${chr} \
      --threads ${task.cpus} \
      --memory ${memMb}
    """
}

/*
 * Inputs on a shared filesystem are passed as values so Nextflow does not copy
 * chromosome-scale files. Every worker must see the same absolute paths.
 */
process PREPROCESS_CHROMOSOME_DIRECT {
    tag "chr${chr}"
    label 'small'
    publishDir "${params.outdir}/01_chromosomes", mode: 'copy', pattern: 'chr*.*'

    input:
    tuple val(chr), val(vcf), val(rsid_map)

    output:
    tuple val(chr), path("chr${chr}.pgen"), path("chr${chr}.pvar"), path("chr${chr}.psam"), emit: pfiles

    script:
    def infoFilter = params.r2 != null ? "--extract-if-info R2 >= ${params.r2}" : (params.aq != null ? "--extract-if-info AQ >= ${params.aq}" : '')
    def memMb = Math.max(1000, task.memory.toMega() - 1000)
    """
    test -r '${vcf}'
    plink2 \
      --vcf '${vcf}' \
      --vcf-half-call missing \
      --snps-only just-acgt \
      --max-alleles 2 \
      --var-filter \
      --mac ${params.mac} \
      --geno ${params.geno} \
      --set-all-var-ids '@:#:${'$'}r:${'$'}a' \
      ${infoFilter} \
      --make-pgen \
      --out filtered \
      --threads ${task.cpus} \
      --memory ${memMb}

    awk '{print ${'$'}2}' '${rsid_map}' > mapped_ids.txt
    plink2 \
      --pfile filtered \
      --update-name '${rsid_map}' 2 1 \
      --extract mapped_ids.txt \
      --maf ${params.maf} \
      --make-pgen \
      --out chr${chr} \
      --threads ${task.cpus} \
      --memory ${memMb}
    """
}

process PREPROCESS_CHROMOSOME_DIRECT_NO_RSID {
    tag "chr${chr}"
    label 'small'
    publishDir "${params.outdir}/01_chromosomes", mode: 'copy', pattern: 'chr*.*'

    input:
    tuple val(chr), val(vcf)

    output:
    tuple val(chr), path("chr${chr}.pgen"), path("chr${chr}.pvar"), path("chr${chr}.psam"), emit: pfiles

    script:
    def infoFilter = params.r2 != null ? "--extract-if-info R2 >= ${params.r2}" : (params.aq != null ? "--extract-if-info AQ >= ${params.aq}" : '')
    def memMb = Math.max(1000, task.memory.toMega() - 1000)
    """
    test -r '${vcf}'
    plink2 \
      --vcf '${vcf}' \
      --vcf-half-call missing \
      --snps-only just-acgt \
      --max-alleles 2 \
      --var-filter \
      --mac ${params.mac} \
      --geno ${params.geno} \
      --set-all-var-ids '@:#:${'$'}r:${'$'}a' \
      ${infoFilter} \
      --maf ${params.maf} \
      --make-pgen \
      --out chr${chr} \
      --threads ${task.cpus} \
      --memory ${memMb}
    """
}

process CONCAT_CHROMOSOMES {
    tag params.cohort
    label 'large'
    publishDir "${params.outdir}/02_concat", mode: 'copy'

    input:
    path chromosome_files

    output:
    tuple path("${params.cohort}.pgen"), path("${params.cohort}.pvar"), path("${params.cohort}.psam"), emit: pfile

    script:
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    find . -maxdepth 1 -name 'chr*.pgen' -print | sed 's/[.]pgen${'$'}//' | sort -V > merge_list.txt
    test -s merge_list.txt
    if [ "${'$'}(wc -l < merge_list.txt)" -eq 1 ]; then
      plink2 \
        --pfile "${'$'}(head -n 1 merge_list.txt)" \
        --make-pgen \
        --out ${params.cohort} \
        --threads ${task.cpus} \
        --memory ${memMb}
    else
      plink2 \
        --pmerge-list merge_list.txt \
        --make-pgen \
        --out ${params.cohort} \
        --threads ${task.cpus} \
        --memory ${memMb}
    fi
    """
}

process MISSINGNESS_QC {
    tag params.cohort
    label 'large'
    publishDir "${params.outdir}/03_missingness", mode: 'copy'

    input:
    tuple path(pgen), path(pvar), path(psam)

    output:
    tuple path("${params.cohort}.pgen"), path("${params.cohort}.pvar"), path("${params.cohort}.psam"), emit: pfile
    path 'fail_variants.txt', emit: fail_variants
    path 'fail_samples.txt', emit: fail_samples

    script:
    def inputPrefix = pgen.baseName
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    plink2 --pfile ${inputPrefix} --missing --out pass1 --threads ${task.cpus} --memory ${memMb}
    awk 'NR>1 && ${'$'}5>${params.variant_miss} {print ${'$'}2}' pass1.vmiss > fail_variants.txt
    plink2 --pfile ${inputPrefix} --exclude fail_variants.txt --missing --out pass2 --threads ${task.cpus} --memory ${memMb}
    awk 'NR>1 && ${'$'}4>${params.sample_miss} {print ${'$'}1, ${'$'}2}' pass2.smiss > fail_samples.txt
    plink2 \
      --pfile ${inputPrefix} \
      --exclude fail_variants.txt \
      --remove fail_samples.txt \
      --make-pgen \
      --out ${params.cohort} \
      --threads ${task.cpus} \
      --memory ${memMb}
    """
}

process SUMMARY_QC {
    tag params.cohort
    label 'large'
    publishDir "${params.outdir}/04_summary", mode: 'copy'

    input:
    tuple path(pgen), path(pvar), path(psam)

    output:
    path "${params.cohort}.*"
    path 'summary_counts.txt'

    script:
    def inputPrefix = pgen.baseName
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    plink2 --pfile ${inputPrefix} --missing --hardy --freq \
      --out ${params.cohort} --threads ${task.cpus} --memory ${memMb}
    printf 'Metric\\tCount\\nSamples\\t%s\\nVariants\\t%s\\n' \
      "${'$'}(awk 'END {print NR-1}' ${params.cohort}.smiss)" \
      "${'$'}(awk 'END {print NR-1}' ${params.cohort}.vmiss)" > summary_counts.txt
    """
}

process SUMMARY_QC_DIRECT {
    tag params.cohort
    label 'large'
    publishDir "${params.outdir}/04_summary", mode: 'copy'

    input:
    tuple val(pgen), val(pvar), val(psam)

    output:
    path "${params.cohort}.*"
    path 'summary_counts.txt'

    script:
    def inputPrefix = pgen.toString().replaceFirst(/[.]pgen$/, '')
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    test -r '${pgen}' && test -r '${pvar}' && test -r '${psam}'
    plink2 --pfile '${inputPrefix}' --missing --hardy --freq \
      --out ${params.cohort} --threads ${task.cpus} --memory ${memMb}
    printf 'Metric\\tCount\\nSamples\\t%s\\nVariants\\t%s\\n' \
      "${'$'}(awk 'END {print NR-1}' ${params.cohort}.smiss)" \
      "${'$'}(awk 'END {print NR-1}' ${params.cohort}.vmiss)" > summary_counts.txt
    """
}

process SCORE_TRAIT {
    tag trait
    label 'scoring'
    publishDir "${params.outdir}/05_scores", mode: 'copy'

    input:
    tuple val(trait), path(weights), val(id_col), val(allele_col), val(effect_col), path(pgen), path(pvar), path(psam)
    path summary_script

    output:
    tuple val(trait), path("${trait}.sscore"), path("${trait}.score_qc.tsv"), emit: scored
    path "${trait}.sscore.vars"

    script:
    def inputPrefix = pgen.baseName
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    plink2 \
      --pfile ${inputPrefix} \
      --score ${weights} ${id_col} ${allele_col} ${effect_col} header center list-variants no-mean-imputation \
      --out ${trait} \
      --threads ${task.cpus} \
      --memory ${memMb}
    requested=${'$'}(awk 'NR > 1 { n++ } END { print n + 0 }' ${weights})
    matched=${'$'}(wc -l < ${trait}.sscore.vars)
    fraction=${'$'}(awk -v m="${'$'}matched" -v r="${'$'}requested" 'BEGIN { if (r == 0) print 0; else printf "%.8f", m / r }')
    awk -v trait='${trait}' -v requested="${'$'}requested" -v matched="${'$'}matched" \
      -v fraction="${'$'}fraction" -v warn='${params.warn_score_variant_match}' \
      -f ${summary_script} ${trait}.sscore > ${trait}.score_qc.tsv
    awk -v f="${'$'}fraction" -v minimum='${params.min_score_variant_match}' \
      'BEGIN { if (f < minimum) { printf "ERROR: score variant match fraction %.4f is below minimum %.4f\\n", f, minimum > "/dev/stderr"; exit 1 } }'
    """
}

process SCORE_TRAIT_DIRECT {
    tag trait
    label 'scoring'
    publishDir "${params.outdir}/05_scores", mode: 'copy'

    input:
    tuple val(trait), path(weights), val(id_col), val(allele_col), val(effect_col),
      val(pgen), val(pvar), val(psam)
    path summary_script

    output:
    tuple val(trait), path("${trait}.sscore"), path("${trait}.score_qc.tsv"), emit: scored
    path "${trait}.sscore.vars"

    script:
    def inputPrefix = pgen.toString().replaceFirst(/[.]pgen$/, '')
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    test -r '${pgen}' && test -r '${pvar}' && test -r '${psam}'
    plink2 \
      --pfile '${inputPrefix}' \
      --score ${weights} ${id_col} ${allele_col} ${effect_col} header center list-variants no-mean-imputation \
      --out ${trait} \
      --threads ${task.cpus} \
      --memory ${memMb}
    requested=${'$'}(awk 'NR > 1 { n++ } END { print n + 0 }' ${weights})
    matched=${'$'}(wc -l < ${trait}.sscore.vars)
    fraction=${'$'}(awk -v m="${'$'}matched" -v r="${'$'}requested" 'BEGIN { if (r == 0) print 0; else printf "%.8f", m / r }')
    awk -v trait='${trait}' -v requested="${'$'}requested" -v matched="${'$'}matched" \
      -v fraction="${'$'}fraction" -v warn='${params.warn_score_variant_match}' \
      -f ${summary_script} ${trait}.sscore > ${trait}.score_qc.tsv
    awk -v f="${'$'}fraction" -v minimum='${params.min_score_variant_match}' \
      'BEGIN { if (f < minimum) { printf "ERROR: score variant match fraction %.4f is below minimum %.4f\\n", f, minimum > "/dev/stderr"; exit 1 } }'
    """
}

process COLLATE_SCORE_RESULTS {
    label 'small'
    publishDir "${params.outdir}/05_scores", mode: 'copy'

    input:
    path sscores
    path score_qcs
    path collate_script

    output:
    path 'combined_scores.tsv', emit: combined
    path 'score_qc_summary.tsv', emit: qc

    script:
    """
    bash ${collate_script} \
      combined_scores.tsv score_qc_summary.tsv \
      ${sscores.join(' ')} -- ${score_qcs.join(' ')}
    """
}

process BUILD_ANALYSIS_DATASET {
    label 'small'
    container params.python_container
    publishDir "${params.outdir}/07_analysis", mode: 'copy'

    input:
    path combined_scores
    path global_pcs
    path ancestry_assignments
    path within_ancestry
    path builder_script

    output:
    path 'analysis_dataset.tsv', emit: dataset
    path 'analysis_dataset_dictionary.tsv', emit: dictionary

    script:
    """
    python ${builder_script} \
      --scores ${combined_scores} \
      --global-pcs ${global_pcs} \
      --ancestry ${ancestry_assignments} \
      --within-dir ${within_ancestry} \
      --num-global-pcs ${params.num_pcs} \
      --num-within-pcs ${params.num_within_ancestry_pcs} \
      --output analysis_dataset.tsv \
      --dictionary analysis_dataset_dictionary.tsv
    """
}

process HARMONIZE_PCA_PANEL {
    tag params.cohort
    label 'small'
    container params.python_container
    publishDir "${params.outdir}/06_pca/harmonization", mode: 'copy', pattern: 'panel_overlap*.tsv'

    input:
    tuple path(pgen), path(pvar), path(psam)
    path panel
    path harmonizer

    output:
    tuple path(pgen), path(pvar), path(psam), path('usable_cohort_ids.txt'),
      path('rename_to_panel_ids.tsv'), path('reference_alleles.tsv'), emit: harmonized_inputs
    path 'panel_overlap.tsv', emit: overlap
    path 'panel_overlap_summary.tsv', emit: summary

    script:
    """
    python ${harmonizer} \
      --panel ${panel} \
      --cohort-pvar ${pvar} \
      --min-overlap ${params.min_pca_variant_overlap} \
      --output-dir .
    """
}

process HARMONIZE_PCA_PANEL_DIRECT {
    tag params.cohort
    label 'small'
    container params.python_container
    publishDir "${params.outdir}/06_pca/harmonization", mode: 'copy', pattern: 'panel_overlap*.tsv'

    input:
    tuple val(pgen), val(pvar), val(psam)
    path panel
    path harmonizer

    output:
    tuple val(pgen), val(pvar), val(psam), path('usable_cohort_ids.txt'),
      path('rename_to_panel_ids.tsv'), path('reference_alleles.tsv'), emit: harmonized_inputs
    path 'panel_overlap.tsv', emit: overlap
    path 'panel_overlap_summary.tsv', emit: summary

    script:
    """
    test -r '${pgen}' && test -r '${pvar}' && test -r '${psam}'
    python ${harmonizer} \
      --panel ${panel} \
      --cohort-pvar '${pvar}' \
      --min-overlap ${params.min_pca_variant_overlap} \
      --output-dir .
    """
}

process VALIDATE_PCA_REFERENCE {
    tag reference_id
    label 'small'
    container params.python_container
    publishDir "${params.outdir}/06_pca/provenance", mode: 'copy', pattern: 'reference_validation.tsv'

    input:
    tuple val(reference_id), val(genome_build), path(panel), path(allele_frequencies),
      path(loadings), path(classifier), path(classifier_metadata), path(checksums)
    path validator

    output:
    tuple val(reference_id), val(genome_build), path(panel), path(allele_frequencies),
      path(loadings), path(classifier), path(classifier_metadata), path(checksums),
      path('reference_validation.tsv'), emit: reference

    script:
    """
    python ${validator} \
      --checksums ${checksums} \
      --artifact ${panel} \
      --artifact ${allele_frequencies} \
      --artifact ${loadings} \
      --artifact ${classifier} \
      --artifact ${classifier_metadata} \
      --output reference_validation.tsv
    """
}

process PREPARE_PCA_PFILE {
    tag params.cohort
    label 'large'

    input:
    tuple path(pgen), path(pvar), path(psam), path(usable_ids), path(rename_map), path(reference_alleles)

    output:
    tuple path('pca_input.pgen'), path('pca_input.pvar'), path('pca_input.psam'), emit: pfile

    script:
    def inputPrefix = pgen.baseName
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    plink2 \
      --pfile ${inputPrefix} \
      --extract ${usable_ids} \
      --update-name ${rename_map} 1 2 \
      --make-pgen \
      --out renamed \
      --threads ${task.cpus} \
      --memory ${memMb}

    plink2 \
      --pfile renamed \
      --ref-allele ${reference_alleles} 2 1 \
      --make-pgen \
      --out pca_input \
      --threads ${task.cpus} \
      --memory ${memMb}
    """
}

process PREPARE_PCA_PFILE_DIRECT {
    tag params.cohort
    label 'large'

    input:
    tuple val(pgen), val(pvar), val(psam), path(usable_ids), path(rename_map), path(reference_alleles)

    output:
    tuple path('pca_input.pgen'), path('pca_input.pvar'), path('pca_input.psam'), emit: pfile

    script:
    def inputPrefix = pgen.toString().replaceFirst(/[.]pgen$/, '')
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    test -r '${pgen}' && test -r '${pvar}' && test -r '${psam}'
    plink2 \
      --pfile '${inputPrefix}' \
      --extract ${usable_ids} \
      --update-name ${rename_map} 1 2 \
      --make-pgen \
      --out renamed \
      --threads ${task.cpus} \
      --memory ${memMb}

    plink2 \
      --pfile renamed \
      --ref-allele ${reference_alleles} 2 1 \
      --make-pgen \
      --out pca_input \
      --threads ${task.cpus} \
      --memory ${memMb}
    """
}

process PROJECT_GLOBAL_PCS {
    tag params.cohort
    label 'large'
    publishDir "${params.outdir}/06_pca/global", mode: 'copy'

    input:
    tuple path(pgen), path(pvar), path(psam)
    path allele_frequencies
    path loadings

    output:
    path 'global_pcs.tsv', emit: scores
    path 'projection_variants.txt', emit: variants

    script:
    def inputPrefix = pgen.baseName
    def lastPc = 5 + params.num_pcs.toInteger()
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    plink2 \
      --pfile ${inputPrefix} \
      --read-freq ${allele_frequencies} \
      --score ${loadings} 2 5 header-read no-mean-imputation variance-standardize list-variants \
      --score-col-nums 6-${lastPc} \
      --out projection \
      --threads ${task.cpus} \
      --memory ${memMb}
    cp projection.sscore global_pcs.tsv
    cp projection.sscore.vars projection_variants.txt
    """
}

process CLASSIFY_ANCESTRY {
    tag params.cohort
    label 'small'
    container params.python_container
    publishDir "${params.outdir}/06_pca/ancestry", mode: 'copy'

    input:
    path global_pcs
    path classifier
    path classifier_metadata
    path classifier_script

    output:
    path 'ancestry_probabilities.tsv', emit: assignments

    script:
    """
    python ${classifier_script} \
      --scores ${global_pcs} \
      --model ${classifier} \
      --classifier-metadata ${classifier_metadata} \
      --output ancestry_probabilities.tsv
    """
}

process WITHIN_ANCESTRY_PCA {
    tag params.cohort
    label 'large'
    publishDir "${params.outdir}/06_pca", mode: 'copy'

    input:
    tuple path(pgen), path(pvar), path(psam)
    path ancestry_assignments
    path within_ancestry_script

    output:
    path 'within_ancestry', emit: groups

    script:
    def inputPrefix = pgen.baseName
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    bash ${within_ancestry_script} \
      ${inputPrefix} \
      ${ancestry_assignments} \
      ${params.min_ancestry_samples} \
      ${params.num_within_ancestry_pcs} \
      ${params.within_ancestry_king_cutoff} \
      ${params.within_ancestry_ld_window} \
      ${params.within_ancestry_ld_step} \
      ${params.within_ancestry_ld_r2} \
      ${task.cpus} \
      ${memMb} \
      within_ancestry
    """
}

process WITHIN_ANCESTRY_PCA_DIRECT {
    tag params.cohort
    label 'large'
    publishDir "${params.outdir}/06_pca", mode: 'copy'

    input:
    tuple val(pgen), val(pvar), val(psam)
    path ancestry_assignments
    path within_ancestry_script

    output:
    path 'within_ancestry', emit: groups

    script:
    def inputPrefix = pgen.toString().replaceFirst(/[.]pgen$/, '')
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    test -r '${pgen}' && test -r '${pvar}' && test -r '${psam}'
    bash ${within_ancestry_script} \
      '${inputPrefix}' \
      ${ancestry_assignments} \
      ${params.min_ancestry_samples} \
      ${params.num_within_ancestry_pcs} \
      ${params.within_ancestry_king_cutoff} \
      ${params.within_ancestry_ld_window} \
      ${params.within_ancestry_ld_step} \
      ${params.within_ancestry_ld_r2} \
      ${task.cpus} \
      ${memMb} \
      within_ancestry
    """
}

workflow PGS_WORKFLOW {
    main:
    def skipRsid = flagEnabled(params.skip_rsid_annotation)
    def directInputsEnabled = flagEnabled(params.direct_inputs) || flagEnabled(params.fsx_direct)
    def directPfileEnabled = params.input_pfile && directInputsEnabled
    def pcaEnabled = flagEnabled(params.run_pca)
    def scoresEnabled = flagEnabled(params.run_scores)
    def summaryQcEnabled = flagEnabled(params.run_summary_qc)
    def globalPcResults = Channel.empty()
    def ancestryResults = Channel.empty()
    def withinAncestryResults = Channel.empty()
    def combinedScoreResults = Channel.empty()
    def analysisDatasetResults = Channel.empty()
    def analysisDictionaryResults = Channel.empty()

    if (pcaEnabled && !params.pca_reference_sheet) {
        error '--pca_reference_sheet is required when --run_pca is true'
    }

    if (!params.input_pfile && (!params.vcfs || (!params.rsid_maps && !skipRsid))) {
        error '--vcfs is required, and --rsid_maps is required unless --skip_rsid_annotation is true; alternatively provide --input_pfile with a QCed PLINK 2 prefix.'
    }

    if (params.input_pfile) {
        def pfilePrefix = params.input_pfile.toString()
        if (directPfileEnabled) {
            qcPfile = Channel.value(tuple(
                "${pfilePrefix}.pgen",
                "${pfilePrefix}.pvar",
                "${pfilePrefix}.psam"
            ))
        } else {
            qcPfile = Channel.value(tuple(
                file("${pfilePrefix}.pgen", checkIfExists: true),
                file("${pfilePrefix}.pvar", checkIfExists: true),
                file("${pfilePrefix}.psam", checkIfExists: true)
            ))
        }
    } else {
        def chromosomes = chromosomeList(params.chromosomes)
        inputStrings = Channel.fromList(chromosomes).map { chr ->
            tuple(
                chr,
                params.vcfs.toString().replace('{chr}', chr.toString()),
                params.rsid_maps ? params.rsid_maps.toString().replace('{chr}', chr.toString()) : ''
            )
        }

        if (directInputsEnabled && skipRsid) {
            directInputs = inputStrings.map { chr, vcf, map -> tuple(chr, vcf) }
            PREPROCESS_CHROMOSOME_DIRECT_NO_RSID(directInputs)
            chromosomePfiles = PREPROCESS_CHROMOSOME_DIRECT_NO_RSID.out.pfiles
        } else if (directInputsEnabled) {
            PREPROCESS_CHROMOSOME_DIRECT(inputStrings)
            chromosomePfiles = PREPROCESS_CHROMOSOME_DIRECT.out.pfiles
        } else {
            stagedInputs = inputStrings.map { chr, vcf, map ->
                tuple(chr, file(vcf, checkIfExists: true), file(map, checkIfExists: true))
            }
            PREPROCESS_CHROMOSOME(stagedInputs)
            chromosomePfiles = PREPROCESS_CHROMOSOME.out.pfiles
        }

        allChromosomeFiles = chromosomePfiles
            .map { chr, pgen, pvar, psam -> [pgen, pvar, psam] }
            .flatten()
            .collect()
        CONCAT_CHROMOSOMES(allChromosomeFiles)
        MISSINGNESS_QC(CONCAT_CHROMOSOMES.out.pfile)
        qcPfile = MISSINGNESS_QC.out.pfile
    }

    if (summaryQcEnabled) {
        if (directPfileEnabled) {
            SUMMARY_QC_DIRECT(qcPfile)
        } else {
            SUMMARY_QC(qcPfile)
        }
    }

    if (pcaEnabled) {
        pcaReference = Channel.fromPath(params.pca_reference_sheet, checkIfExists: true)
            .splitCsv(header: true, sep: '\t', strip: true)
            .map { row ->
                def required = ['reference_id', 'genome_build', 'panel_variants', 'allele_frequencies',
                    'loadings', 'classifier', 'classifier_metadata', 'checksums']
                def missing = required.findAll { !row[it] }
                if (missing) {
                    error "PCA reference sheet missing values: ${missing.join(', ')}"
                }
                if (row.genome_build.toString() != params.genome_build.toString()) {
                    error "Cohort genome build ${params.genome_build} does not match PCA reference build ${row.genome_build}"
                }
                tuple(
                    row.reference_id.toString(),
                    row.genome_build.toString(),
                    file(row.panel_variants.toString(), checkIfExists: true),
                    file(row.allele_frequencies.toString(), checkIfExists: true),
                    file(row.loadings.toString(), checkIfExists: true),
                    file(row.classifier.toString(), checkIfExists: true),
                    file(row.classifier_metadata.toString(), checkIfExists: true),
                    file(row.checksums.toString(), checkIfExists: true)
                )
            }
            .collect(flat: false)
            .map { rows ->
                if (rows.size() != 1) {
                    error 'PCA reference sheet must contain exactly one data row'
                }
                rows[0]
            }

        VALIDATE_PCA_REFERENCE(
            pcaReference,
            Channel.value(file("${projectDir}/bin/validate_pca_reference.py"))
        )
        validatedReference = VALIDATE_PCA_REFERENCE.out.reference
        pcaPanelCh = validatedReference.map { referenceId, build, panelFile, refFrequencies, refLoadings, refClassifier, metadata, checksums, validation -> panelFile }
        pcaFrequenciesCh = validatedReference.map { referenceId, build, panelFile, refFrequencies, refLoadings, refClassifier, metadata, checksums, validation -> refFrequencies }
        pcaLoadingsCh = validatedReference.map { referenceId, build, panelFile, refFrequencies, refLoadings, refClassifier, metadata, checksums, validation -> refLoadings }
        pcaClassifierCh = validatedReference.map { referenceId, build, panelFile, refFrequencies, refLoadings, refClassifier, metadata, checksums, validation -> refClassifier }
        pcaClassifierMetadataCh = validatedReference.map { referenceId, build, panelFile, refFrequencies, refLoadings, refClassifier, metadata, checksums, validation -> metadata }

        if (directPfileEnabled) {
            HARMONIZE_PCA_PANEL_DIRECT(qcPfile, pcaPanelCh, Channel.value(file("${projectDir}/bin/harmonize_pca_panel.py")))
            PREPARE_PCA_PFILE_DIRECT(HARMONIZE_PCA_PANEL_DIRECT.out.harmonized_inputs)
            pcaInput = PREPARE_PCA_PFILE_DIRECT.out.pfile
        } else {
            HARMONIZE_PCA_PANEL(qcPfile, pcaPanelCh, Channel.value(file("${projectDir}/bin/harmonize_pca_panel.py")))
            PREPARE_PCA_PFILE(HARMONIZE_PCA_PANEL.out.harmonized_inputs)
            pcaInput = PREPARE_PCA_PFILE.out.pfile
        }
        PROJECT_GLOBAL_PCS(pcaInput, pcaFrequenciesCh, pcaLoadingsCh)
        globalPcResults = PROJECT_GLOBAL_PCS.out.scores
        CLASSIFY_ANCESTRY(
            PROJECT_GLOBAL_PCS.out.scores,
            pcaClassifierCh,
            pcaClassifierMetadataCh,
            Channel.value(file("${projectDir}/bin/apply_extra_trees.py"))
        )
        ancestryResults = CLASSIFY_ANCESTRY.out.assignments
        if (directPfileEnabled) {
            WITHIN_ANCESTRY_PCA_DIRECT(
                qcPfile,
                CLASSIFY_ANCESTRY.out.assignments,
                Channel.value(file("${projectDir}/bin/run_within_ancestry_pca.sh"))
            )
            withinAncestryResults = WITHIN_ANCESTRY_PCA_DIRECT.out.groups
        } else {
            WITHIN_ANCESTRY_PCA(
                qcPfile,
                CLASSIFY_ANCESTRY.out.assignments,
                Channel.value(file("${projectDir}/bin/run_within_ancestry_pca.sh"))
            )
            withinAncestryResults = WITHIN_ANCESTRY_PCA.out.groups
        }
    }

    if (scoresEnabled) {
        if (!params.score_sheet) {
            error '--score_sheet is required when --run_scores is true'
        }
        weights = Channel.fromPath(params.score_sheet, checkIfExists: true)
            .splitCsv(header: true, sep: '\t', strip: true)
            .map { row ->
                if (!row.trait || !row.weights || !row.id_col || !row.allele_col || !row.effect_col) {
                    error 'Score sheet requires trait, weights, id_col, allele_col, and effect_col columns'
                }
                tuple(
                    row.trait.toString(),
                    file(row.weights.toString(), checkIfExists: true),
                    row.id_col.toString().toInteger(),
                    row.allele_col.toString().toInteger(),
                    row.effect_col.toString().toInteger()
                )
            }
        scoreInputs = weights.combine(qcPfile).map { trait, weight, idCol, alleleCol, effectCol, pgen, pvar, psam ->
            tuple(trait, weight, idCol, alleleCol, effectCol, pgen, pvar, psam)
        }
        if (directPfileEnabled) {
            SCORE_TRAIT_DIRECT(scoreInputs, Channel.value(file("${projectDir}/bin/summarize_score.awk")))
            scoredResults = SCORE_TRAIT_DIRECT.out.scored
        } else {
            SCORE_TRAIT(scoreInputs, Channel.value(file("${projectDir}/bin/summarize_score.awk")))
            scoredResults = SCORE_TRAIT.out.scored
        }
        scoreFiles = scoredResults.map { trait, score, qc -> score }.collect()
        scoreQcFiles = scoredResults.map { trait, score, qc -> qc }.collect()
        COLLATE_SCORE_RESULTS(
            scoreFiles,
            scoreQcFiles,
            Channel.value(file("${projectDir}/bin/collate_scores.sh"))
        )
        combinedScoreResults = COLLATE_SCORE_RESULTS.out.combined
    }

    if (pcaEnabled && scoresEnabled) {
        BUILD_ANALYSIS_DATASET(
            combinedScoreResults,
            globalPcResults,
            ancestryResults,
            withinAncestryResults,
            Channel.value(file("${projectDir}/bin/build_analysis_dataset.py"))
        )
        analysisDatasetResults = BUILD_ANALYSIS_DATASET.out.dataset
        analysisDictionaryResults = BUILD_ANALYSIS_DATASET.out.dictionary
    }

    emit:
    qc_pfile = qcPfile
    combined_scores = combinedScoreResults
    global_pcs = globalPcResults
    ancestry_assignments = ancestryResults
    within_ancestry = withinAncestryResults
    analysis_dataset = analysisDatasetResults
    analysis_dictionary = analysisDictionaryResults
}
