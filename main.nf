nextflow.enable.dsl=2

params.vcfs = null
params.rsid_maps = null
params.chromosomes = '1..22'
params.cohort = 'cohort'
params.genome_build = 'GRCh38'
params.outdir = 'results'
params.score_sheet = null
params.run_scores = false
params.run_pca = false
params.pca_reference_sheet = null
params.r2 = null
params.aq = null
params.mac = 10
params.geno = 0.05
params.maf = 0.01
params.sample_miss = 0.05
params.variant_miss = 0.05
params.num_pcs = 10
params.min_pca_variant_overlap = 0.90
params.min_ancestry_samples = 50
params.python_container = 'python:3.11-slim'
params.fsx_direct = false
params.direct_inputs = false
params.skip_rsid_annotation = false

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

process SCORE_TRAIT {
    tag trait
    label 'large'
    publishDir "${params.outdir}/05_scores", mode: 'copy'

    input:
    tuple val(trait), path(weights), val(id_col), val(allele_col), val(effect_col), path(pgen), path(pvar), path(psam)

    output:
    tuple val(trait), path("${trait}.sscore"), emit: scores
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

workflow {
    def skipRsid = flagEnabled(params.skip_rsid_annotation)
    def directInputsEnabled = flagEnabled(params.direct_inputs) || flagEnabled(params.fsx_direct)
    def pcaEnabled = flagEnabled(params.run_pca)
    def scoresEnabled = flagEnabled(params.run_scores)

    if (pcaEnabled && !params.pca_reference_sheet) {
        error '--pca_reference_sheet is required when --run_pca is true'
    }

    if (!params.vcfs || (!params.rsid_maps && !skipRsid)) {
        error '--vcfs is required, and --rsid_maps is required unless --skip_rsid_annotation is true.'
    }

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
    SUMMARY_QC(qcPfile)

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
            .collect()
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
        panel = validatedReference.map { referenceId, build, panelFile, frequencies, loadings, classifier, metadata, checksums, validation -> panelFile }
        frequencies = validatedReference.map { referenceId, build, panelFile, frequencies, loadings, classifier, metadata, checksums, validation -> frequencies }
        loadings = validatedReference.map { referenceId, build, panelFile, frequencies, loadings, classifier, metadata, checksums, validation -> loadings }
        classifier = validatedReference.map { referenceId, build, panelFile, frequencies, loadings, classifier, metadata, checksums, validation -> classifier }
        classifierMetadata = validatedReference.map { referenceId, build, panelFile, frequencies, loadings, classifier, metadata, checksums, validation -> metadata }

        HARMONIZE_PCA_PANEL(qcPfile, panel, Channel.value(file("${projectDir}/bin/harmonize_pca_panel.py")))
        PREPARE_PCA_PFILE(HARMONIZE_PCA_PANEL.out.harmonized_inputs)
        PROJECT_GLOBAL_PCS(PREPARE_PCA_PFILE.out.pfile, frequencies, loadings)
        CLASSIFY_ANCESTRY(
            PROJECT_GLOBAL_PCS.out.scores,
            classifier,
            classifierMetadata,
            Channel.value(file("${projectDir}/bin/apply_extra_trees.py"))
        )
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
        SCORE_TRAIT(scoreInputs)
    }
}
