nextflow.enable.dsl=2

params.vcfs = null
params.rsid_maps = null
params.chromosomes = '1..22'
params.cohort = 'cohort'
params.outdir = 'results'
params.score_sheet = null
params.run_scores = false
params.run_ancestry = true
params.r2 = null
params.aq = null
params.mac = 10
params.geno = 0.05
params.maf = 0.01
params.sample_miss = 0.05
params.variant_miss = 0.05
params.ancestry_maf = 0.05
params.king_cutoff = 0.0884
params.num_pcs = 10
params.ld_window = 200
params.ld_step = 50
params.ld_r2 = 0.1
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
    def memMb = Math.max(1000, (task.memory.toMega() as int) - 1000)
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
    def memMb = Math.max(1000, (task.memory.toMega() as int) - 1000)
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
    def memMb = Math.max(1000, (task.memory.toMega() as int) - 1000)
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
    def memMb = Math.max(1000, (task.memory.toMega() as int) - 2000)
    """
    find . -maxdepth 1 -name 'chr*.pgen' -print | sed 's/\.pgen${'$'}//' | sort -V > merge_list.txt
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
    def memMb = Math.max(1000, (task.memory.toMega() as int) - 2000)
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
    def memMb = Math.max(1000, (task.memory.toMega() as int) - 2000)
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
    tuple val(trait), path(weights), path(pgen), path(pvar), path(psam)

    output:
    tuple val(trait), path("${trait}.sscore"), emit: scores
    path "${trait}.sscore.vars"

    script:
    def inputPrefix = pgen.baseName
    def memMb = Math.max(1000, (task.memory.toMega() as int) - 2000)
    """
    plink2 \
      --pfile ${inputPrefix} \
      --score ${weights} 2 5 8 header center list-variants no-mean-imputation \
      --out ${trait} \
      --threads ${task.cpus} \
      --memory ${memMb}
    """
}

process ANCESTRY {
    tag params.cohort
    label 'large'
    publishDir "${params.outdir}/06_ancestry", mode: 'copy'

    input:
    tuple path(pgen), path(pvar), path(psam)

    output:
    path 'ancestry/*'

    script:
    def inputPrefix = pgen.baseName
    def memMb = Math.max(1000, (task.memory.toMega() as int) - 2000)
    def lastPc = 5 + params.num_pcs.toInteger()
    """
    mkdir ancestry
    plink2 --pfile ${inputPrefix} --maf ${params.ancestry_maf} --make-pgen \
      --out ancestry/maf --threads ${task.cpus} --memory ${memMb}
    plink2 --pfile ancestry/maf --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} \
      --out ancestry/prune --threads ${task.cpus} --memory ${memMb}
    plink2 --pfile ancestry/maf --extract ancestry/prune.prune.in --make-pgen \
      --out ancestry/pruned --threads ${task.cpus} --memory ${memMb}
    plink2 --pfile ancestry/pruned --make-king-table --out ancestry/king \
      --threads ${task.cpus} --memory ${memMb}
    plink2 --pfile ancestry/pruned --king-cutoff ${params.king_cutoff} --make-just-fam \
      --out ancestry/unrelateds --threads ${task.cpus} --memory ${memMb}
    plink2 --pfile ancestry/pruned --keep ancestry/unrelateds.king.cutoff.in.id \
      --freq counts --pca ${params.num_pcs} allele-wts vcols=chrom,ref,alt \
      --out ancestry/ref --threads ${task.cpus} --memory ${memMb}
    plink2 --pfile ancestry/pruned --read-freq ancestry/ref.acount \
      --score ancestry/ref.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
      --score-col-nums 6-${lastPc} --out ancestry/${params.cohort} \
      --threads ${task.cpus} --memory ${memMb}
    """
}

workflow {
    if (!params.vcfs || (!params.rsid_maps && !params.skip_rsid_annotation)) {
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

    def useDirectInputs = params.direct_inputs || params.fsx_direct
    if (useDirectInputs && params.skip_rsid_annotation) {
        directInputs = inputStrings.map { chr, vcf, map -> tuple(chr, vcf) }
        PREPROCESS_CHROMOSOME_DIRECT_NO_RSID(directInputs)
        chromosomePfiles = PREPROCESS_CHROMOSOME_DIRECT_NO_RSID.out.pfiles
    } else if (useDirectInputs) {
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

    if (params.run_ancestry) {
        ANCESTRY(qcPfile)
    }

    if (params.run_scores) {
        if (!params.score_sheet) {
            error '--score_sheet is required when --run_scores is true'
        }
        weights = Channel.fromPath(params.score_sheet, checkIfExists: true)
            .splitCsv(sep: '\t', strip: true)
            .filter { row -> row && row[0] && !row[0].toString().startsWith('#') }
            .map { row -> tuple(row[0].toString(), file(row[1].toString(), checkIfExists: true)) }
        scoreInputs = weights.combine(qcPfile).map { trait, weight, pgen, pvar, psam ->
            tuple(trait, weight, pgen, pvar, psam)
        }
        SCORE_TRAIT(scoreInputs)
    }
}
