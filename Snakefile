# Snakefile

# List of chromosomes
CHROMS = [str(i) for i in range(21, 23)]

rule all:
    input:
        expand(["03_maf/chr{chr}.pgen",
                "03_maf/chr{chr}.psam",
                "03_maf/chr{chr}.pvar"], chr=CHROMS)

rule FILTER_VCF:
    input:
        vcf = "vcfs/chr{chr}_jointcall_VQSR_combined.vcf.gz"
    output:
        pgen = "01_filter/chr{chr}.pgen",
        psam = "01_filter/chr{chr}.psam",
        pvar = "01_filter/chr{chr}.pvar"
    threads: 24
    resources:
        mem_mb = 16000
    shell:
        """
        plink2 \
          --vcf {input.vcf} \
          --snps-only just-acgt \
          --max-alleles 2 \
          --var-filter \
          --mac 10 \
          --geno 0.05 \
          --set-all-var-ids @:#:\$r:\$a \
          --make-pgen \
          --out 01_filter/chr{wildcards.chr} \
          --threads {threads} \
          --memory {resources.mem_mb}
        """

rule ANNOTATE_RSID:
    input:
        pgen = "01_filter/chr{chr}.pgen",
        psam = "01_filter/chr{chr}.psam",
        pvar = "01_filter/chr{chr}.pvar",
        map = "rsid_maps/chr{chr}.map"
    output:
        pgen = "02_annotate/chr{chr}.pgen",
        psam = "02_annotate/chr{chr}.psam",
        pvar = "02_annotate/chr{chr}.pvar"
    threads: 24
    resources:
        mem_mb = 16000
    shell:
        """
        plink2 \
          --pfile 01_filter/chr{wildcards.chr} \
          --update-name {input.map} 2 1 \
          --extract <(cut -f2 {input.map}) \
          --write-snplist \
          --make-pgen \
          --out 02_annotate/chr{wildcards.chr} \
          --threads {threads} \
          --memory {resources.mem_mb}
        """

rule MAF_FILTER:
    input:
        pgen = "02_annotate/chr{chr}.pgen",
        psam = "02_annotate/chr{chr}.psam",
        pvar = "02_annotate/chr{chr}.pvar"
    output:
        pgen = "03_maf/chr{chr}.pgen",
        psam = "03_maf/chr{chr}.psam",
        pvar = "03_maf/chr{chr}.pvar"
    threads: 24
    resources:
        mem_mb = 16000 
    shell:
        r"""
        plink2 \
          --pfile 02_annotate/chr{wildcards.chr} \
          --maf 0.01 \
          --make-pgen \
          --out 03_maf/chr{wildcards.chr} \
          --threads {threads} \
          --memory {resources.mem_mb}
        """

rule MISSINGNESS:
    input:
        expand("03_maf/chr{chr}.pgen", chr=CHROMS)
    output:
        expand("04_missingness/chr{chr}.pgen", chr=CHROMS)
		smiss = "04_missingness/missingness.smiss"		
    threads: 24
    resources:
        mem_mb = 16000
    shell:
        """
        MERGE_LIST=04_missingness/merge_list.txt
        > $MERGE_LIST

        echo "Building merge list..."
        for CHR in {{1..22}}; do
          PREFIX="03_maf/chr${{CHR}}"
          if [[ -f "${{PREFIX}}.pgen" ]]; then
            echo "${{PREFIX}}" >> $MERGE_LIST
          fi
        done

        echo "Running PLINK2 missingness..."
        plink2 \
          --pmerge-list $MERGE_LIST \
          --missing \
          --out 04_missingness/missingness \
          --threads {threads} \
          --memory {resources.mem_mb}

        echo "Generating fail lists..."
        awk 'NR>1 && $4>0.05 {{print $1}}' 04_missingness/missingness.smiss > 04_missingness/fail_samples.txt
        awk 'NR>1 && $5>0.05 {{print $2}}' 04_missingness/missingness.vmiss > 04_missingness/fail_variants.txt

        echo "Applying missingness filters..."
        for CHR in {{1..22}}; do
          PREFIX="03_maf/chr${{CHR}}"
          if [[ -f "${{PREFIX}}.pgen" ]]; then
            echo "[Chr${{CHR}}] filtering..."
            plink2 \
              --pfile ${{PREFIX}} \
              --remove 04_missingness/fail_samples.txt \
              --exclude 04_missingness/fail_variants.txt \
              --make-pgen \
              --out 04_missingness/chr${{CHR}} \
              --threads {threads} \
              --memory {resources.mem_mb}
          fi
        done
        """
