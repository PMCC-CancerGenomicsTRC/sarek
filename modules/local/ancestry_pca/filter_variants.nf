process ANCESTRY_PCA_FILTER_VARIANTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val   maf
    val   max_missing

    output:
    tuple val(meta), path("filtered.vcf.gz"), path("filtered.vcf.gz.tbi"), emit: vcf
    path  "filter_variants.log",                                            emit: log
    path  "versions.yml",                                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Build autosome region list — handles both UCSC (chr1-chr22) and Ensembl (1-22) naming
    """
    # Detect chromosome naming convention from VCF header
    if bcftools view -h ${vcf} | grep -q '^##contig.*ID=chr'; then
        CHROMS=\$(seq 1 22 | awk '{printf "chr%s,", \$1}' | sed 's/,\$//')
    else
        CHROMS=\$(seq 1 22 | tr '\\n' ',' | sed 's/,\$//')
    fi

    # Count sites before filtering
    BEFORE=\$(bcftools view --no-header ${vcf} | wc -l)
    echo "Sites before filtering: \${BEFORE}" > filter_variants.log

    # Step 1: Filter to PASS sites only
    bcftools view -f PASS ${vcf} -O z -o pass_only.vcf.gz
    bcftools index --tbi pass_only.vcf.gz
    AFTER_PASS=\$(bcftools view --no-header pass_only.vcf.gz | wc -l)
    echo "After PASS filter: \${AFTER_PASS}" >> filter_variants.log

    # Step 2: Apply MAF threshold
    bcftools view --min-af ${maf}:minor pass_only.vcf.gz -O z -o maf_filtered.vcf.gz
    bcftools index --tbi maf_filtered.vcf.gz
    AFTER_MAF=\$(bcftools view --no-header maf_filtered.vcf.gz | wc -l)
    echo "After MAF >= ${maf} filter: \${AFTER_MAF}" >> filter_variants.log

    # Step 3: Apply site missingness threshold
    bcftools view --include "F_MISSING < ${max_missing}" maf_filtered.vcf.gz -O z -o missing_filtered.vcf.gz
    bcftools index --tbi missing_filtered.vcf.gz
    AFTER_MISSING=\$(bcftools view --no-header missing_filtered.vcf.gz | wc -l)
    echo "After missingness < ${max_missing} filter: \${AFTER_MISSING}" >> filter_variants.log

    # Step 4: Restrict to autosomes
    bcftools view --regions "\${CHROMS}" missing_filtered.vcf.gz -O z -o filtered.vcf.gz
    bcftools index --tbi filtered.vcf.gz
    AFTER_AUTOSOMES=\$(bcftools view --no-header filtered.vcf.gz | wc -l)
    echo "After autosome filter: \${AFTER_AUTOSOMES}" >> filter_variants.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        htslib: \$(bcftools --version 2>&1 | grep 'htslib' | sed 's/^.*htslib //')
    END_VERSIONS
    """

    stub:
    """
    echo "stub" | gzip > filtered.vcf.gz
    touch filtered.vcf.gz.tbi
    touch filter_variants.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        htslib: \$(bcftools --version 2>&1 | grep 'htslib' | sed 's/^.*htslib //')
    END_VERSIONS
    """
}
