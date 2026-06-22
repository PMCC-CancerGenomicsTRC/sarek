process ANCESTRY_PCA_CONVERT_VCF_TO_PLINK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.12--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("plink_input.bed"), path("plink_input.bim"), path("plink_input.fam"), emit: plink
    path  "convert_vcf_to_plink.log",                                                            emit: log
    path  "versions.yml",                                                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plink2 \\
        --vcf ${vcf} \\
        --make-bed \\
        --out plink_input \\
        --allow-extra-chr \\
        --chr 1-22 \\
        --threads ${task.cpus}

    # Log sample and variant counts from the PLINK log
    N_SAMPLES=\$(awk '/samples.*remaining after/ {print \$1}' plink_input.log || awk '/^[0-9]+ samples/ {print \$1}' plink_input.log || echo "see plink_input.log")
    N_VARIANTS=\$(awk '/variants.*remaining after/ {print \$1}' plink_input.log || awk '/^[0-9]+ variants/ {print \$1}' plink_input.log || echo "see plink_input.log")
    echo "Samples: \${N_SAMPLES}" > convert_vcf_to_plink.log
    echo "Variants: \${N_VARIANTS}" >> convert_vcf_to_plink.log
    cat plink_input.log >> convert_vcf_to_plink.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -n1 | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch plink_input.bed plink_input.bim plink_input.fam
    touch convert_vcf_to_plink.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -n1 | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
