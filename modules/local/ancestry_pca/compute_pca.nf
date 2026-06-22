process ANCESTRY_PCA_COMPUTE_PCA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.12--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val   n_components

    output:
    tuple val(meta), path("plink2.eigenvec"), emit: eigenvec
    tuple val(meta), path("plink2.eigenval"), emit: eigenval
    path  "compute_pca.log",                  emit: log
    path  "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bfile = bed.baseName
    """
    plink2 \\
        --bfile ${bfile} \\
        --pca ${n_components} \\
        --out plink2 \\
        --allow-extra-chr \\
        --threads ${task.cpus}

    # Log variance explained for first 5 PCs (or fewer if n_components < 5)
    echo "Variance explained (eigenvalues):" > compute_pca.log
    N_REPORT=\$(( ${n_components} < 5 ? ${n_components} : 5 ))
    awk -v n="\${N_REPORT}" 'NR<=n {printf "  PC%d: %s\\n", NR, \$1}' plink2.eigenval >> compute_pca.log
    echo "" >> compute_pca.log
    echo "Full eigenval file: plink2.eigenval" >> compute_pca.log
    cat plink2.log >> compute_pca.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -n1 | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    printf '#FID\\t#IID\\tPC1\\tPC2\\n' > plink2.eigenvec
    printf '1.0\\n0.5\\n' > plink2.eigenval
    touch compute_pca.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -n1 | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
