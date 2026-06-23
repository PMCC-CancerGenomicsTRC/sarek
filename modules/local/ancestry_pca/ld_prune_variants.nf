process ANCESTRY_PCA_LD_PRUNE_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.12--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val   window_kb
    val   step_kb
    val   r2

    output:
    tuple val(meta), path("pruned_input.bed"), path("pruned_input.bim"), path("pruned_input.fam"), emit: plink
    path  "pruned.prune.in",                                                                       emit: prune_in
    path  "ld_prune_variants.log",                                                                 emit: log
    path  "versions.yml",                                                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bfile = bed.baseName
    """
    # Count variants before pruning
    N_BEFORE=\$(wc -l < ${bim})
    echo "Variants before LD pruning: \${N_BEFORE}" > ld_prune_variants.log

    # Step 3a: Generate prune list
    plink2 \\
        --bfile ${bfile} \\
        --indep-pairwise ${window_kb}kb ${step_kb} ${r2} \\
        --out pruned \\
        --allow-extra-chr \\
        --threads ${task.cpus}

    # Step 3b: Extract pruned variants
    plink2 \\
        --bfile ${bfile} \\
        --extract pruned.prune.in \\
        --make-bed \\
        --out pruned_input \\
        --allow-extra-chr \\
        --threads ${task.cpus}

    N_AFTER=\$(wc -l < pruned_input.bim)
    echo "Variants after LD pruning: \${N_AFTER}" >> ld_prune_variants.log
    echo "Pruning parameters: window=${window_kb}kb step=${step_kb} r2=${r2}" >> ld_prune_variants.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -n1 | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch pruned_input.bed pruned_input.bim pruned_input.fam
    touch pruned.prune.in
    touch ld_prune_variants.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -n1 | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
