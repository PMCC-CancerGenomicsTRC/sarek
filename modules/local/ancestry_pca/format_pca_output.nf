process ANCESTRY_PCA_FORMAT_OUTPUT {
    tag "$meta.id"
    label 'process_single'

    // Uses standard POSIX shell tools only — no extra container needed
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'quay.io/biocontainers/ubuntu:20.04' }"

    input:
    tuple val(meta), path(eigenvec)
    tuple val(meta2), path(eigenval)
    val   n_components

    output:
    tuple val(meta), path("pca.tsv"),          emit: pca_tsv
    tuple val(meta), path("pca_variance.tsv"), emit: pca_variance
    path  "pca.log",                           emit: pca_log
    path  "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Format eigenvec: drop #FID column, rename #IID to sample_id
    awk 'BEGIN{OFS="\\t"} NR==1 {
        # Replace "#FID" header marker; first col is #FID, second is #IID
        \$1=""
        sub(/^\\t/, "", \$0)
        # Replace "#IID" with "sample_id"
        sub(/#IID/, "sample_id")
        print
    } NR>1 {
        # Drop first column (FID)
        \$1=""
        sub(/^\\t/, "", \$0)
        print
    }' ${eigenvec} > pca.tsv

    # Build pca_variance.tsv with cumulative variance percent
    awk 'BEGIN {OFS="\\t"; print "PC", "variance_explained_percent", "cumulative_percent"; cum=0}
    {
        total[NR] = \$1
    }
    END {
        sum = 0
        for (i=1; i<=NR; i++) sum += total[i]
        cum = 0
        for (i=1; i<=NR; i++) {
            pct = (sum > 0) ? (total[i] / sum * 100) : 0
            cum += pct
            printf "PC%d\\t%.4f\\t%.4f\\n", i, pct, cum
        }
    }' ${eigenval} > pca_variance.tsv

    # Write summary log
    N_SAMPLES=\$(awk 'NR>1' pca.tsv | wc -l)
    N_PC=${n_components}
    echo "Ancestry PCA summary" > pca.log
    echo "====================" >> pca.log
    echo "Samples:          \${N_SAMPLES}" >> pca.log
    echo "Components:       \${N_PC}" >> pca.log
    echo "" >> pca.log
    echo "Variance explained:" >> pca.log
    cat pca_variance.tsv >> pca.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1)
    END_VERSIONS
    """

    stub:
    """
    printf 'sample_id\\tPC1\\tPC2\\n' > pca.tsv
    printf 'PC\\tvariance_explained_percent\\tcumulative_percent\\n' > pca_variance.tsv
    touch pca.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1)
    END_VERSIONS
    """
}
