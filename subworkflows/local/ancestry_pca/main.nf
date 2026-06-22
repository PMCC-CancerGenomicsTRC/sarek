include { ANCESTRY_PCA_FILTER_VARIANTS      } from '../../../modules/local/ancestry_pca/filter_variants'
include { ANCESTRY_PCA_CONVERT_VCF_TO_PLINK } from '../../../modules/local/ancestry_pca/convert_vcf_to_plink'
include { ANCESTRY_PCA_LD_PRUNE_VARIANTS    } from '../../../modules/local/ancestry_pca/ld_prune_variants'
include { ANCESTRY_PCA_COMPUTE_PCA          } from '../../../modules/local/ancestry_pca/compute_pca'
include { ANCESTRY_PCA_FORMAT_OUTPUT        } from '../../../modules/local/ancestry_pca/format_pca_output'

workflow ANCESTRY_PCA {

    take:
    vcf          // channel: [ meta, vcf.gz, tbi ]
    maf
    max_missing
    ld_window_kb
    ld_step_kb
    ld_r2
    n_components

    main:

    versions = Channel.empty()

    ANCESTRY_PCA_FILTER_VARIANTS(vcf, maf, max_missing)
    versions = versions.mix(ANCESTRY_PCA_FILTER_VARIANTS.out.versions)

    ANCESTRY_PCA_CONVERT_VCF_TO_PLINK(ANCESTRY_PCA_FILTER_VARIANTS.out.vcf)
    versions = versions.mix(ANCESTRY_PCA_CONVERT_VCF_TO_PLINK.out.versions)

    ANCESTRY_PCA_LD_PRUNE_VARIANTS(ANCESTRY_PCA_CONVERT_VCF_TO_PLINK.out.plink, ld_window_kb, ld_step_kb, ld_r2)
    versions = versions.mix(ANCESTRY_PCA_LD_PRUNE_VARIANTS.out.versions)

    ANCESTRY_PCA_COMPUTE_PCA(ANCESTRY_PCA_LD_PRUNE_VARIANTS.out.plink, n_components)
    versions = versions.mix(ANCESTRY_PCA_COMPUTE_PCA.out.versions)

    ANCESTRY_PCA_FORMAT_OUTPUT(ANCESTRY_PCA_COMPUTE_PCA.out.eigenvec, ANCESTRY_PCA_COMPUTE_PCA.out.eigenval, n_components)

    emit:
    pca_tsv      = ANCESTRY_PCA_FORMAT_OUTPUT.out.pca_tsv
    pca_variance = ANCESTRY_PCA_FORMAT_OUTPUT.out.pca_variance
    pca_log      = ANCESTRY_PCA_FORMAT_OUTPUT.out.pca_log
    eigenvec     = ANCESTRY_PCA_COMPUTE_PCA.out.eigenvec
    eigenval     = ANCESTRY_PCA_COMPUTE_PCA.out.eigenval
    versions
}
