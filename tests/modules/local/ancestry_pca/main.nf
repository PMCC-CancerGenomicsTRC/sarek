nextflow_process {

    name "Test ANCESTRY_PCA_FILTER_VARIANTS"
    script "modules/local/ancestry_pca/filter_variants.nf"
    process "ANCESTRY_PCA_FILTER_VARIANTS"

    test("Should filter VCF and produce filtered.vcf.gz") {
        when {
            params {
                pca_maf         = 0.05
                pca_max_missing = 0.01
            }
            process {
                """
                input[0] = [ [id:'test'], file(params.test_data['vcf_gz']), file(params.test_data['vcf_gz_tbi']) ]
                input[1] = 0.05
                input[2] = 0.01
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.vcf.get(0).get(1)).exists()
            assert path(process.out.log.get(0)).text.contains("After MAF")
        }
    }
}

nextflow_process {

    name "Test ANCESTRY_PCA_CONVERT_VCF_TO_PLINK"
    script "modules/local/ancestry_pca/convert_vcf_to_plink.nf"
    process "ANCESTRY_PCA_CONVERT_VCF_TO_PLINK"

    test("Should convert VCF to PLINK binary format") {
        when {
            process {
                """
                input[0] = [ [id:'test'], file(params.test_data['vcf_gz']), file(params.test_data['vcf_gz_tbi']) ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.plink.get(0).get(1)).exists()  // .bed
            assert path(process.out.plink.get(0).get(2)).exists()  // .bim
            assert path(process.out.plink.get(0).get(3)).exists()  // .fam
        }
    }
}

nextflow_process {

    name "Test ANCESTRY_PCA_LD_PRUNE_VARIANTS"
    script "modules/local/ancestry_pca/ld_prune_variants.nf"
    process "ANCESTRY_PCA_LD_PRUNE_VARIANTS"

    test("Should produce LD-pruned PLINK files") {
        when {
            process {
                """
                input[0] = [ [id:'test'], file(params.test_data['plink_bed']), file(params.test_data['plink_bim']), file(params.test_data['plink_fam']) ]
                input[1] = 50
                input[2] = 10
                input[3] = 0.1
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.plink.get(0).get(1)).exists()  // pruned_input.bed
            assert path(process.out.prune_in.get(0)).exists()
        }
    }
}

nextflow_process {

    name "Test ANCESTRY_PCA_COMPUTE_PCA"
    script "modules/local/ancestry_pca/compute_pca.nf"
    process "ANCESTRY_PCA_COMPUTE_PCA"

    test("Should compute PCA and produce eigenvec and eigenval files") {
        when {
            process {
                """
                input[0] = [ [id:'test'], file(params.test_data['plink_bed']), file(params.test_data['plink_bim']), file(params.test_data['plink_fam']) ]
                input[1] = 10
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.eigenvec.get(0).get(1)).exists()
            assert path(process.out.eigenval.get(0).get(1)).exists()
        }
    }
}

nextflow_process {

    name "Test ANCESTRY_PCA_FORMAT_OUTPUT"
    script "modules/local/ancestry_pca/format_pca_output.nf"
    process "ANCESTRY_PCA_FORMAT_OUTPUT"

    test("Should produce pca.tsv, pca_variance.tsv and pca.log") {
        when {
            process {
                """
                input[0] = [ [id:'test'], file(params.test_data['eigenvec']) ]
                input[1] = [ [id:'test'], file(params.test_data['eigenval']) ]
                input[2] = 10
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.pca_tsv.get(0).get(1)).exists()
            assert path(process.out.pca_variance.get(0).get(1)).exists()
            assert path(process.out.pca_tsv.get(0).get(1)).text.contains("sample_id")
            assert !path(process.out.pca_tsv.get(0).get(1)).text.contains("#FID")
            assert path(process.out.pca_variance.get(0).get(1)).text.contains("PC")
        }
    }
}
