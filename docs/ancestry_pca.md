# Ancestry PCA

## Overview

The ancestry PCA module performs principal component analysis (PCA) on germline variants to provide a quick view of population structure and detect sample outliers. It operates on the joint germline VCF produced by multi-sample HaplotypeCaller calling (`--joint_germline`).

The pipeline:
1. Filters VCF to PASS sites, autosomal chromosomes, minimum MAF, and maximum missingness
2. Converts the filtered VCF to PLINK binary format
3. Performs LD pruning to obtain an independent set of variants
4. Computes PCA using PLINK2
5. Formats results to TSV and produces a summary log

---

## Requirements

- A multi-sample germline joint VCF (produced with `--joint_germline`)
- The `--run_ancestry` flag must be set
- BCFtools ≥ 1.15 and PLINK2 ≥ 2.00a (provided via container/conda)

---

## Parameters

| Parameter            | Type    | Default | Description                                                                |
|----------------------|---------|---------|----------------------------------------------------------------------------|
| `run_ancestry`       | boolean | `false` | Enable ancestry PCA analysis on joint germline VCF                         |
| `pca_n_components`   | integer | `10`    | Number of principal components to compute (typically 3–20)                 |
| `pca_maf`            | number  | `0.05`  | Minor allele frequency threshold for PCA variant filtering                 |
| `pca_max_missing`    | number  | `0.01`  | Maximum genotype missingness per site for PCA                              |
| `ld_window_kb`       | integer | `50`    | LD pruning window size in kilobases                                        |
| `ld_step_kb`         | integer | `10`    | LD pruning step size in kilobases                                          |
| `ld_r2`              | number  | `0.1`   | LD pruning r² threshold (strict — use 0.2–0.5 for less stringent pruning) |
| `plot_pca`           | boolean | `false` | Generate PC1 vs PC2 scatter plot (requires R/ggplot2 — reserved for future use) |

---

## Example Command

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --tools haplotypecaller \
    --joint_germline \
    --run_ancestry \
    --genome GRCh38 \
    --outdir results
```

To override default PCA parameters:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --tools haplotypecaller \
    --joint_germline \
    --run_ancestry \
    --pca_maf 0.01 \
    --pca_max_missing 0.05 \
    --pca_n_components 20 \
    --ld_r2 0.2 \
    --genome GRCh38 \
    --outdir results
```

---

## Output Files

All files are written to `<outdir>/ancestry_pca/`.

| File                       | Description                                                            |
|----------------------------|------------------------------------------------------------------------|
| `pca.tsv`                  | PCA coordinates per sample: `sample_id`, `PC1` … `PCn`                |
| `pca_variance.tsv`         | Variance explained per PC: `PC`, `variance_explained_percent`, `cumulative_percent` |
| `pca.log`                  | Summary log: sample count, component count, variance table             |
| `plink2.eigenvec`          | Raw PLINK2 eigenvector file                                            |
| `plink2.eigenval`          | Raw PLINK2 eigenvalue file                                             |
| `filter_variants.log`      | Site counts at each filtering step                                     |
| `convert_vcf_to_plink.log` | Sample and variant counts after VCF → PLINK conversion                 |
| `ld_prune_variants.log`    | Variant counts before/after LD pruning                                 |
| `compute_pca.log`          | Eigenvalue summary and PLINK2 run log                                  |

---

## Assumptions and Notes

- **Germline only**: ancestry PCA is computed from germline variant calls, not somatic calls.
- **Autosomes only**: sex chromosomes are excluded to avoid confounding by sex.
- **PASS sites only**: sites not passing variant calling filters are excluded.
- **MAF rationale**: a default MAF of 5% excludes rare variants that contribute noise rather than population signal. Lower this threshold for small cohorts.
- **LD pruning**: the default r² threshold of 0.1 is stricter than commonly used values (0.2–0.5). This improves PCA accuracy at the cost of retaining fewer variants.
- **Chromosome naming**: both UCSC (`chr1`–`chr22`) and Ensembl (`1`–`22`) naming conventions are supported; the filtering step auto-detects the convention from the VCF header.
- **Multi-sample requirement**: PCA is only meaningful with multiple samples. A warning is emitted if `--joint_germline` is not set.
