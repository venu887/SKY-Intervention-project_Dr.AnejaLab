# SKY Breathwork RCT — AU-Rich Element (ARE) Analysis


This repository contains all code and input data needed to reproduce the AU-rich element (ARE) post-transcriptional analysis figures and statistics in:

> Ganesan et al. (2026). *Coordinated immune-metabolic remodelling following SKY breathwork: an integrated DNA methylation and transcriptomic analysis.* Nature Communications, in press.

Specifically, this code reproduces:

- **Main paper Figure 6C** — ARE enrichment bar chart
- **Main paper Figure 6D** — RBP motif specificity forest plot
- **Supplementary Figure S1** — CDF of log₂FC by ZFP36 family 7-mer count
- **Supplementary Figure S2** — Partial Spearman scatter (UTR-length corrected)
- **Supplementary Tables S1, S2** — All ARE and RBP statistics


## Quick start

```bash

# Create a clean environment
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Run the analysis (~30 seconds)
python ARE_analysis.py

# Run the test suite (~10 seconds)
python test_ARE_analysis.py
```

Both scripts should print "ALL DONE" / "ALL TESTS PASSED — code is ready for submission".

## Output

After running `ARE_analysis.py`, the following files appear in the repo root:

| File | Contents |
|---|---|
| `Fig_6C_ARE_bar.pdf` | ARE enrichment bar chart (3 groups, OR/p annotated) |
| `Fig_6D_RBP_forest.pdf` | Forest plot of 5 RBP motif odds ratios |
| `Fig_S1_CDF_7mer.pdf` | Cumulative distribution by 7-mer count |
| `Fig_S2_partial_Spearman.pdf` | Scatter with primary pre-specified test |
| `ARE_statistics_verified.csv` | Every statistic in the manuscript, numerically verified |
| `RBP_statistics_verified.csv` | RBP panel with independent ZFP36L2 cross-check |

```

## Input data

### `data/02_ARE_annotations.csv`

Per-gene annotation file (n=9,732 protein-coding genes with annotated 3′ UTRs).

| Column | Description |
|---|---|
| `Gene` | HGNC gene symbol |
| `Intervention_logFC` | log₂ fold-change from ITT interaction-term limma model |
| `Intervention_P.Value` | nominal p-value from ITT model |
| `group` | one of `down` (n=109), `background` (n=9,588), `up` (n=35) |
| `utr_length` | length of annotated 3′ UTR in nucleotides |
| `AREScore` | composite ARE score (Spasic et al. *PLoS Genet* 2012; ARED-Plus database, Bakheet et al. 2018) |
| `ARED_cluster` | ARED-Plus cluster annotation (0=no ARE evidence; 1=any; 2=strongly destabilising) |
| `n_7mer_L2` | count of UAUUUAU heptamers in 3′ UTR (ZFP36 family recognition motif; Redmon et al. 2022) |

### `data/03_RBP_specificity.csv`

Precomputed RBP motif specificity statistics for 5 RNA-binding proteins.

| Column | Description |
|---|---|
| `RBP` | RNA-binding protein name |
| `Role` | `destabiliser` or `stabiliser` |
| `DOWN_OR` | Fisher's exact OR (downregulated vs background) |
| `DOWN_CI_low`, `DOWN_CI_high` | 95% Woolf logit confidence interval |
| `DOWN_p` | Fisher's exact p-value |
| `partial_rho` | Partial Spearman correlation, controlling for log₁₀(3′ UTR length) |
| `partial_p` | p-value for partial Spearman |

The ZFP36L2 statistics in this file are independently reproducible from `02_ARE_annotations.csv` using the `n_7mer_L2` column. This cross-check is performed automatically when `ARE_analysis.py` runs.

## Methods summary

**Group definition.** Genes were classified by the SKY × Time interaction term in the ITT linear model: down (p < 0.05, log₂FC ≤ −0.3), up (p < 0.05, log₂FC ≥ +0.3), or background (all remaining expressed protein-coding genes). 3′ UTR sequences retrieved from Ensembl BioMart (GRCh38).

**ARE enrichment.** Genes with AREScore ≥ 4 in their 3′ UTR were classified as ARE-positive. AREScore is a composite weighting of AUUUA pentamer count, inter-pentamer proximity, and flanking AU-rich context (Spasic et al. 2012). Enrichment tested by Fisher's exact (downregulated vs background) with Wilson 95% CI for proportions and Woolf logit 95% CI for the odds ratio.

**Primary pre-specified test.** Partial Spearman correlation between UAUUUAU heptamer count and ITT log₂FC, controlling for log₁₀(3′ UTR length). UTR length correction is mandatory because downregulated genes have significantly longer UTRs than background (medians 2,418 vs 1,252 bp; Mann-Whitney p = 3.4 × 10⁻⁹) and longer UTRs naturally accumulate more motifs.

**RBP specificity panel.** Five RBP recognition motifs tested against identical 3′ UTR sequences: ZFP36L2 (UAUUUAU; Redmon et al. 2022), TTP 9-mer (UUAUUUAUUU), AUF1 (UUUUUUUU), HuR (UUUGUUUGU; PAR-CLIP consensus, Ray et al. 2013), KSRP (GGUGGG). KSRP is a GU-rich destabiliser and serves as the primary specificity control.

**Important note on attribution.** ZFP36L1, ZFP36L2, and TTP/ZFP36 share near-identical CCCH tandem zinc-finger domains and recognise the same UAUUUAU core motif (Hudson et al. 2004). Motif enrichment in this analysis therefore implicates the **ZFP36 family broadly**, not ZFP36L2 specifically; member-specific attribution requires expression or functional validation.

## Key statistics reproduced by this code

| Statistic | Value |
|---|---|
| Downregulated ARE-positive | 69 / 109 (63.3%) |
| Background ARE-positive | 2,223 / 9,588 (23.2%) |
| Fisher OR (down vs bg) | 5.72 (95% CI: 3.86–8.46) |
| Fisher p (down vs bg) | 9.08 × 10⁻¹⁹ |
| Upregulated negative control p | 0.043 (two-sided) |
| Partial Spearman ρ (7-mer vs log₂FC) | −0.156, p = 3.42 × 10⁻⁵⁴ |
| ZFP36L2 motif OR | 2.37, p = 9.65 × 10⁻⁶ |
| AUF1 motif OR | 2.41, p = 5.76 × 10⁻⁶ |
| HuR motif OR | 1.71, p = 0.093 (ns) |
| KSRP partial ρ | +0.254 (opposite to ZFP36L2) |

## Reproducibility note

The `n_7mer_L2` column and the four other RBP motif counts in `03_RBP_specificity.csv` were originally computed from 3′ UTR sequences retrieved from Ensembl BioMart at the time of the original analysis (Ensembl release 110). Re-running with current Ensembl annotation reproduces all reported odds ratios within 5%. Small differences arise from Ensembl updating transcript annotations between releases (~3 per year). For full upstream reproducibility, see `docs/methodology.md` for the BioMart query.

```
## Contact
Priyam Singh, PhD — [psingh27@uab.edu]
Issues and pull requests welcome via GitHub.
