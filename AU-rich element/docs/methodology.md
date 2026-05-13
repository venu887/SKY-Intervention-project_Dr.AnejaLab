# Detailed Methodology — ARE Analysis

## 1. Gene classification

Genes were classified into three groups based on the SKY × Time interaction term in the limma `voom` linear model fit to the RNA-seq count data.

- **Downregulated** (n = 109): nominal p < 0.05, log₂ fold-change ≤ −0.3, protein-coding, with annotated 3′ UTR
- **Upregulated negative control** (n = 35): nominal p < 0.05, log₂ fold-change ≥ +0.3, protein-coding, with annotated 3′ UTR
- **Background** (n = 9,588): all remaining expressed protein-coding genes with annotated 3′ UTR

Total: 9,732 genes after filtering for protein-coding biotype and successful 3′ UTR retrieval.

## 2. 3′ UTR retrieval

Sequences retrieved from Ensembl BioMart (GRCh38, Ensembl release 110 used in the original analysis) using the following query parameters:

- Dataset: `hsapiens_gene_ensembl`
- Attributes: `hgnc_symbol`, `ensembl_transcript_id`, `3utr`
- Filter: `transcript_biotype = "protein_coding"`

Representative transcript per gene selected as the transcript with the highest AUUUA pentamer count in its 3′ UTR (Fallmann et al. AREsite2, *Nucleic Acids Res* 2016). Minimum UTR length of 10 nt enforced (Gruber et al. 2011).

## 3. ARE scoring

AREScore computed using the algorithm of Spasic et al. (*PLoS Genet* 2012), accessed via the ARED-Plus database (Bakheet et al. *Nucleic Acids Res* 2018). AREScore is a composite weighting of:

1. AUUUA pentamer count
2. Inter-pentamer proximity
3. Flanking AU-rich context

Genes with AREScore ≥ 4 are classified as ARE-positive. This threshold corresponds to genes with strong destabilisation signal and is robust across thresholds 2–8 (OR range 4.21–5.89, all p < 10⁻¹³).

## 4. Statistical tests

### 4.1 Fisher's exact (primary group test)

Two-by-two contingency: ARE-positive count vs ARE-negative count, downregulated vs background. Two-sided p-value reported. Confidence interval on OR computed by Woolf logit method:

```
SE(log OR) = sqrt(1/a + 1/b + 1/c + 1/d)
95% CI    = exp(log OR ± 1.96 × SE)
```

### 4.2 Wilson 95% CI for proportions (error bars)

Standard Wilson score interval used for percentage error bars:

```
denom  = 1 + z²/n
centre = (p + z²/(2n)) / denom
margin = z × sqrt(p(1-p)/n + z²/(4n²)) / denom
CI     = [centre - margin, centre + margin]
```

### 4.3 Mann-Whitney U test (UTR length)

Two-sided test comparing 3′ UTR lengths between downregulated (n=109) and background (n=9,588). Required because downregulated genes have significantly longer UTRs (p = 3.4 × 10⁻⁹) which would inflate any motif-based enrichment if uncorrected.

### 4.4 Partial Spearman correlation (primary pre-specified test)

Following Mukherjee et al. (*Genome Biol* 2014), partial Spearman correlation between motif count and log₂FC, controlling for log₁₀(3′ UTR length).

Implementation: rank-transform x, y, and covariate; residualise x and y on the covariate using ordinary least squares on ranks; compute Spearman correlation of residuals.

```python
def partial_spearman(x, y, covariate):
    rx, ry, rc = rankdata(x), rankdata(y), rankdata(covariate)
    res_x = rx - linfit(rx, rc)
    res_y = ry - linfit(ry, rc)
    return spearmanr(res_x, res_y)
```

This is mathematically equivalent to the partial correlation method used in `ppcor::pcor.test()` in R for ranked data; numerical agreement is to within ~4% on test datasets due to subtle differences in rank-tie handling.

## 5. RBP motif specificity panel

Five RBP recognition motifs tested against identical 3′ UTR sequences to ensure the comparison is UTR-length-neutral:

| RBP | Motif (RNA) | Class | Citation |
|---|---|---|---|
| ZFP36L2 | UAUUUAU | Destabiliser | Redmon et al. 2022 |
| TTP 9-mer | UUAUUUAUUU | Destabiliser | TTP/ZFP36 paralog |
| AUF1 | UUUUUUUU | Destabiliser | Standard AUF1 motif |
| HuR | UUUGUUUGU | Stabiliser | Ray et al. 2013 (PAR-CLIP) |
| KSRP | GGUGGG | GU-rich destabiliser | Specificity control |

KSRP is the **primary specificity control** because it is a destabiliser like ZFP36L2 but recognises a GU-rich rather than AU-rich motif. If the ARE enrichment in downregulated genes were a generic AU/GU-rich signal rather than ZFP36-specific, KSRP would also be enriched. Instead KSRP shows the opposite direction (partial ρ = +0.254 vs ZFP36L2 partial ρ = −0.156), confirming AU-rich class specificity.

## 6. Important caveat on family attribution

ZFP36L1, ZFP36L2, and TTP/ZFP36 share near-identical CCCH tandem zinc-finger domains and bind the same UAUUUAU core motif (Hudson et al. *Nat Struct Mol Biol* 2004). The motif enrichment analysis cannot distinguish which specific family member is responsible for the observed signal in this dataset.

The manuscript therefore frames the finding as **ZFP36 family-mediated** decay rather than ZFP36L2-specific. HIF1A is reported as a validated target of ZFP36L1 (Ying et al. *Cancer Res* 2020), not ZFP36L2.

## 7. Differences from manuscript values

Values produced by this code may differ slightly from manuscript values due to:

1. **Ensembl annotation drift.** Ensembl releases new annotations ~3 times per year. Re-running the BioMart query in a different release will return different transcript isoforms for some genes, producing slightly different motif counts. All reported OR values reproduce within 5%.

2. **Partial Spearman implementation.** The manuscript values were computed using `ppcor::pcor.test()` in R; this code uses a Python implementation of the rank-residual approach. Both are valid; numerical agreement is within ~4%.

These differences do not affect any directional conclusions or significance calls.

## 8. References

1. Bailey TL, et al. *MEME SUITE: tools for motif discovery and searching.* Nucleic Acids Res. 2009;37:W202–W208.
2. Bakheet T, et al. *ARED-Plus.* Nucleic Acids Res. 2018;46(D1):D218–D220.
3. Hudson BP, et al. *Recognition of mRNA AU-rich element by ZFP36L2.* Nat Struct Mol Biol. 2004;11(3):257–264.
4. Mukherjee N, et al. *Global target mRNA specification by ZFP36.* Genome Biol. 2014;15(1):R12.
5. Ray D, et al. *A compendium of RNA-binding motifs.* Nature. 2013;499:172–177.
6. Redmon SN, et al. *Sequence and tissue targeting specificity of ZFP36L2.* Nucleic Acids Res. 2022;50(7):4068–4082.
7. Spasic M, et al. *Genome-wide assessment of AU-rich elements by AREScore.* PLoS Genet. 2012;8(1):e1002433.
8. Ying Y, et al. *ZFP36L1 suppresses hypoxia and cell-cycle signaling.* Cancer Res. 2020;80(2):219–232.
