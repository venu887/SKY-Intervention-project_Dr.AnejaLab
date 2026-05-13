"""
Sanity checks for ARE_analysis.py
===================================
Run this BEFORE submitting code with manuscript.

Usage:
    python3 test_ARE_analysis.py

All tests must print PASS. Any FAIL means a bug in the analysis.

Tests are grouped into 4 levels:
  Level 1 — Data integrity (is the input file what we expect?)
  Level 2 — Statistical correctness (do the numbers match expected values?)
  Level 3 — Directional logic (do results make biological sense?)
  Level 4 — Robustness (do edge cases break anything?)
"""

import sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr, norm

# ── import the functions we are testing ───────────────────────────────────
sys.path.insert(0, '.')
from ARE_analysis import (
    load_data, run_are_statistics, partial_spearman,
    wilson_ci, or_confidence_interval
)

PASS = "\033[92mPASS\033[0m"
FAIL = "\033[91mFAIL\033[0m"
n_pass = 0
n_fail = 0

def check(name, condition, detail=""):
    global n_pass, n_fail
    if condition:
        print(f"  {PASS}  {name}")
        n_pass += 1
    else:
        print(f"  {FAIL}  {name}  {detail}")
        n_fail += 1

def approx(a, b, tol=0.01):
    """Check two values are within relative tolerance."""
    if b == 0:
        return abs(a) < tol
    return abs(a - b) / abs(b) < tol

def section(title):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


# ── Load data ─────────────────────────────────────────────────────────────
are, rbp = load_data("data/02_ARE_annotations.csv", "data/03_RBP_specificity.csv")
down = are[are['group'] == 'down']
bg   = are[are['group'] == 'background']
up   = are[are['group'] == 'up']
stats = run_are_statistics(are)


# ══════════════════════════════════════════════════════════════════════════
# LEVEL 1 — DATA INTEGRITY
# ══════════════════════════════════════════════════════════════════════════
section("LEVEL 1: Data integrity")

# 1.1 Group sizes
check("Group size: downregulated n=109",    len(down) == 109)
check("Group size: background n=9,588",     len(bg)   == 9588)
check("Group size: upregulated n=35",       len(up)   == 35)
check("Total genes: 9,732",                 len(are)  == 9732)

# 1.2 No overlapping groups
down_genes = set(down['Gene'])
bg_genes   = set(bg['Gene'])
up_genes   = set(up['Gene'])
check("No genes in both down and background",
      len(down_genes & bg_genes) == 0)
check("No genes in both up and background",
      len(up_genes & bg_genes) == 0)
check("No genes in both down and up",
      len(down_genes & up_genes) == 0)

# 1.3 All genes have a 3' UTR length
check("All genes have UTR length > 0",
      (are['utr_length'] > 0).all())
check("No NaN in AREScore",
      are['AREScore'].isna().sum() == 0)
check("No NaN in n_7mer_L2",
      are['n_7mer_L2'].isna().sum() == 0)
check("No NaN in Intervention_logFC",
      are['Intervention_logFC'].isna().sum() == 0)

# 1.4 Correct directional logFC assignments
check("All downregulated genes have negative logFC",
      (down['Intervention_logFC'] < 0).all())
check("All upregulated genes have positive logFC",
      (up['Intervention_logFC'] > 0).all())

# 1.5 AREScore is non-negative
check("AREScore is non-negative everywhere",
      (are['AREScore'] >= 0).all())
check("7-mer count is non-negative everywhere",
      (are['n_7mer_L2'] >= 0).all())

# 1.6 Threshold choice matters — check ARE+ counts
d_pos = (down['AREScore'] >= 4).sum()
b_pos = (bg['AREScore']   >= 4).sum()
check("Downregulated ARE+ count = 69",   d_pos == 69)
check("Background ARE+ count = 2,223",   b_pos == 2223)

# 1.7 ARE+ rate in downregulated > background (basic sanity)
check("Down ARE+ rate > background ARE+ rate",
      d_pos/len(down) > b_pos/len(bg))

# 1.8 ARE+ rate in upregulated < downregulated (negative control works)
u_pos = (up['AREScore'] >= 4).sum()
check("Up ARE+ rate < down ARE+ rate",
      u_pos/len(up) < d_pos/len(down))


# ══════════════════════════════════════════════════════════════════════════
# LEVEL 2 — STATISTICAL CORRECTNESS
# ══════════════════════════════════════════════════════════════════════════
section("LEVEL 2: Statistical correctness")

# 2.1 Fisher exact values match manuscript
check("Fisher OR (down/bg) = 5.72 ± 1%",
      approx(stats['or_down_bg'], 5.72),
      f"got {stats['or_down_bg']:.4f}")
check("Fisher CI lower (down/bg) = 3.86 ± 1%",
      approx(stats['ci_lo_down_bg'], 3.86),
      f"got {stats['ci_lo_down_bg']:.4f}")
check("Fisher CI upper (down/bg) = 8.46 ± 1%",
      approx(stats['ci_hi_down_bg'], 8.46),
      f"got {stats['ci_hi_down_bg']:.4f}")
check("Fisher p (down/bg) < 1e-15",
      stats['p_down_bg'] < 1e-15,
      f"got {stats['p_down_bg']:.2e}")
check("Fisher p (down/bg) ≈ 9.08e-19 ± 5%",
      approx(stats['p_down_bg'], 9.08e-19, tol=0.05),
      f"got {stats['p_down_bg']:.2e}")

check("Fisher OR (up/bg) = 0.31 ± 1%",
      approx(stats['or_up_bg'], 0.31),
      f"got {stats['or_up_bg']:.4f}")
check("Fisher p (up/bg) < 0.05",
      stats['p_up_bg'] < 0.05,
      f"got {stats['p_up_bg']:.4f}")

# 2.2 Wilson CI sanity
check("Down CI lower < down pct < down CI upper",
      stats['d_ci_lo'] < stats['d_pct'] < stats['d_ci_hi'])
check("BG CI lower < bg pct < bg CI upper",
      stats['b_ci_lo'] < stats['b_pct'] < stats['b_ci_hi'])
check("Down CI width reasonable (5-20 pp)",
      5 < (stats['d_ci_hi'] - stats['d_ci_lo']) < 20)

# 2.3 UTR length statistics
check("Down UTR median = 2,418 bp ± 1%",
      approx(stats['d_utr_median_bp'], 2418),
      f"got {stats['d_utr_median_bp']:.0f}")
check("BG UTR median = 1,252 bp ± 1%",
      approx(stats['b_utr_median_bp'], 1252),
      f"got {stats['b_utr_median_bp']:.0f}")
check("UTR Mann-Whitney p < 1e-5",
      stats['mannwhitney_p_utr'] < 1e-5,
      f"got {stats['mannwhitney_p_utr']:.2e}")
check("Down UTR longer than background",
      stats['d_utr_median_bp'] > stats['b_utr_median_bp'])

# 2.4 Partial Spearman direction
check("Partial rho (7-mer vs logFC) is negative",
      stats['partial_rho_7mer'] < 0,
      f"got {stats['partial_rho_7mer']:.4f}")
check("Partial rho (AREScore vs logFC) is negative",
      stats['partial_rho_are'] < 0,
      f"got {stats['partial_rho_are']:.4f}")
check("Partial p (7-mer) < 1e-40",
      stats['partial_p_7mer'] < 1e-40,
      f"got {stats['partial_p_7mer']:.2e}")

# 2.5 CDF Fisher
check("CDF Fisher OR (>=2 vs <2) ≈ 3.37 ± 2%",
      approx(stats['cdf_or_ge2_vs_lt2'], 3.37, tol=0.02),
      f"got {stats['cdf_or_ge2_vs_lt2']:.4f}")
check("CDF Fisher p (>=2 vs <2) < 1e-5",
      stats['cdf_p_ge2_vs_lt2'] < 1e-5,
      f"got {stats['cdf_p_ge2_vs_lt2']:.2e}")

# 2.6 RBP panel values from CSV
zbp = rbp[rbp['RBP'] == 'ZFP36L2'].iloc[0]
ksrp = rbp[rbp['RBP'] == 'KSRP'].iloc[0]
hur  = rbp[rbp['RBP'] == 'HuR'].iloc[0]
auf1 = rbp[rbp['RBP'] == 'AUF1'].iloc[0]

check("ZFP36L2 OR ≈ 2.37",        approx(zbp['DOWN_OR'], 2.37),
      f"got {zbp['DOWN_OR']:.4f}")
check("ZFP36L2 p < 0.05",         zbp['DOWN_p'] < 0.05)
check("ZFP36L2 partial rho < 0",  zbp['partial_rho'] < 0)
check("AUF1 OR ≈ 2.41",           approx(auf1['DOWN_OR'], 2.41),
      f"got {auf1['DOWN_OR']:.4f}")
check("AUF1 p < 0.05",            auf1['DOWN_p'] < 0.05)
check("HuR p > 0.05 (ns)",        hur['DOWN_p'] > 0.05,
      f"got {hur['DOWN_p']:.4f}")
check("KSRP partial rho > 0",     ksrp['partial_rho'] > 0,
      f"got {ksrp['partial_rho']:.4f}  — must be opposite to ZFP36L2")
check("KSRP OR < 1",              ksrp['DOWN_OR'] < 1,
      f"got {ksrp['DOWN_OR']:.4f}  — must be depleted not enriched")


# ══════════════════════════════════════════════════════════════════════════
# LEVEL 3 — DIRECTIONAL LOGIC (BIOLOGICAL SANITY)
# ══════════════════════════════════════════════════════════════════════════
section("LEVEL 3: Directional logic")

# 3.1 ARE enrichment goes in the right direction
check("Down ARE+ rate > bg (primary finding)",
      stats['d_pct'] > stats['b_pct'])
check("Up ARE+ rate < bg (negative control depleted)",
      stats['u_pct'] < stats['b_pct'])
check("OR down/bg > 1",  stats['or_down_bg'] > 1)
check("OR up/bg < 1",    stats['or_up_bg']   < 1)

# 3.2 More 7-mers = more downregulated (dose-response)
ge2  = are[are['n_7mer_L2'] >= 2]['Intervention_logFC'].mean()
one  = are[are['n_7mer_L2'] == 1]['Intervention_logFC'].mean()
zero = are[are['n_7mer_L2'] == 0]['Intervention_logFC'].mean()
check("Mean logFC: >=2 7-mers < 1 7-mer < 0 7-mers (dose-response)",
      ge2 < one < zero,
      f"ge2={ge2:.3f}, one={one:.3f}, zero={zero:.3f}")

# 3.3 Destabilisers enriched, stabiliser not
dest_ors = rbp[rbp['Role']=='destabiliser']['DOWN_OR'].values
stab_ors = rbp[rbp['Role']=='stabiliser']['DOWN_OR'].values
check("Mean destabiliser OR > 1",
      dest_ors.mean() > 1,
      f"mean={dest_ors.mean():.3f}")
check("Stabiliser (HuR) OR < destabiliser mean",
      stab_ors.mean() < dest_ors.mean())

# 3.4 KSRP is the specificity control — must go opposite direction
check("KSRP partial rho opposite sign to ZFP36L2 partial rho",
      np.sign(ksrp['partial_rho']) != np.sign(zbp['partial_rho']))

# 3.5 Partial rho < raw Spearman (UTR correction reduces correlation)
raw_rho, _ = spearmanr(are['n_7mer_L2'], are['Intervention_logFC'])
check("Raw Spearman |rho| > partial |rho| (UTR correction reduces estimate)",
      abs(raw_rho) > abs(stats['partial_rho_7mer']),
      f"raw={raw_rho:.4f}, partial={stats['partial_rho_7mer']:.4f}")

# 3.6 UTR length confound is real
check("Downregulated genes have longer UTRs than background",
      stats['d_utr_median_bp'] > stats['b_utr_median_bp'])
check("UTR confound is significant (p<0.001)",
      stats['mannwhitney_p_utr'] < 0.001)

# 3.7 HIF1A and CXCL3 are in the downregulated set
hif = are[are['Gene'] == 'HIF1A']
cxcl = are[are['Gene'] == 'CXCL3']
check("HIF1A is in dataset",   len(hif) > 0)
check("CXCL3 is in dataset",   len(cxcl) > 0)
check("HIF1A group = down",    hif.iloc[0]['group'] == 'down' if len(hif) else False)
check("CXCL3 group = down",    cxcl.iloc[0]['group'] == 'down' if len(cxcl) else False)
check("HIF1A AREScore = 8",    hif.iloc[0]['AREScore'] == 8 if len(hif) else False)
check("CXCL3 AREScore = 5",    cxcl.iloc[0]['AREScore'] == 5 if len(cxcl) else False)
check("HIF1A n_7mer = 3",      hif.iloc[0]['n_7mer_L2'] == 3 if len(hif) else False)
check("CXCL3 n_7mer = 3",      cxcl.iloc[0]['n_7mer_L2'] == 3 if len(cxcl) else False)


# ══════════════════════════════════════════════════════════════════════════
# LEVEL 4 — ROBUSTNESS (EDGE CASES)
# ══════════════════════════════════════════════════════════════════════════
section("LEVEL 4: Robustness")

# 4.1 Wilson CI with extreme proportions
lo, hi = wilson_ci(0, 100)
check("Wilson CI with k=0 doesn't crash and lower=0",
      lo == 0.0 and hi > 0)
lo, hi = wilson_ci(100, 100)
check("Wilson CI with k=n doesn't crash and upper=100",
      hi > 95 and lo < 100)

# 4.2 Partial Spearman with known input
# If x and covariate are uncorrelated, partial rho ≈ raw rho
np.random.seed(42)
x   = np.random.randn(500)
y   = -0.3 * x + np.random.randn(500) * 0.9
cov = np.random.randn(500)   # independent covariate
rho_partial, _ = partial_spearman(x, y, cov)
rho_raw, _     = spearmanr(x, y)
check("Partial Spearman ≈ raw Spearman when covariate is independent",
      abs(rho_partial - rho_raw) < 0.05,
      f"partial={rho_partial:.4f}, raw={rho_raw:.4f}")

# 4.3 Partial Spearman removes confound correctly
# Build dataset where correlation is entirely due to a covariate
np.random.seed(0)
cov = np.random.randn(1000)
x   = cov + np.random.randn(1000) * 0.1      # x almost = covariate
y   = -cov + np.random.randn(1000) * 0.1     # y almost = -covariate
# Raw rho(x,y) should be strongly negative (both driven by cov)
rho_raw2, _ = spearmanr(x, y)
# Partial rho(x,y | cov) should be near zero (no direct x->y relationship)
rho_partial2, _ = partial_spearman(x, y, cov)
check("Partial Spearman removes covariate-driven correlation",
      abs(rho_raw2) > 0.5 and abs(rho_partial2) < 0.2,
      f"raw={rho_raw2:.3f}, partial={rho_partial2:.3f}")

# 4.4 Fisher exact symmetry
or1, p1 = fisher_exact([[69, 40], [2223, 7365]])
or2, p2 = fisher_exact([[2223, 7365], [69, 40]])
check("Fisher exact: swapping rows inverts OR",
      approx(or1, 1/or2, tol=0.001))
check("Fisher exact: swapping rows gives same p-value",
      approx(p1, p2, tol=0.001))

# 4.5 OR CI contains OR
lo, hi = or_confidence_interval(69, 40, 2223, 7365)
or_val, _ = fisher_exact([[69, 40], [2223, 7365]])
check("OR confidence interval contains the OR",
      lo < or_val < hi,
      f"OR={or_val:.4f}, CI=({lo:.4f},{hi:.4f})")

# 4.6 AREScore threshold sensitivity
# Check that enrichment holds across thresholds 2-8
print("\n  AREScore threshold sensitivity (should all show OR > 1, p < 0.05):")
all_sig = True
for t in range(2, 9):
    dp = (down['AREScore'] >= t).sum()
    dn = len(down) - dp
    bp = (bg['AREScore']   >= t).sum()
    bn = len(bg)   - bp
    if dp == 0 or bp == 0:
        continue
    or_t, p_t = fisher_exact([[dp, dn], [bp, bn]])
    sig = p_t < 0.05 and or_t > 1
    all_sig = all_sig and sig
    print(f"    threshold={t}: OR={or_t:.2f}, p={p_t:.1e}  {'OK' if sig else 'FAIL'}")
check("ARE enrichment significant at all thresholds 2-8", all_sig)


# ══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════
total = n_pass + n_fail
print(f"\n{'='*60}")
print(f"  RESULTS: {n_pass}/{total} tests passed, {n_fail} failed")
if n_fail == 0:
    print("  \033[92mALL TESTS PASSED — code is ready for submission\033[0m")
else:
    print("  \033[91mFAILURES DETECTED — do not submit until resolved\033[0m")
print(f"{'='*60}\n")

sys.exit(0 if n_fail == 0 else 1)
