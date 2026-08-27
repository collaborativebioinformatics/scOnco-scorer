# Methods — scOnco-scorer

A purity-resolved benchmark of single-cell malignancy calling on the FDA SEQC2 HCC1395 reference.

---

## 1 · Design

**Question.** How well do single-cell CNV callers separate cancer from normal cells as tumour purity drops?

**Reference.** FDA SEQC2 **HCC1395** (triple-negative breast cancer) versus **HCC1395BL** (matched normal B-lymphoblastoid line from the same donor). GRCh38 throughout. Because the two lines come from one donor, per-cell ground-truth labels are known.

**Tumour-purity series.** In-silico mixtures spanning **1, 5, 10, 20, 40, 80 %** tumour fraction. Each mixture is built within a single site: a disjoint subset of normal-B cells is held out as the caller reference, and the remaining cells are mixed (tumour-A + scored-B) at the target fraction. Real **5 %** and **10 %** physical mixtures serve as an orthogonal check.

**Two arms, kept separate.**

| Arm | Platform | Callers |
|---|---|---|
| Expression | 10x | inferCNV, CopyKAT, SCEVAN |
| Allele | Fluidigm C1 | Numbat |

Both arms are scored against the same DNA CNA truth.

**Novelty.** Purity-resolved evaluation **+** a calibrated per-cell malignancy probability **+** scoring against the full DNA CNA truth.

---

## 2 · Truth and scoring

**CNA truth.** The Masood 2024 SEQC2 consensus call set: **346 gain, 33 loss, 320 LOH**, derived from six callers plus orthogonal validation. Scored as **state per region**, not integer copy number.

**Cross-check.** Per-replicate ascatNgs segments provide a single-caller comparison.

**Masking.** chr **6p / 16q / X** are masked, since they are aberrant in the normal line itself.

**Per-cell output.** Each caller returns a **calibrated malignancy probability** with an **indeterminate band** defined by two thresholds fixed before the test:

```
P(malignant)   0 ─────────────┬──────────────┬───────────────→ 1
                   normal      │ indeterminate │   malignant
```

Cells that fall in the band are called indeterminate rather than forced to a label; this is what keeps the score reliable under low purity.

**Metrics.** AUROC and calibration on the in-silico purity sweep. A label-shuffle control yields AUC ≈ 0.5.

---

## 3 · Caller pipelines

**Expression arm (10x).**

```
mixture matrix → inferCNV / CopyKAT / SCEVAN → per-cell CNV signal → calibrated P(malignant)
```

Held-out normal-B cells act as the internal reference for each caller.

**Allele arm (Fluidigm C1).**

```
37 C1 libraries → align (→ BAMs) → pileup + phase → Numbat → per-unit posterior
```

37 pooled libraries at run level: 20 tumour (A) + 17 normal (B). Allele-based, independent of the expression arm.

**Both arms** return a calibrated per-cell / per-unit malignancy call, scored against the same CNA truth across the purity sweep.

---

*Repository: github.com/collaborativebioinformatics/scOnco-scorer*
