# Methods — scOnco-scorer

A purity-resolved benchmark of single-cell CNV callers on the FDA SEQC2 HCC1395 reference, returning a calibrated per-cell malignancy probability scored against a DNA copy-number truth.

---

## 1 · Design

**Question.** How well do single-cell CNV callers separate cancer from normal cells as tumour purity drops, and is their per-cell probability calibrated?

**Reference.** FDA SEQC2 **HCC1395** (triple-negative breast cancer) versus **HCC1395BL** (matched normal B-lymphoblastoid line from the same donor). GRCh38 throughout. Because both lines come from one donor, per-cell ground-truth labels are known by construction.

**Two arms, kept separate.**

| Arm | Platform | Callers |
|---|---|---|
| Expression | 10x | inferCNV, CopyKAT, SCEVAN |
| Allele | Fluidigm C1, SMART-Seq v4 | Numbat |

Both arms are scored against the same DNA CNA truth. The arms use different cells and are never fused per cell.

**Novelty.** Purity-resolved evaluation **+** a calibrated per-cell malignancy probability with an indeterminate band **+** scoring against the full DNA CNA truth, plus five analyses this titration uniquely enables (Section 5).

---

## 2 · In-silico purity series and mixture construction

**Purity ladder.** Tumour fraction **1, 5, 10, 20, 40, 80 %**.

**Site-matched, within-capture reference (the central confound fix).** Mixtures are built **within one site**: the LLU-B capture is partitioned by cell into a **reference half** and a disjoint **scored half**. Each mixture is `LLU-A + scored-half-B`, and `ref-half-B` is the caller reference. A, scored-B, and ref-B all come from LLU, so the only systematic tumour-versus-normal difference is biology, not batch. The NCI site is run as an **independent replicate**, never as a cross-site reference; a cross-site reference injects a batch axis the callers read as copy number.

**Cohort.** 3000 cells per context, sampled **without replacement** within a context. Absolute tumour-cell support: 30 / 150 / 300 / 600 / 1200 / 2400 across the ladder. A context is **refused rather than downsampled** if the required tumour or normal count exceeds the pure pool remaining after the reference holdout.

**Replicates and nested titration.** At least **10 independently seeded** mixture replicates per site and fraction. Within a replicate, one random ordering of A and B cells is drawn and every fraction is built from **prefixes** of that ordering, so the same cell is scored across fractions and cohort-composition effects are observable.

**Physical mixtures.** Real **5 %** and **10 %** mixtures provide an orthogonal check with no per-cell source truth; their labels are recovered by the demux in Section 4, and physical and in-silico results are never pooled into one metric.

---

## 3 · Truth, masking, and identifier convention

**CNA truth.** GRCh38 copy-number **segment** truth from the SEQC2 bulk WGS (ASCAT/consensus), as gain / loss / neutral state per region; the Masood 2024 consensus (**346 gain, 33 loss, 320 LOH**) is scored as state per region, not integer copy number. Per-replicate **ascatNgs** segments give a single-caller cross-check. The somatic SNV/indel VCFs are **not** CNA truth; they are point mutations and are used only for the somatic-SNV classifier (Section 5), not for copy-number scoring.

**Masking (genome-internal negative control).** chr **6p, 16q, X** are masked from all scoring. HCC1395BL is EBV-transformed and itself carries an unbalanced t(6;16) (net loss of 6p and 16q) and loss of X, so on those arms the normal reference is itself lost and tumour-versus-normal separation is uninterpretable. SEQC2 excludes these arms from the DNA truth for this reason. Their differential-signal rate is reported separately as a per-caller false-positive readout.

**Identifiers and build.** One convention: **version-stripped Ensembl** gene IDs everywhere (strip `\.\d+$` at load, collapse duplicate rows by summing). Gene positions come from the **refdata-cellranger-GRCh38-1.2.0** (Ensembl v84) GTF that produced the count matrices. CopyKAT `id.type="E"`, `genome="hg20"` (its name for GRCh38); inferCNV four-column gene-order file on the same IDs; SCEVAN `norm_cell` on the same IDs. Every tool emits input / mapped / unmapped / duplicate-collapsed / final gene counts and **fails loudly** if mapped genes fall below threshold; unstripped versions cause silent mass gene loss and a shifted baseline. Build is asserted GRCh38 before any intersection; `1` vs `chr1` is harmonized only after that assertion, with REF alleles checked against the FASTA.

---

## 4 · Caller pipelines

**Expression arm (10x).**

```
mixture matrix -> inferCNV / CopyKAT / SCEVAN -> per-cell CNV signal -> calibrated P(malignant)
```

Callers run with `ref-half-B` as the internal reference. inferCNV is also run reference-free as a separate condition, since a fixed reference suppresses its miscentring failure mode. CopyKAT returns a binary label, so a continuous score (mean absolute CNA per cell from `CNAmat`) is derived for AUC. SCEVAN runs `FIXED_NORMAL_CELLS=FALSE`; retaining its continuous score needs a pinned wrapper around its preprocessing/classify steps.

**Allele arm (Fluidigm C1).**

```
37 C1 libraries (20 A + 17 B) -> STAR align (-> per-cell BAMs) -> pileup + phase (--smartseq) -> Numbat -> per-unit posterior
```

The C1 featureCounts matrix is Numbat's **expression** input (`count_mat`); the **allele** input (`df_allele`) comes from the C1 BAMs via `pileup_and_phase.R --smartseq` against the bundled 1000G hg38 panel. `lambdas_ref` is built from the **pure-B reference subset only**, no tumour cells. Settings: `--smartseq` (plate-based, no CB/UMI), `gamma=5` (non-UMI SMART-Seq; the default 20 is for 10x), `genome="hg38"`. LOH handling: prefer known clonal-LOH/consensus segments from the bulk WGS via `segs_loh` / `segs_consensus_fix`; use `call_clonal_loh=TRUE` only as fallback when DNA segments are not wired in; never supply both. A 20-A + 20-B pilot runs before all cells, since the SRA download is the main wall-clock uncertainty; the C1 arm is pure-line only, so its purity evaluation comes from in-silico pooling of C1 A and B, not from physical mixtures.

**Both arms** return a calibrated per-cell / per-unit malignancy call, scored against the same CNA truth across the purity sweep.

---

## 5 · Calibration, controls, and analyses

**Calibrated probability with a learned indeterminate band.** Each caller's score is calibrated (Platt / beta / isotonic, chosen by inner CV on log loss / Brier) under **nested, source-cell-grouped** cross-validation: complete source cells are held out across every context they appear in, calibration is fit on training cells only, and evaluation is once on outer held-out cells. Two thresholds are **prespecified** before final evaluation: a lower threshold for a target NPV and an upper threshold for a target PPV; cells between them are **indeterminate**. Calibration is proven against the pure-A/B labels; DNA truth separately tests genomic fidelity.

**Uncertainty.** Because cells recur across fractions and replicates, resampling is by **source-cell clusters stratified by true label** (row-level bootstrap is invalid), at least 2000 replicates, aggregated over mixture replicates. Every per-fraction metric carries a 95 % cluster-bootstrap CI; a detectability crossover counts only when adjacent-fraction CIs do not overlap.

**Controls (the gate).** Three must pass or downstream numbers are void: a held-out pure-B half scored against the B reference returns near **0 % aneuploid**; the **label-shuffle** drops AUC to ~0.5; the chr6p/16q/X **false-positive spectrum** stays low and flat.

**Analyses this titration enables.** (A) Per-arm, per-caller **detectability floor**: lowest fraction at which an arm's tumour-vs-normal AUC ≥ 0.8, predicted to track the arm's WGS log-ratio amplitude. (B) **Reliability-weighted fusion** of expression and allele posteriors (running Numbat on the 10x cells so four posteriors share a cell), judged on held-out Brier/ECE against the best single caller. (C) Independent **somatic-SNV classifier** at the ~40k SEQC2 somatic sites as a high-precision anchor, with pure-B alt-rate at the error floor as its falsifier. (D) Binomial **clonal-versus-subclonal** test per event against a clonal null with sensitivity estimated from a known-clonal event, BH-corrected. (E) **Calibration drift**: ECE per fraction, locating the calibration floor separately from the accuracy floor.

**Determinism.** One master `--seed` enters mixture building and is written to a seed ledger; each fraction gets a derived seed (`seed + round(f*1000)`); the reference/scored-B split, every bootstrap, every calibration split, and every label-shuffle take the recorded seed. Rerunning from the ledger reproduces mixtures and scores exactly.

---

*Repository: github.com/collaborativebioinformatics/scOnco-scorer*
