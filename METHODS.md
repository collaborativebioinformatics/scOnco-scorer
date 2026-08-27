<div align="center">

# scOnco-scorer · Methods

**Purity-resolved benchmarking of single-cell copy-number callers on the FDA SEQC2 HCC1395/HCC1395BL reference system**

![Genome](https://img.shields.io/badge/genome-GRCh38-4C78A8)
![Truth](https://img.shields.io/badge/truth-SEQC2%20DNA%20CNV-6F42C1)
![10x](https://img.shields.io/badge/expression-10x%20Genomics-1F8A70)
![C1](https://img.shields.io/badge/allele-Fluidigm%20C1%20SMART--Seq%20v4-C76D2D)
![Reproducibility](https://img.shields.io/badge/reproducibility-seeded%20%2B%20ledgered-2E8B57)

</div>

> [!IMPORTANT]
> **Calibration and genomic fidelity are evaluated against different truth layers.**
> Per-cell malignancy calibration is evaluated against known source labels (HCC1395 versus HCC1395BL). Copy-number event fidelity is evaluated independently against SEQC2 DNA CNV/LOH truth.

> [!CAUTION]
> **The 10x and Fluidigm C1 arms contain different cells.** They are compared at the benchmark level and are never fused as if they were measurements of the same cell.

---

## 1 · Benchmark question and design

### Primary question

How well do single-cell CNV callers distinguish tumour from matched-normal cells as tumour purity falls, and how well calibrated are their per-cell malignancy probabilities?

### Reference system

The benchmark uses the FDA SEQC2 paired reference cell lines:

- **HCC1395 (sample A):** triple-negative breast cancer cell line.
- **HCC1395BL (sample B):** matched EBV-transformed B-lymphoblastoid line from the same donor.
- **Genome build:** GRCh38 throughout.

The two lines share the donor germline, while **source labels are known by construction for the pure captures and in-silico mixtures**. The shared donor origin is therefore a matched-normal property; it is not the reason the single-cell source labels are known.

### Two benchmark arms

| Arm | Platform | Primary callers | Signal |
|---|---|---|---|
| **Expression** | 10x Genomics scRNA-seq | inferCNV · CopyKAT · SCEVAN | expression-derived CNV |
| **Allele** | Fluidigm C1 · SMART-Seq v4 | Numbat | expression + phased allele imbalance |

The arms are analysed independently and compared against the same DNA-derived genomic truth framework. They do **not** share individual cells.

### Benchmark contributions

The benchmark is designed to quantify:

1. **Purity-dependent discrimination** of malignant versus normal cells.
2. **Per-cell probability calibration**, including an explicit indeterminate interval.
3. **Genomic fidelity** against orthogonally supported SEQC2 CNV/LOH truth.
4. **Context sensitivity**, detectability limits, and calibration drift across tumour fractions.
5. **Independent allele/SNV evidence** where read-level data permit it.

---

## 2 · In-silico purity series

### Purity ladder

The primary 10x purity series uses tumour fractions:

**1%, 5%, 10%, 20%, 40%, 80%**

Each context contains **3,000 scored cells**, sampled without replacement within that context:

| Tumour fraction | HCC1395 cells | HCC1395BL cells |
|---:|---:|---:|
| 1% | 30 | 2,970 |
| 5% | 150 | 2,850 |
| 10% | 300 | 2,700 |
| 20% | 600 | 2,400 |
| 40% | 1,200 | 1,800 |
| 80% | 2,400 | 600 |

A context is **refused rather than silently resized or sampled with replacement** if the required A or B pool is not available after the reference holdout.

### Site-matched reference design

The primary benchmark never uses a cross-site normal reference.

For each site independently:

1. Partition the pure-B capture into a **reference-B subset** and a disjoint **scored-B subset**.
2. Construct each mixture from `A + scored-B`.
3. Supply only the same-site `reference-B` subset as the caller reference.
4. Keep the NCI and LLU series as separate replication axes.

This prevents the caller from confusing a sequencing-site/capture effect with copy-number signal.

### Replicates and nested titration

At least **10 independently seeded mixture replicates** are generated per site. Within each replicate, A and B cells are randomly ordered once, and fractions are constructed from deterministic prefixes of those orderings. This allows the same source cell to recur across fractions and makes cohort-composition sensitivity measurable.

> [!NOTE]
> Repeated appearances of the same source cell are treated as repeated observations of one biological unit during uncertainty estimation; they are not bootstrapped as independent rows.

### Physical mixtures

Real **5%** and **10%** tumour mixtures are analysed separately from the in-silico series. Their nominal mixture fraction is known, but per-cell source identity is not directly observed.

Because HCC1395 and HCC1395BL come from the same donor, **standard de novo donor demultiplexing is not used as source truth**. Physical-mixture source assignment is instead based on fixed tumour-informative allele evidence, including the SEQC2 somatic-SNV set and copy-number/LOH-aware allelic imbalance. Cells with insufficient evidence remain unassigned rather than being forced into A or B.

Physical-mixture results and in-silico source-labelled results are reported separately.

---

## 3 · DNA truth, masking, and coordinate conventions

### Primary CNA/LOH truth

The primary genomic truth is the **GRCh38 SEQC2 HCC1395 CNV benchmark call set** from Masood *et al.* (2024), constructed from consensus NGS evidence and orthogonal validation. The published benchmark contains:

- **346 high-confidence gain intervals** spanning 1,525.6 Mb;
- **33 high-confidence loss intervals** spanning 87.9 Mb;
- **320 high-confidence LOH intervals** spanning 1,490.4 Mb.

The benchmark evaluates **regional state concordance** rather than pretending these consensus intervals provide single-cell integer copy number.

Per-replicate **ascatNgs** segments from SEQC2 are retained as a secondary single-caller sensitivity analysis. HCC1395 also contains documented clonal copy-neutral LOH on **17q**, which is relevant to the allele arm and LOH handling.

> [!WARNING]
> Somatic SNV and indel VCFs are **not CNA truth**. They are point-mutation resources and are used only for fixed-site allelic classification / orthogonal somatic evidence, never as copy-number segments.

### HCC1395BL abnormal regions

HCC1395BL is not copy-number-normal across the whole genome. SEQC2 cytogenetic and array analyses show:

- loss of **chr6p**;
- loss of **chr16q**;
- loss of **chrX**;
- an unbalanced chr6/chr16 rearrangement underlying the 6p/16q losses.

These regions are therefore **excluded from tumour-versus-normal CNA scoring** because the matched-normal baseline is itself abnormal. SEQC2 likewise removed these normal-LOH regions from high-confidence somatic benchmarking regions/call sets.

These arms may be reported as a **masked-region diagnostic**, but they are **not treated as a false-positive ground truth** and are not required to be “silent”: tumour-versus-normal interpretation there is intrinsically ambiguous.

### Gene identifiers

The internal gene-ID convention is **version-stripped Ensembl gene IDs**:

```text
ENSG00000141510.17  ->  ENSG00000141510
```

At matrix load:

1. restore the gene dimension explicitly from `dimnames(x)[[1]]` when needed;
2. strip Ensembl version suffixes with `\.\d+$`;
3. collapse duplicate Ensembl IDs by **summing count rows**;
4. record input, mapped, unmapped, duplicate-collapsed, and final gene counts.

Tool-specific conventions are kept explicit:

| Tool | Identifier / genome convention |
|---|---|
| CopyKAT | `id.type="E"`, `genome="hg20"` — CopyKAT's GRCh38 label |
| inferCNV | four-column gene-order file keyed to the same version-stripped Ensembl IDs |
| SCEVAN | same Ensembl ID space; same-site B cells supplied through `norm_cell` |
| Numbat | custom/reference GTF must match the expression-matrix gene identifiers used by the run |

### Genome-build and contig validation

Every genomic resource is asserted to be **GRCh38 before intersection**. The pipeline does not infer genome build from whether chromosomes are named `1` or `chr1`.

`1`/`chr1` normalization is performed only after:

- FASTA sequence-dictionary inspection;
- coordinate-build validation;
- representative VCF REF-allele checks against the FASTA when applicable.

No resource is silently lifted over.

---

## 4 · Caller pipelines

### 4.1 · Expression arm — 10x

```text
site-matched count matrix
        |
        +--> inferCNV
        +--> CopyKAT
        +--> SCEVAN
        |
        v
caller-specific continuous CNV score
        |
        v
held-out calibration -> P(malignant)
```

All primary runs include a disjoint same-site B reference.

#### inferCNV

The primary inferCNV run uses the same-site B reference. A **reference-free sensitivity condition** is retained separately and is never pooled with the referenced benchmark results.

#### CopyKAT

CopyKAT is run with:

```r
id.type = "E"
genome = "hg20"
cell.line = "no"
norm.cell.names = <same-site reference-B cells present in rawmat>
```

Because CopyKAT's native tumour/normal output is categorical, discrimination and calibration use a continuous CNV-burden score derived from `CNAmat` after the benchmark mask is applied.

#### SCEVAN

SCEVAN receives the same-site normal cells through `norm_cell` and uses:

```r
FIXED_NORMAL_CELLS = FALSE
```

A continuous SCEVAN score is retained through a pinned-version wrapper around the relevant preprocessing/classification steps rather than inferred from a binary label after the fact.

---

### 4.2 · Allele arm — Fluidigm C1 / Numbat

The SEQC2 full-length C1 dataset contains **146 one-cell libraries**:

- **80 HCC1395 (A)**;
- **66 HCC1395BL (B)**;
- SMART-Seq v4;
- 150 bp paired-end sequencing.

```text
C1 FASTQ
  |
  v
STAR -> one coordinate-sorted BAM per cell
  |
  v
pileup_and_phase.R --smartseq
  |
  +--> phased allele counts (df_allele)
  |
featureCounts matrix --------------------+
  |                                      |
  +--> count_mat                         |
                                         v
                                   Numbat inference
                                         |
                                         +--> p_cnv
                                         +--> p_cnv_x
                                         +--> p_cnv_y
                                         +--> clone assignments
                                         +--> phylogeny
```

#### Expression versus allele inputs

The deposited C1 **featureCounts matrices are the expression input** (`count_mat`). They are not allele counts.

The allele input (`df_allele`) is generated from the per-cell BAMs with Numbat's `pileup_and_phase.R` in `--smartseq` mode against the hg38 1000 Genomes SNP/phasing resources bundled with the Numbat container.

SMART-Seq mode uses ordered text files containing BAM paths and cell names; the orders must match. No 10x CB/UMI tags are assumed.

#### B-reference split

The 66 B cells are split deterministically into:

- **reference-B:** used only to construct `lambdas_ref`;
- **scored-B:** included with A cells as the diploid negative-control population.

No A cell contributes to `lambdas_ref`.

Reference-B cells are disjoint from scored-B cells. The role assignment and seed are written to the seed ledger.

#### Numbat settings

```r
genome = "hg38"
gamma  = 5
```

`gamma=5` is used because SMART-Seq is non-UMI and has noisier allele counts than 10x; Numbat recommends a smaller gamma for non-UMI protocols than the 10x default/recommendation of 20.

Numbat's native outputs are retained separately:

- `p_cnv` — posterior probability of belonging to an aneuploid clone using joint evidence;
- `p_cnv_x` — expression-only aneuploid-clone posterior;
- `p_cnv_y` — allele-only aneuploid-clone posterior.

These native posteriors are evaluated directly and may also enter the same held-out calibration framework used for the expression callers.

#### LOH handling

Preference order:

1. **Known clonal-LOH segments** via `segs_loh`, when a validated LOH track is available; or
2. **fixed consensus CNV segments** via `segs_consensus_fix`, when a complete validated bulk profile is used; otherwise
3. `call_clonal_loh=TRUE` as the fallback.

`segs_loh` and automatic clonal-LOH calling are never supplied together. Segment mode, source file, coordinate convention, and checksum are written to the run metadata.

#### Pilot gate

A **20-A + 20-scored-B** pilot is run before expansion. The pilot must satisfy all of the following before the full C1 arm proceeds:

- non-empty phased allele counts are produced;
- the selected cells are represented in Numbat output;
- scored-B cells are predominantly diploid / low-aneuploidy posterior;
- A cells show the expected broad tumour aneuploidy signal;
- output clone/posterior objects and the phylogeny are non-empty;
- genome/contig compatibility checks pass.

The raw C1 archive contains 146 libraries, but the **Numbat target population is A + scored-B**; reference-B cells remain disjoint and are used to construct `lambdas_ref` rather than being treated as scored targets.

> [!NOTE]
> C1 provides an orthogonal, allele-rich validation arm. Its sample size is much smaller than the 3,000-cell 10x contexts, so C1 admixture results are not interpreted as equal-powered replicas of the 10x 1%–80% purity surface.

---

### 4.3 · Physical-mixture source assignment

For real 5%/10% mixtures, source assignment is treated as a probabilistic classification problem rather than ordinary donor demultiplexing.

Fixed candidate sites are drawn from high-confidence SEQC2 somatic variants after excluding unsuitable loci. Per-cell evidence is aggregated across informative sites, with local CNA/LOH state considered where relevant. Pure A and pure B captures define the positive and negative evidence distributions. Cells below minimum allele support remain **unassigned**.

This source-assignment layer is an external check on the physical mixtures; it is not substituted for the known A/B labels in the in-silico benchmark.

---

## 5 · Calibration, uncertainty, and controls

### Calibration

Each caller contributes a continuous score. Calibration models are evaluated under **nested, source-cell-grouped cross-validation**:

1. Complete source cells are assigned to outer folds so the same biological cell cannot occur in both training and test data through different mixture contexts.
2. Calibration is fit only inside the training partition.
3. Platt, beta, and isotonic calibration are compared by inner cross-validation using log loss and/or Brier score.
4. The selected calibrator is evaluated once on the held-out outer fold.

Calibration truth is the **A/B source label**, not the bulk DNA segment track.

### Indeterminate interval

Two decision thresholds are fixed from training/calibration data before final test evaluation:

- a lower threshold chosen for the prespecified negative-predictive-value target;
- an upper threshold chosen for the prespecified positive-predictive-value target.

Cells between the two thresholds are reported as **indeterminate**, not forced into a binary class.

### Uncertainty

Because source cells recur across mixture fractions and replicates, ordinary row-level bootstrap resampling is invalid.

Uncertainty is estimated with at least **2,000 source-cell cluster-bootstrap replicates**, stratified by true A/B label and aggregated over mixture replicates. Every principal per-fraction metric receives a 95% interval.

A detectability threshold is declared only when the **prespecified performance criterion is met by its lower 95% confidence bound** and the replicate-consistency requirement is satisfied. Non-overlap of adjacent confidence intervals is not used as a surrogate significance test.

### Hard controls

The following are benchmark gates:

| Control | Expected result |
|---|---|
| **B-versus-B** | held-out scored-B against disjoint reference-B has a low false aneuploid rate |
| **Label shuffle** | discrimination collapses to chance (`AUC ≈ 0.5`) |
| **Reference-split stability** | changing the seeded B reference split does not qualitatively rewrite the dominant CNV pattern |
| **Build / ID QC** | no systematic coordinate mismatch or catastrophic gene loss |

chr6p, chr16q, and chrX are **masked diagnostic regions**, not false-positive ground truth, and therefore are not included as a hard “must be silent” control.

---

## 6 · Secondary analyses enabled by the benchmark

### A · Per-arm detectability floor

For each caller and DNA-truth arm/region, determine the lowest tumour fraction at which tumour-versus-normal discrimination reaches the prespecified threshold with adequate lower-confidence-bound support and replicate consistency.

A prespecified biological hypothesis is that detectability improves with larger-magnitude WGS copy-number deviation. Arms that never satisfy the criterion are reported as **not detectable within the tested purity range**.

### B · Same-cell expression/allele fusion

Cross-platform C1-to-10x fusion is prohibited because the cells are different.

An optional **same-cell 10x** fusion experiment may instead run Numbat on the same 10x BAMs and combine its calibrated output with the calibrated expression-caller scores. Any gain is judged on held-out Brier score, log loss, and ECE against the best single caller.

### C · Independent somatic-SNV classifier

The approximately 40,000 SEQC2 somatic SNVs provide an orthogonal tumour-specific signal. Fixed-site pileup is used as a high-precision classifier/anchor, with pure-B alternate evidence defining the error floor. Cells without enough informative reads remain unclassified.

### D · Clonal-versus-subclonal inference

For events with adequate informative-cell support, event prevalence among malignant cells is estimated after accounting for per-cell detection sensitivity and background error. Multiple-event testing is corrected with Benjamini-Hochberg. Low-support events remain indeterminate rather than receiving a forced clonality label.

### E · Calibration drift

Calibration error is measured separately at each tumour fraction. The **calibration floor** is reported separately from the discrimination/detectability floor because a caller may retain ranking performance while becoming systematically over- or under-confident.

### F · Context Sensitivity Index

Because the nested design can expose the same source cell to different surrounding tumour fractions, within-cell score variability is quantified across contexts. Large context-dependent shifts identify callers whose per-cell score depends strongly on cohort composition rather than only on the cell itself.

---

## 7 · Determinism and provenance

Every stochastic operation receives an explicit seed and is written to a machine-readable ledger.

Seeded operations include:

- reference-B / scored-B partition;
- mixture replicate construction;
- nested A/B orderings;
- C1 pilot subsampling;
- bootstrap resampling;
- outer and inner calibration folds;
- label-shuffle controls.

The run ledger records, at minimum:

```text
context_id
site
platform
fraction
replicate
operation
seed
input_file
input_checksum
reference_build
caller
caller_version
container_digest
```

Rerunning from the recorded inputs, checksums, versions, and seeds must reproduce the mixture composition and benchmark scores.

---

## 8 · Required outputs

### Per-cell outputs

- true source label where known;
- caller raw score;
- calibrated `P(malignant)`;
- final decision: normal / indeterminate / malignant;
- mixture context and replicate;
- reference split identifier.

### Genomic outputs

- arm/segment call after masking;
- truth state;
- direction concordance;
- informative-gene / informative-SNP support;
- caller-specific uncertainty measures where available.

### Numbat outputs

- `p_cnv`;
- `p_cnv_x`;
- `p_cnv_y`;
- clone assignments and clone posteriors;
- consensus CNV segments;
- expression, allele, and joint posteriors;
- phylogeny / mutation graph.

### Benchmark summaries

- ROC-AUC and PR-AUC;
- Brier score;
- log loss;
- expected calibration error (ECE);
- sensitivity / specificity at prespecified operating points;
- indeterminate fraction;
- per-arm detectability floor;
- calibration drift by purity;
- context-sensitivity metrics;
- cluster-bootstrap 95% intervals.

---

## 9 · References

1. **Chen W, et al.** A multicenter study benchmarking single-cell RNA sequencing technologies using reference samples. *Nature Biotechnology* (2021). https://doi.org/10.1038/s41587-020-00748-9
2. **Chen X, et al.** A multi-center cross-platform single-cell RNA sequencing reference dataset. *Scientific Data* **8**, 39 (2021). https://doi.org/10.1038/s41597-021-00809-x
3. **Fang LT, et al.** Establishing community reference samples, data and call sets for benchmarking cancer mutation detection using whole-genome sequencing. *Nature Biotechnology* **39**, 1151–1160 (2021). https://doi.org/10.1038/s41587-021-00993-6
4. **Talsania K, et al.** Structural variant analysis of a cancer reference cell line sample using multiple sequencing technologies. *Genome Biology* **23**, 255 (2022). https://doi.org/10.1186/s13059-022-02816-6
5. **Masood D, et al.** Evaluation of somatic copy number variation detection by NGS technologies and bioinformatics tools on a hyper-diploid cancer genome. *Genome Biology* **25**, 163 (2024). https://doi.org/10.1186/s13059-024-03294-8
6. **Gao T, et al.** Haplotype-aware analysis of somatic copy number variations from single-cell transcriptomes. *Nature Biotechnology* (2023). https://doi.org/10.1038/s41587-022-01468-y
7. **Numbat user guide.** SMART-Seq preprocessing, `gamma`, clonal-LOH handling, and fixed consensus CNV inputs. https://github.com/kharchenkolab/numbat/blob/main/vignettes/numbat.Rmd

---

<div align="center">

**Repository:** `collaborativebioinformatics/scOnco-scorer`

*Reproducible by construction · calibrated by source-cell truth · validated against orthogonal DNA CNV/LOH evidence*

</div>
