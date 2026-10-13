# Identifying Direct Transcriptional Targets of Mef2d in the Mouse Retina

Final project for **BCH 419** (Bioinformatics). This project integrates **ChIP-seq** and **RNA-seq** data to find genes that are both *bound* by the transcription factor **Mef2d** and *differentially expressed* when Mef2d is knocked out. Genes meeting both criteria are strong candidates for direct Mef2d targets.

## Overview

Mef2d (Myocyte Enhancer Factor 2D) is a MADS-box transcription factor. Binding sites alone don't prove a gene is regulated, and expression changes alone don't prove the effect is direct. Intersecting the two narrows the list down to genes where Mef2d binds near the promoter *and* loss of Mef2d changes expression.

The pipeline has three parts:

1. **ChIP-seq peak processing (Python + bedtools):** Reformat MACS peak summits, keep reproducible peaks found in both wild-type replicates, and remove peaks that also appear in the knockout control.
2. **Peak-to-gene annotation (GREAT + Python):** Assign peaks to nearby genes with GREAT and keep genes whose nearest Mef2d peak is within ±2 kb.
3. **Differential expression and integration (R):** Quantify expression with Salmon, test WT vs. KO with DESeq2, run GO enrichment with clusterProfiler, and intersect significant genes with the Mef2d-bound gene list.

## Key Results

| Step | Count |
|---|---|
| Peaks shared between WT replicates | 1,711 |
| High-confidence peaks (absent from KO) | 1,464 |
| Genes with a Mef2d peak within ±2 kb | 245 |
| Significantly differentially expressed genes (padj ≤ 0.01, \|log2FC\| ≥ 0.5) | 401 |
| **Bound *and* differentially expressed** | **24** |

Of the 24 putative direct targets, **23 are downregulated** in the knockout, consistent with Mef2d acting mainly as a transcriptional **activator**. Many of the top hits are photoreceptor and phototransduction genes, including *Arr3*, *Sag*, *Gngt1*, *Rd3l*, *Rgs20*, and *Wdr17*. *Arr3* shows the strongest effect (log2FC ≈ −8.6). The only upregulated gene is *Gins2*.

The full table is in [`Mef2d_final.csv`](Mef2d_final.csv).

## Repository Structure

```
.
├── README.md
├── Required Material/
│   ├── Sumrall_Aidan_fp_Part1.py     # Filter/reformat MACS summit BED files
│   ├── Sumrall_Aidan_fp_Part2.py     # Parse GREAT output -> Mef2d target genes
│   ├── Sumrall_Aidan_fp.R            # DESeq2, GO enrichment, ChIP/RNA integration
│   └── Mef2d_truePeaks.bed           # Final high-confidence Mef2d peaks
└── Miscellaneous/
    ├── Mef2d_R1_summits_great.bed    # WT replicate 1 peaks (summit ±500 bp)
    ├── Mef2d_R2_summits_great.bed    # WT replicate 2 peaks
    ├── Mef2d_KO_summits_great.bed    # Knockout control peaks
    ├── Mef2d_common_summits.bed      # Peaks shared by R1 and R2
    ├── greatExportAll.tsv.txt        # GREAT region-gene associations
    ├── Mef2d_targets.csv             # Genes within ±2 kb of a Mef2d peak
    ├── WT-KO-All.csv                 # Full DESeq2 results
    ├── WT-KO-All_PreMef2d.csv        # Full DESeq2 results (earlier run, identical)
    ├── W_WT-KO-Sig_PreMef2d.csv      # Significant DE genes (before ChIP integration)
    ├── W_WT-KO-Sig.csv               # Significant DE genes after integration
    ├── Mef2d_final.csv               # Final list of putative direct targets
    └── References.txt                # Documentation and tutorials used
```

## Pipeline Details

### Part 1: Peak processing (`fp_Part1.py`)

Reads a MACS `*_summits.bed` file and:

- keeps only standard chromosomes (`chr1`–`chr19`, `chrX`)
- extends each summit by 500 bp on each side to create ~1 kb regions
- rounds peak scores to integers so the file is compatible with GREAT

The script is run once per sample (R1, R2, KO) by changing the input/output file names.

The processed files were then combined with **bedtools**: peaks overlapping between R1 and R2 were kept (`Mef2d_common_summits.bed`), and any of those overlapping a KO peak were removed to give `Mef2d_truePeaks.bed`. Commands along these lines:

```bash
bedtools intersect -u -a Mef2d_R1_summits_great.bed -b Mef2d_R2_summits_great.bed > Mef2d_common_summits.bed
bedtools intersect -v -a Mef2d_common_summits.bed  -b Mef2d_KO_summits_great.bed > Mef2d_truePeaks.bed
```

### Part 2: Gene annotation (`fp_Part2.py`)

`Mef2d_truePeaks.bed` was submitted to [GREAT](http://great.stanford.edu/) (v4.0.4, mm10, basal + extension rule). The script parses the exported region-gene associations and:

- keeps genes within **±2,000 bp** of a peak
- if a gene is linked to multiple peaks, keeps the **closest** one
- writes `Mef2d_targets.csv` (`gene, peak, distance`), sorted alphabetically

### Part 3: Differential expression and integration (`fp.R`)

- Imports Salmon `quant.sf` files for 2 WT and 2 KO samples with **tximport**, using Ensembl GRCm38 (release 100) annotation
- Filters genes with fewer than 10 total reads and runs **DESeq2** (KO vs. WT), with **apeglm** log-fold-change shrinkage
- Calls genes significant at **padj ≤ 0.01** and **|log2FC| ≥ 0.5**
- Intersects significant genes with `Mef2d_targets.csv`, attaching peak name and distance, and writes `Mef2d_final.csv`
- Splits results into up- and downregulated genes and runs **GO Biological Process** enrichment with **clusterProfiler** (q < 0.1), using all tested genes as the background

## Requirements

**Python 3** (standard library only)

**R / Bioconductor packages:**
`ensembldb`, `tximport`, `DESeq2`, `vsn`, `apeglm`, `clusterProfiler`, `enrichplot`, `org.Mm.eg.db`, `ggplot2`, `dplyr`, `scales`, `gplots`

**External tools:** MACS (peak calling), bedtools, Salmon, GREAT (web tool)

## Running It

The scripts were written to run on the course's HPC cluster and expect specific relative paths:

- `Part1.py` reads `Mef2d_KO/Mef2d_KO_summits.bed`; edit the input/output names for each sample.
- `Part2.py` expects `greatExportAll.tsv.txt` in the working directory.
- `fp.R` loads the Ensembl database from a cluster path (`/projects/academic/courses/bch519/...`) and Salmon output from `./WT1`, `./WT2`, `./KO1`, `./KO2`. Update these paths to run it elsewhere.

Raw sequencing reads and Salmon quantification files are not included in this repository because of their size.

```bash
python3 "Required Material/Sumrall_Aidan_fp_Part1.py"
# bedtools steps + GREAT submission
python3 "Required Material/Sumrall_Aidan_fp_Part2.py"
Rscript "Required Material/Sumrall_Aidan_fp.R"
```

## Notes

- In `fp.R`, the significant-gene table `sig` is overwritten by the Mef2d intersection before being saved, so `W_WT-KO-Sig.csv` contains the same 24 genes as `Mef2d_final.csv`. The full significant list from before integration is in `W_WT-KO-Sig_PreMef2d.csv`.
- The R script saves to `W:WT-KO-Sig.csv`; the colon was converted to an underscore when the file was copied to Windows.
- With only two biological replicates per condition, these results are exploratory and best treated as candidates for follow-up validation.

## References

A full list of documentation and tutorials used is in [`Miscellaneous/References.txt`](Miscellaneous/References.txt). Key tools:

- Love, Huber & Anders (2014). DESeq2. *Genome Biology*.
- Patro et al. (2017). Salmon. *Nature Methods*.
- Soneson, Love & Robinson (2015). tximport. *F1000Research*.
- Wu et al. (2021). clusterProfiler 4.0. *The Innovation*.
- McLean et al. (2010). GREAT. *Nature Biotechnology*.
- Quinlan & Hall (2010). BEDTools. *Bioinformatics*.
- Zhu, Ibrahim & Love (2019). apeglm. *Bioinformatics*.

## Author

Aidan Sumrall
