# GeneTrack Studio

A Shiny GUI for building publication-quality genomic track figures using [plotgardener](https://phanstiellab.github.io/plotgardener/).

## Requirements

R packages:
- shiny
- plotgardener
- GenomicFeatures, GenomicRanges, IRanges, S4Vectors
- TxDb.Hsapiens.UCSC.hg38.refGene
- org.Hs.eg.db
- AnnotationDbi, BiocGenerics

## Running

```r
# From the genetrack-studio directory:
shiny::runApp(".")

# Or from R anywhere:
shiny::runApp("/path/to/genetrack-studio")
```

## Supported Track Types (Phase 1)

- **Signal** — bigWig files (CUT&RUN, ChIP-seq, ATAC-seq, RNA-seq coverage)
- **Hi-C** — .hic contact matrices (Micro-C, Hi-C) rendered as triangles
- **Ranges** — BED files (CpG islands, peaks, any interval data)
- **Genes** — Custom IGV-style gene annotation from refGene

## Usage

1. Set your genomic region (chrom, start, end)
2. Add tracks using the + buttons
3. For each track, provide the file path and adjust settings
4. Click **Render Figure**
5. Download the PNG

## Roadmap

- Phase 2: Track reordering, save/load configs, live preview
- Phase 3: More track types (loops, GWAS, bedGraph), multi-region panels
- Phase 4: Docker packaging, Electron desktop wrapper
