# utils.R — shared helpers

id_to_symbol <- function(entrez_ids) {
  syms <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = as.character(entrez_ids),
    column  = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  syms <- unname(syms)
  syms[is.na(syms)] <- as.character(entrez_ids)[is.na(syms)]
  syms
}

make_mapper <- function(x, width, start, end) {
  function(pos) x + (pos - start) / (end - start) * width
}

pack_rows <- function(starts, ends) {
  ord <- order(starts, ends)
  starts <- starts[ord]; ends <- ends[ord]; idx <- ord
  row_end <- numeric(0); row_of <- integer(length(starts))
  for (i in seq_along(starts)) {
    placed <- FALSE
    for (r in seq_along(row_end)) {
      if (starts[i] > row_end[r]) {
        row_of[i] <- r; row_end[r] <- ends[i]; placed <- TRUE; break
      }
    }
    if (!placed) { row_end <- c(row_end, ends[i]); row_of[i] <- length(row_end) }
  }
  out <- integer(length(row_of)); out[idx] <- row_of; out
}

get_gene_representative_models <- function(txdb, chrom, start, end) {
  view_gr <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start, end))
  genes_gr <- suppressMessages(GenomicFeatures::genes(txdb))
  genes_gr <- IRanges::subsetByOverlaps(genes_gr, view_gr)
  if (length(genes_gr) == 0)
    return(list(genes = genes_gr, cds = GenomicRanges::GRanges(), utr = GenomicRanges::GRanges()))

  tx_by_gene <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
  gene_ids <- intersect(names(genes_gr), names(tx_by_gene))
  genes_gr <- genes_gr[gene_ids]; tx_by_gene <- tx_by_gene[gene_ids]
  ex_by_tx <- GenomicFeatures::exonsBy(txdb, by = "tx")
  cds_by_tx <- GenomicFeatures::cdsBy(txdb, by = "tx")

  rep_tx_id <- character(0)
  for (gid in gene_ids) {
    txs <- tx_by_gene[[gid]]
    if (length(txs) == 0) next
    tx_ids <- as.character(S4Vectors::mcols(txs)$tx_id)
    tx_ids <- tx_ids[tx_ids %in% names(ex_by_tx)]
    if (length(tx_ids) == 0) next
    exon_bp <- sapply(tx_ids, function(tid) sum(BiocGenerics::width(GenomicRanges::reduce(ex_by_tx[[tid]]))))
    cds_bp <- sapply(tx_ids, function(tid) {
      if (tid %in% names(cds_by_tx)) sum(BiocGenerics::width(GenomicRanges::reduce(cds_by_tx[[tid]]))) else 0
    })
    best <- tx_ids[order(exon_bp, cds_bp, decreasing = TRUE)][1]
    rep_tx_id[gid] <- best
  }

  cds_gr <- GenomicRanges::GRanges(); utr_gr <- GenomicRanges::GRanges()
  for (gid in names(rep_tx_id)) {
    tid <- rep_tx_id[[gid]]
    ex_g <- GenomicRanges::reduce(ex_by_tx[[tid]])
    if (length(ex_g) == 0) next
    if (tid %in% names(cds_by_tx)) {
      cds_g <- GenomicRanges::reduce(cds_by_tx[[tid]])
      utr_g <- GenomicRanges::reduce(GenomicRanges::setdiff(ex_g, cds_g))
    } else { cds_g <- GenomicRanges::GRanges(); utr_g <- ex_g }
    if (length(cds_g) > 0) { cds_g$gene_id <- gid; cds_gr <- c(cds_gr, cds_g) }
    if (length(utr_g) > 0) { utr_g$gene_id <- gid; utr_gr <- c(utr_gr, utr_g) }
    gi <- match(gid, names(genes_gr))
    if (!is.na(gi))
      GenomicRanges::ranges(genes_gr)[gi] <- IRanges::IRanges(min(BiocGenerics::start(ex_g)), max(BiocGenerics::end(ex_g)))
  }
  genes_gr <- IRanges::subsetByOverlaps(genes_gr, view_gr)
  cds_gr <- IRanges::subsetByOverlaps(cds_gr, view_gr)
  utr_gr <- IRanges::subsetByOverlaps(utr_gr, view_gr)
  list(genes = genes_gr, cds = cds_gr, utr = utr_gr)
}
