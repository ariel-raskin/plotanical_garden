# render_figure.R — main plotgardener orchestrator

source("R/utils.R", local = TRUE)

draw_elbow_arrow <- function(x_tss, y_base, dir, col = "black", lwd = 0.28,
                              elbow_h = 0.05, elbow_w = 0.12, arrow_mm = 0.75,
                              y_offset = 0.00, page_height) {
  y0 <- y_base + y_offset
  plotgardener::plotSegments(x0 = x_tss, y0 = y0, x1 = x_tss, y1 = y0 + elbow_h,
                             linecolor = col, lwd = lwd)
  if (dir == "+") { x_end <- x_tss + elbow_w }
  else if (dir == "-") { x_end <- x_tss - elbow_w }
  else return(invisible(NULL))

  plotgardener::plotSegments(x0 = x_tss, y0 = y0 + elbow_h, x1 = x_end,
                             y1 = y0 + elbow_h, linecolor = col, lwd = lwd)
  a_in <- arrow_mm / 25.4; h_in <- a_in * 0.70
  y_tip_top <- y0 + elbow_h
  to_bottom <- function(y_top) page_height - y_top
  if (dir == "+") {
    grid::grid.polygon(
      x = grid::unit(c(x_end, x_end - a_in, x_end - a_in), "inches"),
      y = grid::unit(c(to_bottom(y_tip_top), to_bottom(y_tip_top + h_in), to_bottom(y_tip_top - h_in)), "inches"),
      gp = grid::gpar(col = col, fill = col, lwd = lwd))
  } else {
    grid::grid.polygon(
      x = grid::unit(c(x_end, x_end + a_in, x_end + a_in), "inches"),
      y = grid::unit(c(to_bottom(y_tip_top), to_bottom(y_tip_top + h_in), to_bottom(y_tip_top - h_in)), "inches"),
      gp = grid::gpar(col = col, fill = col, lwd = lwd))
  }
  invisible(NULL)
}

draw_gene_track <- function(chrom, start, end, x, y, width, height,
                             txdb, page_height,
                             exon_blue = "#0B3D91", cds_h = 0.115, utr_h = 0.055,
                             baseline_col = "black", baseline_lwd = 0.22,
                             label_fs = 7, label_col = "black", max_rows = 6,
                             min_label_bp = 12000, label_y_offset = -0.15,
                             arrows = TRUE, arrow_head_mm = 0.75, arrow_y_offset = 0.07) {
  models <- get_gene_representative_models(txdb, chrom, start, end)
  genes_gr <- models$genes; cds_gr <- models$cds; utr_gr <- models$utr
  if (length(genes_gr) == 0) {
    plotgardener::plotText("No genes in region", x = x, y = y + height/2,
                           just = c("left","center"), fontsize = 8)
    return(invisible(NULL))
  }
  gene_ids <- names(genes_gr)
  syms <- id_to_symbol(gene_ids)
  g_st <- pmax(GenomicRanges::start(genes_gr), start)
  g_en <- pmin(GenomicRanges::end(genes_gr), end)
  rows <- pack_rows(g_st, g_en)
  nrows <- min(max(rows), max_rows)
  map_x <- make_mapper(x, width, start, end)
  row_spacing <- height / nrows
  ord <- order(g_st, g_en)
  for (k in ord) {
    r <- rows[k]; if (r > nrows) next
    y_mid <- y + (r - 0.5) * row_spacing
    gid <- gene_ids[k]; sym <- syms[k]
    strand_k <- as.character(GenomicRanges::strand(genes_gr[k]))
    gs <- g_st[k]; ge <- g_en[k]
    plotgardener::plotSegments(x0 = map_x(gs), y0 = y_mid, x1 = map_x(ge), y1 = y_mid,
                               linecolor = baseline_col, lwd = baseline_lwd)
    if (arrows && strand_k %in% c("+","-")) {
      x_tss <- if (strand_k == "+") map_x(gs) else map_x(ge)
      draw_elbow_arrow(x_tss, y_mid, strand_k, baseline_col, baseline_lwd,
                       arrow_mm = arrow_head_mm, y_offset = arrow_y_offset,
                       page_height = page_height)
    }
    utr_i <- utr_gr[utr_gr$gene_id == gid]
    if (length(utr_i) > 0) for (j in seq_along(utr_i)) {
      u <- utr_i[j]; us <- max(GenomicRanges::start(u), start); ue <- min(GenomicRanges::end(u), end)
      plotgardener::plotRect(x = map_x(us), y = y_mid - utr_h/2,
        width = max(0.0005, map_x(ue) - map_x(us)), height = utr_h,
        fill = exon_blue, border = NA, lwd = 0, just = c("left","top"))
    }
    cds_i <- cds_gr[cds_gr$gene_id == gid]
    if (length(cds_i) > 0) for (j in seq_along(cds_i)) {
      cc <- cds_i[j]; cs <- max(GenomicRanges::start(cc), start); ce <- min(GenomicRanges::end(cc), end)
      plotgardener::plotRect(x = map_x(cs), y = y_mid - cds_h/2,
        width = max(0.0005, map_x(ce) - map_x(cs)), height = cds_h,
        fill = exon_blue, border = NA, lwd = 0, just = c("left","top"))
    }
    if ((ge - gs) >= min_label_bp) {
      plotgardener::plotText(label = sym, x = map_x(gs), y = y_mid + label_y_offset,
        just = c("left","top"), fontsize = label_fs, fontcolor = label_col)
    }
  }
  invisible(NULL)
}

#' Render the full figure to a temp PNG and return the path
#' @param tracks list of track configs (list with type, file, label, color, height, ...)
#' @param region list(chrom, start, end)
#' @param page list(width, height, dpi)
#' @param txdb a TxDb object
#' @return path to rendered PNG
render_figure <- function(tracks, region, page, txdb) {
  outfile <- tempfile(fileext = ".png")
  PAGE_W <- page$width %||% 9
  PAGE_H <- page$height %||% 6
  DPI    <- page$dpi %||% 600

  png_args <- list(filename = outfile, width = PAGE_W, height = PAGE_H,
                   units = "in", res = DPI)
  if (capabilities("cairo")) png_args$type <- "cairo-png"
  do.call(png, png_args)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)

  plotgardener::pageCreate(width = PAGE_W, height = PAGE_H,
                            default.units = "inches", showGuides = FALSE)

  # Layout constants
  TRACK_W   <- PAGE_W * 0.58
  X_RIGHT   <- PAGE_W - 0.50
  LABEL_GAP <- 1.05
  w  <- TRACK_W
  x0 <- X_RIGHT - w
  x_lab <- x0 - LABEL_GAP

  # Compute vertical positions: stack tracks top-down with gaps
  y_cursor <- 0.55
  GAP <- 0.10

  for (tk in tracks) {
    h <- tk$height %||% 0.75
    tk$y_pos  <- y_cursor
    tk$h_calc <- h

    tryCatch({
      if (tk$type == "hic") {
        hic <- plotgardener::plotHicTriangle(
          data = tk$file, chrom = region$chrom,
          chromstart = region$start, chromend = region$end,
          assembly = "hg38", norm = tk$norm %||% "KR",
          x = x0, y = y_cursor, width = w, height = h,
          just = c("left","top"))
        plotgardener::annoGenomeLabel(plot = hic, scale = "Kb",
          x = x0, y = y_cursor - 0.18, just = c("left","top"))
        plotgardener::annoHeatmapLegend(plot = hic,
          x = x0 + w + 0.20, y = y_cursor, width = 0.15, height = h - 0.35,
          just = c("left","top"))

      } else if (tk$type == "signal") {
        sig_range <- if (!is.null(tk$range)) tk$range else {
          plotgardener::calcSignalRange(data = tk$file, chrom = region$chrom,
            chromstart = region$start, chromend = region$end, assembly = "hg38")
        }
        col <- tk$color %||% "#1f77b4"
        sig <- plotgardener::plotSignal(
          data = tk$file, chrom = region$chrom,
          chromstart = region$start, chromend = region$end,
          x = x0, y = y_cursor, width = w, height = h,
          just = c("left","top"), range = sig_range,
          fill = col, linecolor = col, label = NULL)
        ymax_lbl <- format(round(sig_range[2], 2), trim = TRUE)
        plotgardener::annoYaxis(plot = sig, at = sig_range[2], label = ymax_lbl,
          axisLine = TRUE, fontsize = 4, main = TRUE)

      } else if (tk$type == "ranges") {
        bed <- utils::read.table(tk$file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
        colnames(bed)[1:3] <- c("chrom","start","end")
        plotgardener::plotRanges(data = bed, chrom = region$chrom,
          chromstart = region$start, chromend = region$end,
          x = x0, y = y_cursor, width = w, height = h,
          just = c("left","top"))

      } else if (tk$type == "genes") {
        draw_gene_track(region$chrom, region$start, region$end,
          x = x0, y = y_cursor, width = w, height = h,
          txdb = txdb, page_height = PAGE_H,
          exon_blue = tk$color %||% "#0B3D91",
          label_fs = tk$label_fs %||% 7, max_rows = tk$max_rows %||% 5)
      }
    }, error = function(e) {
      plotgardener::plotText(label = paste("Error:", conditionMessage(e)),
        x = x0, y = y_cursor + h/2, just = c("left","center"), fontsize = 7,
        fontcolor = "red")
    })

    # Left label
    plotgardener::plotText(label = tk$label %||% tk$type,
      x = x_lab, y = y_cursor + h/2,
      just = c("right","center"), fontsize = 9)

    y_cursor <- y_cursor + h + GAP
  }

  dev.off()
  on.exit(NULL)  # cancel the on.exit since we already closed
  outfile
}
