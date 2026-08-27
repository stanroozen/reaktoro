lib <- file.path(Sys.getenv('USERPROFILE'), 'Documents', 'R', 'win-library', '4.5')
.libPaths(c(lib, .libPaths()))
library(CHNOSZ)

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep('^--file=', args_all)]
script_path <- sub('^--file=', '', script_arg[1])
OUTDIR <- normalizePath(dirname(script_path), winslash = '/', mustWork = TRUE)

nearest_idx <- function(vals, x) {
  which.min(abs(vals - x))
}

write_case <- function(case_name, pH_vals, Eh_vals, basis_species, species_names, out_file, T = 25, P = 1, basis_logas = NULL) {
  reset()
  OBIGT(no.organics = TRUE)

  basis(basis_species)
  if (!is.null(basis_logas)) {
    for (nm in names(basis_logas)) basis(nm, basis_logas[[nm]])
  }
  species(species_names)

  a <- affinity(pH = c(min(pH_vals), max(pH_vals), length(pH_vals)),
                Eh = c(min(Eh_vals), max(Eh_vals), length(Eh_vals)),
                T = T, P = P)
  d <- diagram(a, plot.it = FALSE)

  # a$values is list of arrays, one per species
  vals_arr <- lapply(a$values, function(x) {
    d <- dim(x)
    if (is.null(d)) {
      matrix(x, nrow = length(a$vals[[1]]), ncol = length(a$vals[[2]]))
    } else if (length(d) == 3) {
      x[, , 1]
    } else if (length(d) == 2) {
      x
    } else {
      stop('Unsupported dimensions in affinity values')
    }
  })

  rows <- list()
  k <- 1
  for (pH in pH_vals) {
    ix <- nearest_idx(a$vals[[1]], pH)
    for (Eh in Eh_vals) {
      iy <- nearest_idx(a$vals[[2]], Eh)

      spvals <- sapply(vals_arr, function(v) v[ix, iy])
      o <- order(spvals, decreasing = TRUE)

      top1_i <- o[1]
      top2_i <- if (length(o) > 1) o[2] else NA
      pred_idx <- top1_i
      pred_name <- species_names[top1_i]

      rows[[k]] <- data.frame(
        case = case_name,
        pH = pH,
        Eh_V = Eh,
        pred_idx = pred_idx,
        pred_species_chnosz = pred_name,
        top1_species_chnosz = species_names[top1_i],
        top1_affinity_A2303RT = spvals[top1_i],
        top2_species_chnosz = if (!is.na(top2_i)) species_names[top2_i] else '',
        top2_affinity_A2303RT = if (!is.na(top2_i)) spvals[top2_i] else NA_real_,
        stringsAsFactors = FALSE
      )
      k <- k + 1
    }
  }

  out <- do.call(rbind, rows)
  write.csv(out, out_file, row.names = FALSE)
  cat('Wrote:', out_file, '\n')
}

# --- Fe Pourbaix (CHNOSZ demo-like) ---
pH_p <- c(-2, 1, 4, 7, 10, 13, 16)
Eh_p <- c(-2, -1, 0, 1, 2)
basis_p <- c('Fe+2', 'H2O', 'H+', 'e-')
species_p <- c('Fe+2','Fe+3','FeOH+','FeOH+2','HFeO2-','HFeO2','FeO+','FeO2-',
               'hematite','magnetite','goethite','iron')
out_p <- file.path(OUTDIR, 'diagnostics_chnosz_pourbaix_points.csv')
write_case('Fe_Pourbaix', pH_p, Eh_p, basis_p, species_p, out_p)

# --- Fe Mosaic (CHNOSZ demo-like, simplified) ---
pH_m <- c(0, 2, 5, 8, 11, 14)
Eh_m <- c(-1, -0.5, 0, 0.5, 1)
basis_m <- c('FeO', 'SO4-2', 'CO3-2', 'H2O', 'H+', 'e-')
species_m <- c('Fe+2','Fe+3','HFeO2-','pyrite','pyrrhotite','siderite','hematite','magnetite')
out_m <- file.path(OUTDIR, 'diagnostics_chnosz_mosaic_points.csv')
write_case('Fe_Mosaic', pH_m, Eh_m, basis_m, species_m, out_m, basis_logas = list('SO4-2' = -6, 'CO3-2' = 0))
