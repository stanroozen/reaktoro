lib <- file.path(Sys.getenv('USERPROFILE'), 'Documents', 'R', 'win-library', '4.5')
.libPaths(c(lib, .libPaths()))
library(CHNOSZ)

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep('^--file=', args_all)]
script_path <- if (length(script_arg) > 0) sub('^--file=', '', script_arg[1]) else getwd()
OUTDIR <- normalizePath(dirname(script_path), winslash = '/', mustWork = FALSE)

env_int <- function(name, default) {
  v <- Sys.getenv(name, unset = as.character(default))
  n <- suppressWarnings(as.integer(v))
  if (is.na(n)) default else n
}

n_pH_p <- max(20, env_int('CHNOSZ_POURBAIX_NPH', 120))
n_Eh_p <- max(20, env_int('CHNOSZ_POURBAIX_NEH', 100))
n_pH_m <- max(20, env_int('CHNOSZ_MOSAIC_NPH', 100))
n_Eh_m <- max(20, env_int('CHNOSZ_MOSAIC_NEH', 80))

# ------------------------------------------------------------------
# CHNOSZ Pourbaix equivalent (matched ranges)
# ------------------------------------------------------------------
reset()
OBIGT(no.organics = TRUE)

basis(c('Fe+2', 'H2O', 'H+', 'e-'))
species_p <- c(
  'Fe+2', 'Fe+3', 'FeO+', 'FeO2-',
  'FeOH+', 'FeOH+2', 'HFeO2', 'HFeO2-',
  'goethite', 'hematite', 'iron', 'magnetite'
)
species(species_p)

a_p <- affinity(
  pH = c(-2.0, 16.0, n_pH_p),
  Eh = c(-2.0, 2.0, n_Eh_p),
  T = 25,
  P = 1
)
d_p <- diagram(a_p, plot.it = FALSE)

out_p <- file.path(OUTDIR, 'CHNOSZ_Pourbaix_Fe.png')
png(out_p, width = 1600, height = 1200, res = 180)
par(mar = c(4.2, 4.2, 2.8, 1.2))
diagram(d_p, lwd = 1.8, xlab = 'pH', ylab = 'Eh')
try(water.lines(lty = 4, lwd = 1.2, col = 'black'), silent = TRUE)
abline(h = 0, v = 7, lty = 3, lwd = 1.0, col = 'gray45')
legend(
  'bottomright',
  legend = c(
    'Predominance boundary (solid)',
    'Water-stability lines (dash-dot)',
    'Reference guides: Eh = 0, pH = 7'
  ),
  lty = c(1, 4, 3),
  lwd = c(1.8, 1.2, 1.0),
  col = c('black', 'black', 'gray45'),
  bg = 'white',
  cex = 0.85
)
title(main = 'CHNOSZ Fe-O-H Pourbaix (Matched Setup)')
dev.off()

# ------------------------------------------------------------------
# CHNOSZ Mosaic equivalent (matched ranges)
# Uses Fe-S-C-O-H basis from diagnostics script conventions.
# ------------------------------------------------------------------
reset()
OBIGT(no.organics = TRUE)

basis(c('FeO', 'SO4-2', 'CO3-2', 'H2O', 'H+', 'e-'))
basis('SO4-2', -6)
basis('CO3-2', 0)

species_m <- c('Fe+2', 'Fe+3', 'HFeO2-', 'pyrite', 'pyrrhotite', 'siderite', 'hematite', 'magnetite')
species(species_m)

a_m <- affinity(
  pH = c(0.0, 14.0, n_pH_m),
  Eh = c(-1.0, 1.0, n_Eh_m),
  T = 25,
  P = 1
)
d_m <- diagram(a_m, plot.it = FALSE)

out_m <- file.path(OUTDIR, 'CHNOSZ_Mosaic_Fe.png')
png(out_m, width = 1600, height = 1200, res = 180)
par(mar = c(4.2, 4.2, 2.8, 1.2))
diagram(d_m, lwd = 1.8, xlab = 'pH', ylab = 'Eh')
try(water.lines(lty = 4, lwd = 1.2, col = 'black'), silent = TRUE)
abline(h = 0, v = 7, lty = 3, lwd = 1.0, col = 'gray45')
legend(
  'bottomright',
  legend = c(
    'Predominance boundary (solid)',
    'Water-stability lines (dash-dot)',
    'Reference guides: Eh = 0, pH = 7'
  ),
  lty = c(1, 4, 3),
  lwd = c(1.8, 1.2, 1.0),
  col = c('black', 'black', 'gray45'),
  bg = 'white',
  cex = 0.85
)
title(main = 'CHNOSZ Fe-S-C-O-H Mosaic (Matched Setup)')
dev.off()

out_txt <- file.path(OUTDIR, 'CHNOSZ_matched_setup.txt')
writeLines(c(
  'Generated CHNOSZ matched equivalents',
  sprintf('- Pourbaix ranges: pH [-2,16], Eh [-2,2], %dx%d points, T=25C, P=1 bar', n_pH_p, n_Eh_p),
  sprintf('- Mosaic ranges: pH [0,14], Eh [-1,1], %dx%d points, T=25C, P=1 bar', n_pH_m, n_Eh_m),
  '- Mosaic basis settings: loga(SO4-2) = -6, loga(CO3-2) = 0',
  '- Line semantics: solid=predominance boundaries, dash-dot=water-stability lines, dotted=Eh0/pH7 guides',
  sprintf('- Output: %s', normalizePath(out_p, winslash = '/', mustWork = FALSE)),
  sprintf('- Output: %s', normalizePath(out_m, winslash = '/', mustWork = FALSE))
), out_txt)

cat('Wrote:', out_p, '\n')
cat('Wrote:', out_m, '\n')
cat('Wrote:', out_txt, '\n')
