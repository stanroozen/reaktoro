lib <- file.path(Sys.getenv('USERPROFILE'), 'Documents', 'R', 'win-library', '4.5')
.libPaths(c(lib, .libPaths()))
library(CHNOSZ)

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- args_all[grep('^--file=', args_all)]
script_path <- sub('^--file=', '', script_arg[1])
OUTDIR <- normalizePath(dirname(script_path), winslash = '/', mustWork = TRUE)

reset()
OBIGT(no.organics = TRUE)

# ── Open-system basis definition ──────────────────────────────────────────────
# Setting Fe+2 as a basis species with log10(a) = 0 (a = 1) is the CHNOSZ way
# of specifying an open system: iron is not conserved; instead its chemical
# potential (activity) is fixed.  This is equivalent to Reaktoro's
#   specs.lgActivity("Fe+2")  +  conditions.lgActivity("Fe+2", 0.0)
LOGA_FE2 <- 0  # log10(a(Fe+2)) = 0  → a = 1 (CHNOSZ default for basis species)

# Strict test: only species/minerals shared with Reaktoro list.
basis(c('Fe+2', 'H2O', 'H+', 'e-'))
basis('Fe+2', LOGA_FE2)   # explicitly set log10(a(Fe+2)) = 0 (matches Reaktoro constraint)
species_shared <- c('Fe+2','Fe+3','FeO+','FeO2-','FeOH+','FeOH+2','HFeO2','HFeO2-',
			  'goethite','hematite','iron','magnetite')
species(species_shared)

a <- affinity(pH = c(-2, 16, 20), Eh = c(-2, 2, 20), T = 25, P = 1)
d <- diagram(a, plot.it = FALSE)

outfile <- file.path(OUTDIR, 'CHNOSZ_Pourbaix_Fe_identical_test.png')
png(outfile, width = 1600, height = 1200, res = 180)
diagram(d, lwd = 2, xlab = 'pH', ylab = 'Eh')
title(main = sprintf('CHNOSZ Fe Pourbaix — open system, log10(a(Fe2+)) = %g', LOGA_FE2),
	sub  = '(Identical shared species test vs Reaktoro)')
dev.off()

cat('Wrote:', outfile, '\n')
