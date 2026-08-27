#!/usr/bin/env Rscript

suppressMessages(library(CHNOSZ))

`%||%` <- function(a, b) if (is.null(a)) b else a

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- args_full[grep("^--file=", args_full)]
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else NULL
BASE <- if (!is.null(script_path) && nzchar(script_path)) {
  normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

T_C <- 350.0
P_BAR <- 2000.0
LOGFO2_MIN <- -48.0
LOGFO2_MAX <- -18.0
PH_MIN <- 3.0
PH_MAX <- 8.5
N_FO2 <- 140
N_PH <- 120
LOGA_FE2 <- 0.0

reset()
suppressMessages(T.units("C"))
suppressMessages(P.units("bar"))

basis(c("Fe+2", "H2O", "H+", "O2"))
basis("Fe+2", LOGA_FE2)

species(c(
  "Fe+2", "Fe+3", "FeOH+", "FeOH+2", "HFeO2", "FeO2-",
  "hematite", "magnetite", "goethite", "iron"
))

a <- affinity(
  O2 = c(LOGFO2_MIN, LOGFO2_MAX, N_FO2),
  pH = c(PH_MIN, PH_MAX, N_PH),
  T = T_C,
  P = P_BAR
)

out_png <- file.path(BASE, "CHNOSZ_DeepFluid_LogfO2_pH_Fe.png")
out_txt <- file.path(BASE, "deepfluid_logfo2_ph_setup_chnosz.txt")

png(out_png, width = 1700, height = 1200, res = 180)
par(mar = c(4.2, 4.2, 2.8, 1.2))
diagram(
  a,
  fill = "terrain",
  names = TRUE,
  cex.names = 0.7,
  lwd = 1.2,
  xlab = expression(log[10] * f[O[2]]),
  ylab = "pH"
)
abline(v = c(-24, -30, -36), lty = 3, col = "gray40", lwd = 1)
text(x = c(-24, -30, -36), y = rep(PH_MAX - 0.2, 3), labels = c("HM-ish", "FMQ-ish", "IW-ish"), cex = 0.7, col = "gray30", pos = 4)
title(main = sprintf("CHNOSZ deep-fluid Fe stability: log10(fO2) vs pH\nT = %.0f C, P = %.0f bar", T_C, P_BAR))
dev.off()

writeLines(c(
  "Deep-fluid potential diagram setup (CHNOSZ)",
  sprintf("- T = %.1f C", T_C),
  sprintf("- P = %.1f bar", P_BAR),
  sprintf("- log10(fO2) range: [%.1f, %.1f] bar", LOGFO2_MIN, LOGFO2_MAX),
  sprintf("- pH range: [%.1f, %.1f]", PH_MIN, PH_MAX),
  sprintf("- N(logfO2) = %d", N_FO2),
  sprintf("- N(pH) = %d", N_PH),
  sprintf("- Fixed iron basis activity: log10(a(Fe+2)) = %.1f", LOGA_FE2),
  "- Diagram type: logfO2-pH (deep-fluid style potential space)",
  sprintf("- Output: %s", normalizePath(out_png, winslash = "/", mustWork = FALSE))
), out_txt)

cat(sprintf("Wrote: %s\n", out_png))
cat(sprintf("Wrote: %s\n", out_txt))
