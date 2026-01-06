## library() calls go here
library(conflicted)
library(dotenv)
library(targets)
library(tarchetypes)
paint::mask_print()
ggplot2::theme_set(cowplot::theme_cowplot())
library(furrr)
options(
	parallelly.fork.enable = TRUE,
	future.globals.maxSize = 2000 * 1024^2,
	tflow.report_dir = "docs"
)
plan(multicore(workers = 4))
library(ggplot2)
theme_set(cowplot::theme_cowplot())
library(ggforce)
library(quarto)
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ComplexHeatmap))
library(biomaRt)
library(categoryCompare2)
library(patchwork)
library(scico)
library(metabolomicsUtilities)
library(glue)
#tar_option_set(error = "continue")
