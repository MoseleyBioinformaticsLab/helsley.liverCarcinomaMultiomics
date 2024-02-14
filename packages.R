## library() calls go here
library(conflicted)
library(dotenv)
library(targets)
library(tarchetypes)
paint::mask_print()
ggplot2::theme_set(cowplot::theme_cowplot())
library(furrr)
options(parallelly.fork.enable = TRUE,
				tflow.report_dir = "docs")