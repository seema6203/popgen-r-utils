# ---------------------------------------------------------------------------
# End-to-end SSR diversity analysis with popgenutils
#
# Replaces the usual pile of one-off scripts: read the data once, then run
# every analysis off the same object. Nothing here hardcodes a path outside
# this file, so it is safe to run from any machine.
# ---------------------------------------------------------------------------

library(popgenutils)
library(adegenet)

# --- 1. Settings -----------------------------------------------------------

input_file <- "data/genotypes.csv"   # GenAlEx-formatted CSV
output_dir <- "results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. Read the data ------------------------------------------------------

gen <- read_genalex_csv(input_file)
print(gen)

groups <- pop(gen)   # used for colours throughout

# --- 3. Diversity summary --------------------------------------------------

div <- diversity_summary(gen)
print(div)
write.csv(div, file.path(output_dir, "diversity_summary.csv"),
          row.names = FALSE)

# --- 4. Population differentiation -----------------------------------------

fst_wc  <- pairwise_fst(gen, method = "WC")
fst_nei <- pairwise_fst(gen, method = "nei")

write.csv(fst_wc,  file.path(output_dir, "pairwise_fst_weir_cockerham.csv"))
write.csv(fst_nei, file.path(output_dir, "pairwise_fst_nei.csv"))

# --- 5. AMOVA --------------------------------------------------------------

amova_res <- run_amova(gen, ~Pop, nrepet = 999, seed = 1999)
print(amova_res)

capture.output(print(amova_res),
               file = file.path(output_dir, "amova.txt"))

# --- 6. Neighbour-joining tree ---------------------------------------------

# Tip: use bootstrap = 0 while you are still adjusting the figure, then raise
# it to 1000 for the final run. Bootstrapping is by far the slow step.
nj_fit <- nj_bootstrap(gen, method = "nei", bootstrap = 1000,
                       support_cutoff = 70)
print(nj_fit)

save_figure(
  file.path(output_dir, "nj_tree.pdf"),
  plot_nj(nj_fit, groups = groups, type = "unrooted"),
  width = 10, height = 10
)

# --- 7. Principal coordinates analysis -------------------------------------

pc <- run_pcoa(gen, method = "nei", nf = 3, correction = "cailliez")
print(pc)

# Axis labels carry the real variance explained -- no hardcoded percentages.
pcoa_plot <- plot_pcoa(pc, groups = groups, axes = c(1, 2), ellipse = TRUE)

save_figure(file.path(output_dir, "pcoa.pdf"), pcoa_plot,
            width = 8, height = 7)

write.csv(pc$coords, file.path(output_dir, "pcoa_coordinates.csv"),
          row.names = FALSE)

message("Done. Output written to: ", normalizePath(output_dir))
