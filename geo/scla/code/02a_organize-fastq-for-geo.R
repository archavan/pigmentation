# Arun Chavan
# Started: 2026-11-10

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)

basedir <- here("geo/scla")
scdir <- here("species/scla/bulk/PIG022/refseq_v1")
fqdir <- path.expand("~/palmer_scratch/data/bulk-rna/PIG022_scla_gut/01_fastq")

resdir <- path.expand("~/palmer_scratch/geo/pigmentation/scla/fq")
fs::dir_create(here(resdir))

###############################################################################
# files info ==================================================================
###############################################################################

fq <- read_tsv(here(basedir, "results/files-info.tsv"))

md <- read_csv(here(scdir, "data/samples.csv")) |>
  select(sample_id, seqfile_id, animal, sample, rep) |>
  left_join(fq, by = "seqfile_id")

###############################################################################
# exclude bad and unused files ================================================
###############################################################################

## Either empty or with a very small number of reads (bad BCL convert runs)

reads_cutoff <- 5000

keep_runs <- md |>
  summarize(mates_above_cutoff = all(nreads > reads_cutoff), .by = prefix) |>
  filter(mates_above_cutoff) |>
  pull(prefix)

md <- md |>
  filter(prefix %in% keep_runs) |>
  mutate(run = as.integer(as.factor(prefix)), .by = sample_id)

## replicate 3 (animal D) were outliers excluded from all analyses.
md <- md |>
  filter(animal != "D")

## sanity checks

md |>
  summarize(nreads = list(nreads), .by = prefix) |>
  mutate(mates_equal = map_lgl(nreads, ~ identical(.x[1], .x[2]))) |>
  pull(mates_equal) |>
  all() |>
  stopifnot()

###############################################################################
# copy files with a GEO-ready name to a new dir ===============================
###############################################################################

md <- md |>
  select(-prefix) |>
  mutate(fq_new = paste0(sample_id, "_run", run, "_", mate, ".fastq.gz"))

stopifnot(!any(duplicated(md$fq_new)))

copy_safely <- function(.source, .dest) {
  md5_source <- unname(tools::md5sum(.source))
  already_ok <- fs::file_exists(.dest) &&
    identical(md5_source, unname(tools::md5sum(.dest)))

  if (!already_ok) {
    tmp <- paste0(.dest, ".part")
    fs::file_copy(.source, tmp, overwrite = TRUE)
    if (!identical(md5_source, unname(tools::md5sum(tmp)))) {
      fs::file_delete(tmp)
      stop("md5 mismatch: ", fs::path_file(.source))
    }
    fs::file_move(tmp, .dest)
    message("Copied to: ", fs::path_file(.dest))
  } else {
    message("Already okay: ", fs::path_file(.dest))
  }
  invisible(.dest)
}

walk2(
  md$orig,
  md$fq_new,
  function(x, y) {
    copy_safely(here(fqdir, x), here(resdir, y))
  },
  .progress = TRUE
)

###############################################################################
# prepare file paths to save to a file ========================================
###############################################################################

fq_new_wide <- md |>
  select(sample_id, run, mate, fq_new) |>
  pivot_wider(
    id_cols = c(sample_id, run),
    names_from = mate,
    values_from = fq_new
  )

fq_all_tidy <- md |>
  select(
    sample_id,
    run,
    mate,
    fq_orig = orig,
    fq_new,
    md5
  )

###############################################################################
# save files ==================================================================
###############################################################################

write_csv(fq_new_wide, here(basedir, "results/fq_renamed_wide.csv"))
write_csv(fq_all_tidy, here(basedir, "results/fq_all_tidy.csv"))

# end #########################################################################
