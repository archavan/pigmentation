# Arun Chavan
# Started: 2026-11-10

###############################################################################
# setup =======================================================================
###############################################################################

library(tidyverse)
library(here)

basedir <- here("geo/scla")
fqdir <- path.expand("~/palmer_scratch/data/bulk-rna/PIG022_scla_gut/01_fastq")

###############################################################################
# get file info ===============================================================
###############################################################################
extract_sample_name <- function(x) {
  gsub(".+\\/(.+)_S\\d+_L\\d+_R[12]_001.fastq.gz", "\\1", x)
}

extract_prefix <- function(x) {
  gsub("(.+)_R[12]_001.fastq.gz", "\\1", x)
}

extract_mate <- function(x) {
  gsub(".+(R[12])_001.fastq.gz", "\\1", x)
}

fq <- tibble(
  orig = list.files(fqdir, pattern = "*.fastq.gz", recursive = TRUE)
) |>
  mutate(seqfile_id = extract_sample_name(orig)) |>
  relocate(seqfile_id) |>
  mutate(prefix = extract_prefix(orig)) |>
  mutate(mate = extract_mate(orig))

###############################################################################
# get MD5 checksums and file sizes ============================================
###############################################################################

fq <- fq |>
  mutate(filesize = file.size(here(fqdir, orig))) |>
  mutate(md5 = tools::md5sum(here(fqdir, orig)))

###############################################################################
# count reads =================================================================
###############################################################################

count_reads <- function(.file) {
  out <- suppressWarnings(system2(
    path.expand("~/applications/bioawk/bioawk"),
    args = c("-c", "fastx", "'END{print NR}'", shQuote(.file)),
    stdout = TRUE,
    stderr = FALSE
  ))
  status <- attr(out, "status")
  if (!is.null(status) || length(out) != 1) {
    return(NA_integer_)
  }
  as.integer(out)
}

fq <- fq |>
  mutate(
    nreads = map_int(orig, ~ count_reads(here(fqdir, .x)), .progress = TRUE)
  )

###############################################################################
# write file ==================================================================
###############################################################################
write_tsv(fq, here(basedir, "results/files-info.tsv"))

# end =========================================================================
