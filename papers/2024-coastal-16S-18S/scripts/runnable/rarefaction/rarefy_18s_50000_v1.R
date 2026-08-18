#!/usr/bin/env Rscript

# Isolated 18S rarefaction worker.
#
# This script deliberately performs only count-table validation and
# vegan::rrarefy.  Linked taxonomy/FASTA packaging and package provenance are
# handled by the guarded Python orchestrator with the same version suffix.

WORKER_VERSION <- "1.0.0"
EXPECTED_FEATURES <- 15601L
EXPECTED_SAMPLES <- 75L
DEFAULT_DEPTH <- 50000L
DEFAULT_SEED <- 1L
RNG_CONTRACT <- list(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

usage <- function(status = 0L) {
  message(paste(
    "Usage:",
    "  Rscript rarefy_18s_50000_v1.R --dry-run --input COUNTS.tsv [--depth 50000 --seed 1]",
    "  Rscript rarefy_18s_50000_v1.R --run --input COUNTS.tsv --output RAREFIED.tsv [--depth 50000 --seed 1]",
    "  Rscript rarefy_18s_50000_v1.R --check --input COUNTS.tsv --expected RAREFIED.tsv [--depth 50000 --seed 1]",
    "",
    "The input and output orientation is ASV x sample. Internally,",
    "vegan::rrarefy receives sample x ASV. Sampling is without replacement.",
    sep = "\n"
  ))
  quit(save = "no", status = status, runLast = FALSE)
}

fail <- function(message) {
  stop(message, call. = FALSE)
}

require_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    fail(message)
  }
}

parse_integer <- function(value, label) {
  parsed <- suppressWarnings(as.integer(value))
  require_true(
    length(parsed) == 1L && !is.na(parsed) && as.character(parsed) == value,
    paste0(label, " must be an integer; observed: ", value)
  )
  parsed
}

parse_args <- function(args) {
  if (length(args) == 0L) {
    usage(2L)
  }
  if ("--help" %in% args || "-h" %in% args) {
    usage(0L)
  }

  result <- list(
    run = FALSE,
    dry_run = FALSE,
    check = FALSE,
    input = NULL,
    output = NULL,
    expected = NULL,
    depth = DEFAULT_DEPTH,
    seed = DEFAULT_SEED
  )
  flag_names <- c("--run", "--dry-run", "--check")
  value_names <- c("--input", "--output", "--expected", "--depth", "--seed")

  index <- 1L
  while (index <= length(args)) {
    token <- args[[index]]
    if (token %in% flag_names) {
      if (token == "--run") result$run <- TRUE
      if (token == "--dry-run") result$dry_run <- TRUE
      if (token == "--check") result$check <- TRUE
      index <- index + 1L
      next
    }
    require_true(token %in% value_names, paste0("Unknown argument: ", token))
    require_true(index < length(args), paste0("Missing value for ", token))
    value <- args[[index + 1L]]
    require_true(!startsWith(value, "--"), paste0("Missing value for ", token))
    if (token == "--input") result$input <- value
    if (token == "--output") result$output <- value
    if (token == "--expected") result$expected <- value
    if (token == "--depth") result$depth <- parse_integer(value, "--depth")
    if (token == "--seed") result$seed <- parse_integer(value, "--seed")
    index <- index + 2L
  }

  mode_count <- sum(c(result$run, result$dry_run, result$check))
  require_true(mode_count == 1L, "Select exactly one of --run, --dry-run, or --check")
  require_true(!is.null(result$input), "--input is required")
  require_true(result$depth == DEFAULT_DEPTH, "This version is frozen at depth 50000")
  require_true(result$seed == DEFAULT_SEED, "This version is frozen at seed 1")
  if (result$run) require_true(!is.null(result$output), "--output is required with --run")
  if (result$check) require_true(!is.null(result$expected), "--expected is required with --check")
  result
}

read_count_table <- function(path) {
  require_true(file.exists(path), paste0("Missing count table: ", path))
  require_true(!dir.exists(path), paste0("Count-table path is a directory: ", path))

  table <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character",
    na.strings = character(0)
  )
  require_true(ncol(table) >= 2L, paste0("Malformed count table: ", path))
  require_true(nrow(table) == EXPECTED_FEATURES, paste0(
    "Expected ", EXPECTED_FEATURES, " ASVs; observed ", nrow(table), " in ", path
  ))
  require_true(ncol(table) - 1L == EXPECTED_SAMPLES, paste0(
    "Expected ", EXPECTED_SAMPLES, " samples; observed ", ncol(table) - 1L, " in ", path
  ))

  id_header <- names(table)[[1L]]
  feature_ids <- table[[1L]]
  sample_ids <- names(table)[-1L]
  require_true(nzchar(id_header), "Empty feature-ID header")
  require_true(all(nzchar(feature_ids)), "Empty feature ID")
  require_true(all(nzchar(sample_ids)), "Empty sample ID")
  require_true(!anyDuplicated(feature_ids), "Duplicate feature IDs")
  require_true(!anyDuplicated(sample_ids), "Duplicate sample IDs")
  require_true(all(startsWith(feature_ids, "18S_ASV_")), "Unexpected non-18S feature ID")
  require_true(
    !any(grepl("^SS_BS(?:_|$)", sample_ids, perl = TRUE)),
    "SS_BS is forbidden from this 75-sample package"
  )

  values_character <- unlist(table[-1L], use.names = FALSE)
  values_numeric <- suppressWarnings(as.numeric(values_character))
  require_true(length(values_numeric) == length(values_character), "Count conversion changed cell count")
  require_true(!anyNA(values_numeric), "Non-numeric or missing count cell")
  require_true(all(is.finite(values_numeric)), "Non-finite count cell")
  require_true(all(values_numeric >= 0), "Negative count cell")
  require_true(all(values_numeric == floor(values_numeric)), "Non-integer count cell")
  require_true(all(values_numeric <= .Machine$integer.max), "Count exceeds R integer range")

  counts <- matrix(
    as.integer(values_numeric),
    nrow = nrow(table),
    ncol = ncol(table) - 1L,
    dimnames = list(feature_ids, sample_ids)
  )
  list(id_header = id_header, counts = counts)
}

validate_depth <- function(counts, depth) {
  sample_depths <- colSums(counts)
  require_true(all(sample_depths >= depth), paste0(
    "At least one sample is below ", depth, " reads: ",
    paste(names(sample_depths)[sample_depths < depth], collapse = ", ")
  ))
  sample_depths
}

configure_rng <- function(seed) {
  RNGkind(
    kind = RNG_CONTRACT$kind,
    normal.kind = RNG_CONTRACT$normal.kind,
    sample.kind = RNG_CONTRACT$sample.kind
  )
  set.seed(seed)
}

rarefy_counts <- function(counts, depth, seed) {
  configure_rng(seed)
  sample_by_feature <- t(counts)
  rare_sample_by_feature <- vegan::rrarefy(sample_by_feature, sample = depth)
  rarefied <- t(rare_sample_by_feature)
  dimnames(rarefied) <- dimnames(counts)
  storage.mode(rarefied) <- "integer"

  require_true(identical(dim(rarefied), dim(counts)), "Rarefaction changed table dimensions")
  require_true(identical(rownames(rarefied), rownames(counts)), "Rarefaction changed feature order")
  require_true(identical(colnames(rarefied), colnames(counts)), "Rarefaction changed sample order")
  require_true(all(rarefied >= 0L), "Rarefaction generated negative counts")
  require_true(all(rarefied <= counts), "Rarefaction generated a cell larger than its source")
  require_true(all(colSums(rarefied) == depth), "Rarefied sample total is not exactly the selected depth")
  rarefied
}

write_count_table <- function(path, id_header, counts) {
  require_true(!file.exists(path), paste0("Refusing to overwrite existing output: ", path))
  parent <- dirname(path)
  require_true(dir.exists(parent), paste0("Output parent does not exist: ", parent))

  output <- data.frame(
    feature_id = rownames(counts),
    as.data.frame(counts, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(output)[[1L]] <- id_header
  write.table(
    output,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "",
    eol = "\n"
  )
}

print_contract <- function(mode, data, sample_depths, depth, seed) {
  cat(paste0("worker_version=", WORKER_VERSION, "\n"))
  cat(paste0("mode=", mode, "\n"))
  cat(paste0("features=", nrow(data$counts), "\n"))
  cat(paste0("samples=", ncol(data$counts), "\n"))
  cat(paste0("source_min_depth=", min(sample_depths), "\n"))
  cat(paste0("source_max_depth=", max(sample_depths), "\n"))
  cat(paste0("depth=", depth, "\n"))
  cat(paste0("seed=", seed, "\n"))
  cat("algorithm=vegan::rrarefy\n")
  cat("replacement=false\n")
  cat(paste0("rng_kind=", RNG_CONTRACT$kind, "\n"))
  cat(paste0("rng_normal_kind=", RNG_CONTRACT$normal.kind, "\n"))
  cat(paste0("rng_sample_kind=", RNG_CONTRACT$sample.kind, "\n"))
  cat(paste0("r_version=", R.version.string, "\n"))
  cat(paste0("vegan_version=", as.character(utils::packageVersion("vegan")), "\n"))
  cat(paste0("vegan_built=", utils::packageDescription("vegan")$Built, "\n"))
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  requireNamespace("vegan", quietly = TRUE)
  data <- read_count_table(opts$input)
  sample_depths <- validate_depth(data$counts, opts$depth)
  mode <- if (opts$run) "run" else if (opts$dry_run) "dry-run" else "check"
  print_contract(mode, data, sample_depths, opts$depth, opts$seed)

  if (opts$dry_run) {
    cat("status=dry_run_validated\n")
    return(invisible(NULL))
  }

  rarefied <- rarefy_counts(data$counts, opts$depth, opts$seed)
  if (opts$run) {
    write_count_table(opts$output, data$id_header, rarefied)
    cat(paste0("rarefied_total=", sum(rarefied), "\n"))
    cat(paste0("post_rarefaction_zero_features=", sum(rowSums(rarefied) == 0L), "\n"))
    cat("status=run_complete\n")
    return(invisible(NULL))
  }

  expected <- read_count_table(opts$expected)
  require_true(identical(data$id_header, expected$id_header), "Expected-table ID header differs")
  require_true(identical(dimnames(data$counts), dimnames(expected$counts)), "Expected-table IDs/order differ")
  mismatches <- sum(rarefied != expected$counts)
  require_true(mismatches == 0L, paste0("Deterministic recomputation mismatches: ", mismatches))
  cat("deterministic_recomputation_mismatches=0\n")
  cat("status=check_passed\n")
}

tryCatch(
  main(),
  error = function(error) {
    message(paste0("ERROR: ", conditionMessage(error)))
    quit(save = "no", status = 1L, runLast = FALSE)
  }
)
