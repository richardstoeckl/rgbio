.rgbio_validate_location <- function(loc) {
  if (!is.character(loc) || length(loc) != 1 || is.na(loc) || !nzchar(loc)) {
    stop("Each feature 'location' must be a non-empty character string.", call. = FALSE)
  }
  loc_trim <- gsub("\\s+", "", loc)
  loc_no_keywords <- gsub("(?i)(complement|join|order|one-of)", "", loc_trim, perl = TRUE)
  loc_no_allowed <- gsub("[0-9<>,.()^]", "", loc_no_keywords)
  if (grepl("[A-Za-z]", loc_no_allowed)) {
    stop("Location string contains unsupported tokens.", call. = FALSE)
  }
}

.rgbio_ensure_qualifiers <- function(qualifiers) {
  if (!is.list(qualifiers)) {
    stop("Features 'qualifiers' column must be a list.", call. = FALSE)
  }

  for (i in seq_along(qualifiers)) {
    q <- qualifiers[[i]]
    if (!is.character(q)) {
      stop("Each element of 'qualifiers' must be a named character vector.", call. = FALSE)
    }
    q_names <- names(q)
    if (is.null(q_names) || any(is.na(q_names)) || any(q_names == "")) {
      stop("Each element of 'qualifiers' must have non-empty names.", call. = FALSE)
    }
  }
}

.rgbio_parse_feature_location <- function(location) {
  if (is.na(location) || !nzchar(location)) {
    return(list(start = NA_integer_, end = NA_integer_, strand = "*"))
  }

  strand <- if (grepl("complement", location, fixed = TRUE)) "-" else "+"
  values <- as.integer(unlist(regmatches(location, gregexpr("[0-9]+", location))))
  values <- values[!is.na(values)]

  if (length(values) == 0) {
    return(list(start = NA_integer_, end = NA_integer_, strand = strand))
  }

  list(start = min(values), end = max(values), strand = strand)
}

.rgbio_rust_records <- function(path) {
  res <- read_genbank_rs(path)
  if (is.character(res) && length(res) == 1) {
    stop(res, call. = FALSE)
  }
  if (!is.list(res) || length(res) == 0) {
    stop("No records found in GenBank file.", call. = FALSE)
  }

  lapply(res, function(rec) {
    feat <- rec$features
    features_df <- data.frame(
      key = feat$key,
      location = feat$location,
      qualifiers = I(feat$qualifiers),
      stringsAsFactors = FALSE
    )

    sequence <- rec$sequence
    if (is.character(sequence) && length(sequence) == 1 && !is.na(sequence)) {
      sequence <- toupper(sequence)
    }

    list(
      metadata = rec$metadata,
      features = features_df,
      sequence = sequence
    )
  })
}

.rgbio_write_record <- function(path, sequence, features, metadata, append) {
  .Call(wrap__write_genbank_rs, path, sequence, features, metadata, append)
}

.rgbio_validate_record <- function(record) {
  if (!is.list(record)) {
    stop("Parsed record is invalid.", call. = FALSE)
  }
  if (!is.list(record$metadata)) {
    stop("Record metadata is invalid.", call. = FALSE)
  }
  if (!is.data.frame(record$features)) {
    stop("Record features are invalid.", call. = FALSE)
  }
  if (!all(c("key", "location", "qualifiers") %in% names(record$features))) {
    stop("Record features must contain columns key, location, qualifiers.", call. = FALSE)
  }
  if (!is.character(record$sequence) || length(record$sequence) != 1) {
    stop("Record sequence is invalid.", call. = FALSE)
  }
}

.rgbio_select_records <- function(records_raw, selector) {
  if (is.null(selector)) {
    return(records_raw)
  }

  n <- length(records_raw)
  idx <- integer(0)

  if (is.numeric(selector)) {
    if (length(selector) == 0) {
      return(list())
    }
    if (any(is.na(selector)) || any(selector <= 0) || any(selector != as.integer(selector))) {
      stop("'records' numeric selectors must be positive integers.", call. = FALSE)
    }
    idx <- unique(as.integer(selector))
    if (any(idx > n)) {
      stop("'records' index out of bounds.", call. = FALSE)
    }
  } else if (is.character(selector)) {
    ids <- vapply(
      records_raw,
      function(rec) {
        acc <- rec$metadata$accession
        if (is.character(acc) && length(acc) == 1 && nzchar(acc)) {
          return(acc)
        }
        name <- rec$metadata$name
        if (is.character(name) && length(name) == 1 && nzchar(name)) {
          return(name)
        }
        ""
      },
      character(1)
    )
    idx <- which(ids %in% selector)
    if (length(idx) == 0 && length(selector) > 0) {
      stop("No records matched supplied accession/name selectors.", call. = FALSE)
    }
  } else {
    stop("'records' must be NULL, integer vector, or character vector.", call. = FALSE)
  }

  records_raw[idx]
}

.rgbio_record_ids <- function(records_raw) {
  out <- character(length(records_raw))
  for (i in seq_along(records_raw)) {
    accession <- records_raw[[i]]$metadata$accession
    name <- records_raw[[i]]$metadata$name
    out[i] <- if (is.character(accession) && length(accession) == 1 && nzchar(accession)) {
      accession
    } else if (is.character(name) && length(name) == 1 && nzchar(name)) {
      name
    } else {
      paste0("record_", i)
    }
  }
  out
}

.rgbio_build_tidy <- function(records_raw, include_sequences, include_features, include_metadata) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for tidy format.", call. = FALSE)
  }

  record_ids <- .rgbio_record_ids(records_raw)
  output <- list()

  if (include_sequences) {
    output$sequences <- tibble::tibble(
      record_id = record_ids,
      sequence = vapply(records_raw, function(rec) rec$sequence, character(1)),
      length = nchar(vapply(records_raw, function(rec) rec$sequence, character(1))),
      topology = vapply(
        records_raw,
        function(rec) {
          top <- rec$metadata$topology
          if (is.character(top) && length(top) == 1 && nzchar(top)) top else NA_character_
        },
        character(1)
      )
    )
  }

  if (include_features) {
    all_rows <- vector("list", length(records_raw))
    for (i in seq_along(records_raw)) {
      feat <- records_raw[[i]]$features
      if (nrow(feat) == 0) {
        all_rows[[i]] <- NULL
        next
      }
      parsed <- lapply(feat$location, .rgbio_parse_feature_location)
      all_rows[[i]] <- tibble::tibble(
        record_id = rep(record_ids[i], nrow(feat)),
        type = feat$key,
        start = vapply(parsed, function(x) x$start, integer(1)),
        end = vapply(parsed, function(x) x$end, integer(1)),
        strand = vapply(parsed, function(x) x$strand, character(1)),
        qualifiers = feat$qualifiers
      )
    }
    output$features <- if (length(all_rows) == 0 || all(vapply(all_rows, is.null, logical(1)))) {
      tibble::tibble(
        record_id = character(),
        type = character(),
        start = integer(),
        end = integer(),
        strand = character(),
        qualifiers = list()
      )
    } else {
      do.call(rbind, all_rows[!vapply(all_rows, is.null, logical(1))])
    }
  }

  if (include_metadata) {
    output$metadata <- tibble::tibble(
      record_id = record_ids,
      name = vapply(records_raw, function(rec) rec$metadata$name %||% "", character(1)),
      definition = vapply(records_raw, function(rec) rec$metadata$definition %||% "", character(1)),
      accession = vapply(records_raw, function(rec) rec$metadata$accession %||% "", character(1)),
      version = vapply(records_raw, function(rec) rec$metadata$version %||% "", character(1)),
      keywords = lapply(records_raw, function(rec) rec$metadata$keywords %||% character()),
      source = vapply(records_raw, function(rec) rec$metadata$source %||% "", character(1)),
      organism = vapply(records_raw, function(rec) rec$metadata$organism %||% "", character(1)),
      molecule_type = vapply(records_raw, function(rec) rec$metadata$molecule_type %||% "", character(1)),
      topology = vapply(records_raw, function(rec) rec$metadata$topology %||% "", character(1)),
      division = vapply(records_raw, function(rec) rec$metadata$division %||% "", character(1)),
      date = as.character(vapply(records_raw, function(rec) rec$metadata$date %||% NA_character_, character(1))),
      references = lapply(records_raw, function(rec) rec$metadata$references %||% list())
    )
  }

  if (length(output) == 1) {
    output[[1]]
  } else {
    output
  }
}

.rgbio_build_bioconductor <- function(records_raw, include_sequences, include_features, include_metadata) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required for bioconductor format.", call. = FALSE)
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required for bioconductor format.", call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required for bioconductor format.", call. = FALSE)
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("Package 'IRanges' is required for bioconductor format.", call. = FALSE)
  }

  record_ids <- .rgbio_record_ids(records_raw)
  output <- list()

  if (include_sequences) {
    seq_vec <- vapply(records_raw, function(rec) rec$sequence, character(1))
    dna <- Biostrings::DNAStringSet(seq_vec)
    names(dna) <- record_ids
    output$sequences <- dna
  }

  if (include_features) {
    seqnames <- character(0)
    starts <- integer(0)
    ends <- integer(0)
    strands <- character(0)
    types <- character(0)
    qualifiers <- list()
    locations <- character(0)

    for (i in seq_along(records_raw)) {
      feat <- records_raw[[i]]$features
      if (nrow(feat) == 0) {
        next
      }
      parsed <- lapply(feat$location, .rgbio_parse_feature_location)
      seqnames <- c(seqnames, rep(record_ids[i], nrow(feat)))
      starts <- c(starts, vapply(parsed, function(x) ifelse(is.na(x$start), 1L, x$start), integer(1)))
      ends <- c(ends, vapply(parsed, function(x) ifelse(is.na(x$end), 1L, x$end), integer(1)))
      strands <- c(strands, vapply(parsed, function(x) x$strand, character(1)))
      types <- c(types, feat$key)
      qualifiers <- c(qualifiers, feat$qualifiers)
      locations <- c(locations, feat$location)
    }

    ranges <- IRanges::IRanges(start = starts, end = ends)
    gr <- GenomicRanges::GRanges(seqnames = seqnames, ranges = ranges, strand = strands)
    S4Vectors::mcols(gr)$type <- types
    S4Vectors::mcols(gr)$qualifiers <- qualifiers
    S4Vectors::mcols(gr)$location <- locations
    output$features <- gr
  }

  if (include_metadata) {
    output$metadata <- S4Vectors::DataFrame(
      record_id = record_ids,
      name = vapply(records_raw, function(rec) rec$metadata$name %||% "", character(1)),
      definition = vapply(records_raw, function(rec) rec$metadata$definition %||% "", character(1)),
      accession = vapply(records_raw, function(rec) rec$metadata$accession %||% "", character(1)),
      version = vapply(records_raw, function(rec) rec$metadata$version %||% "", character(1)),
      keywords = I(lapply(records_raw, function(rec) rec$metadata$keywords %||% character())),
      source = vapply(records_raw, function(rec) rec$metadata$source %||% "", character(1)),
      organism = vapply(records_raw, function(rec) rec$metadata$organism %||% "", character(1)),
      molecule_type = vapply(records_raw, function(rec) rec$metadata$molecule_type %||% "", character(1)),
      topology = vapply(records_raw, function(rec) rec$metadata$topology %||% "", character(1)),
      division = vapply(records_raw, function(rec) rec$metadata$division %||% "", character(1)),
      date = vapply(records_raw, function(rec) as.character(rec$metadata$date %||% NA_character_), character(1)),
      references = I(lapply(records_raw, function(rec) rec$metadata$references %||% list()))
    )
  }

  if (length(output) == 1) {
    output[[1]]
  } else {
    output
  }
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Read a GenBank file
#'
#' Reads one or more GenBank records and returns selected components in either
#' Bioconductor or tidy format.
#'
#' @param file Character path to a GenBank file.
#' @param format Output format. One of "bioconductor" or "tidy".
#' @param sequences Logical; include sequence data.
#' @param features Logical; include feature annotations.
#' @param metadata Logical; include record metadata.
#' @param records Integer indices or character accession/name selectors.
#' @param validate Logical; validate parsed records.
#' @return A variable object based on selected components. Returns either a
#' single object or a named list with `sequences`, `features`, and/or `metadata`.
#' @export
read_gbk <- function(
  file,
  format = "bioconductor",
  sequences = TRUE,
  features = TRUE,
  metadata = TRUE,
  records = NULL,
  validate = TRUE
) {
  if (!is.character(file) || length(file) != 1 || is.na(file) || !nzchar(file)) {
    stop("'file' must be a non-empty character path.", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }
  if (!is.character(format) || length(format) != 1 || !format %in% c("bioconductor", "tidy")) {
    stop("'format' must be one of: 'bioconductor', 'tidy'.", call. = FALSE)
  }
  if (!is.logical(sequences) || length(sequences) != 1 || is.na(sequences)) {
    stop("'sequences' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(features) || length(features) != 1 || is.na(features)) {
    stop("'features' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(metadata) || length(metadata) != 1 || is.na(metadata)) {
    stop("'metadata' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(records) && !is.numeric(records) && !is.character(records)) {
    stop("'records' must be NULL, integer vector, or character vector.", call. = FALSE)
  }
  if (!is.logical(validate) || length(validate) != 1 || is.na(validate)) {
    stop("'validate' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!any(c(sequences, features, metadata))) {
    stop("At least one of 'sequences', 'features', or 'metadata' must be TRUE.", call. = FALSE)
  }

  path <- normalizePath(file, mustWork = TRUE)
  records_raw <- .rgbio_rust_records(path)
  records_raw <- .rgbio_select_records(records_raw, records)

  if (validate) {
    for (record in records_raw) {
      .rgbio_validate_record(record)
    }
  }

  if (format == "tidy") {
    return(.rgbio_build_tidy(records_raw, sequences, features, metadata))
  }

  .rgbio_build_bioconductor(records_raw, sequences, features, metadata)
}

.rgbio_normalize_sequences <- function(sequences) {
  if (inherits(sequences, "DNAStringSet")) {
    seq_vec <- as.character(sequences)
    seq_names <- names(sequences)
  } else if (is.character(sequences)) {
    seq_vec <- sequences
    seq_names <- names(sequences)
  } else {
    stop("'sequences' must be a DNAStringSet or named character vector.", call. = FALSE)
  }

  if (length(seq_vec) == 0) {
    stop("'sequences' cannot be empty.", call. = FALSE)
  }
  if (anyNA(seq_vec) || any(!nzchar(seq_vec))) {
    stop("All sequences must be non-empty character strings.", call. = FALSE)
  }

  if (is.null(seq_names) || any(!nzchar(seq_names))) {
    seq_names <- paste0("record_", seq_along(seq_vec))
  }

  list(values = toupper(seq_vec), names = seq_names)
}

.rgbio_normalize_features <- function(features, sequence_name) {
  if (is.null(features)) {
    return(data.frame(
      key = character(),
      location = character(),
      qualifiers = I(list()),
      stringsAsFactors = FALSE
    ))
  }

  if (inherits(features, "GRanges")) {
    mcols <- S4Vectors::mcols(features)
    if (!all(c("type", "qualifiers") %in% names(mcols))) {
      stop("GRanges features must contain mcols 'type' and 'qualifiers'.", call. = FALSE)
    }

    starts <- as.integer(GenomicRanges::start(features))
    ends <- as.integer(GenomicRanges::end(features))
    strands <- as.character(GenomicRanges::strand(features))
    locations <- ifelse(strands == "-", paste0("complement(", starts, "..", ends, ")"), paste0(starts, "..", ends))

    return(data.frame(
      key = as.character(mcols$type),
      location = locations,
      qualifiers = I(as.list(mcols$qualifiers)),
      stringsAsFactors = FALSE
    ))
  }

  if (!is.data.frame(features)) {
    stop("'features' must be NULL, GRanges, tibble, or data.frame.", call. = FALSE)
  }
  if (!all(c("type", "start", "end", "strand", "qualifiers") %in% names(features))) {
    stop("Feature table must contain columns: type, start, end, strand, qualifiers.", call. = FALSE)
  }

  location <- ifelse(
    as.character(features$strand) == "-",
    paste0("complement(", as.integer(features$start), "..", as.integer(features$end), ")"),
    paste0(as.integer(features$start), "..", as.integer(features$end)
  ))

  data.frame(
    key = as.character(features$type),
    location = location,
    qualifiers = I(as.list(features$qualifiers)),
    stringsAsFactors = FALSE
  )
}

.rgbio_normalize_metadata <- function(metadata, sequence_name) {
  if (is.null(metadata)) {
    return(list(name = sequence_name, definition = sequence_name, accession = sequence_name))
  }

  if (inherits(metadata, "DataFrame") || is.data.frame(metadata)) {
    if (nrow(metadata) == 0) {
      return(list(name = sequence_name, definition = sequence_name, accession = sequence_name))
    }
    metadata <- as.list(metadata[1, , drop = FALSE])
  }

  if (!is.list(metadata)) {
    stop("'metadata' must be NULL, DataFrame, data.frame, or list.", call. = FALSE)
  }

  if (is.null(metadata$name) || !nzchar(metadata$name)) {
    metadata$name <- sequence_name
  }
  if (is.null(metadata$definition) || !nzchar(metadata$definition)) {
    metadata$definition <- sequence_name
  }
  if (is.null(metadata$accession) || !nzchar(metadata$accession)) {
    metadata$accession <- sequence_name
  }
  if (is.null(metadata$molecule_type) || !nzchar(metadata$molecule_type)) {
    metadata$molecule_type <- "DNA"
  }

  metadata
}

.rgbio_features_for_record <- function(features, record_id, record_index, n_records) {
  if (is.null(features)) {
    return(NULL)
  }

  if (inherits(features, "GRanges")) {
    seqnames_vec <- as.character(GenomicRanges::seqnames(features))
    if (length(seqnames_vec) == 0) {
      return(features)
    }
    selected <- which(seqnames_vec == record_id)
    if (length(selected) == 0 && n_records > 1) {
      return(features[0])
    }
    if (length(selected) == 0) {
      return(features)
    }
    return(features[selected])
  }

  if (!is.data.frame(features)) {
    stop("'features' must be NULL, GRanges, tibble, or data.frame.", call. = FALSE)
  }

  if ("record_id" %in% names(features)) {
    subset_rows <- which(as.character(features$record_id) == record_id)
    return(features[subset_rows, setdiff(names(features), "record_id"), drop = FALSE])
  }

  if (n_records > 1) {
    return(features[0, , drop = FALSE])
  }

  features
}

.rgbio_metadata_for_record <- function(metadata, record_id, record_index, n_records) {
  if (is.null(metadata)) {
    return(NULL)
  }

  if (is.list(metadata) && !inherits(metadata, "data.frame") && !inherits(metadata, "DataFrame")) {
    if (!is.null(names(metadata)) && record_id %in% names(metadata) && is.list(metadata[[record_id]])) {
      return(metadata[[record_id]])
    }
    if (length(metadata) == n_records && all(vapply(metadata, is.list, logical(1)))) {
      return(metadata[[record_index]])
    }
    return(metadata)
  }

  if (inherits(metadata, "DataFrame") || is.data.frame(metadata)) {
    if (nrow(metadata) == 0) {
      return(metadata)
    }

    if ("record_id" %in% names(metadata)) {
      selected <- which(as.character(metadata$record_id) == record_id)
      if (length(selected) == 0) {
        if (n_records > 1) {
          stop("No metadata row matched record_id: ", record_id, call. = FALSE)
        }
        selected <- 1L
      }
      return(metadata[selected[1], setdiff(names(metadata), "record_id"), drop = FALSE])
    }

    if (nrow(metadata) == n_records) {
      return(metadata[record_index, , drop = FALSE])
    }

    return(metadata[1, , drop = FALSE])
  }

  metadata
}

#' Write GenBank records
#'
#' Writes GenBank records from sequence, feature, and metadata components.
#'
#' @param file Character path to output file.
#' @param sequences DNAStringSet or named character vector.
#' @param features GRanges with `type` and `qualifiers` in `mcols()`, or tidy
#' feature table with columns `type`, `start`, `end`, `strand`, `qualifiers`.
#' @param metadata DataFrame, data.frame, or list with record metadata.
#' @param append Logical; append to file.
#' @param validate Logical; validate inputs.
#' @param line_width Integer sequence line width.
#' @return Logical TRUE on success.
#' @export
write_gbk <- function(
  file,
  sequences,
  features = NULL,
  metadata = NULL,
  append = FALSE,
  validate = TRUE,
  line_width = 80
) {
  if (!is.character(file) || length(file) != 1 || is.na(file) || !nzchar(file)) {
    stop("'file' must be a non-empty character path.", call. = FALSE)
  }
  if (!is.logical(append) || length(append) != 1 || is.na(append)) {
    stop("'append' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(validate) || length(validate) != 1 || is.na(validate)) {
    stop("'validate' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(line_width) || length(line_width) != 1 || is.na(line_width) || line_width < 1) {
    stop("'line_width' must be a positive integer.", call. = FALSE)
  }

  if (append && !file.exists(file)) {
    stop("append = TRUE requires an existing GenBank file.", call. = FALSE)
  }
  if (append) {
    existing_path <- normalizePath(file, mustWork = TRUE)
    existing <- read_genbank_rs(existing_path)
    if (is.character(existing) && length(existing) == 1) {
      stop("Cannot append: existing file is not a valid GenBank file.", call. = FALSE)
    }
  }
  if (line_width != 80) {
    warning("'line_width' is currently ignored by the Rust backend; using backend default.", call. = FALSE)
  }

  seq_norm <- .rgbio_normalize_sequences(sequences)

  path <- normalizePath(file, mustWork = FALSE)
  n_records <- length(seq_norm$values)

  for (i in seq_len(n_records)) {
    record_id <- seq_norm$names[[i]]

    feature_i <- .rgbio_features_for_record(features, record_id, i, n_records)
    metadata_i <- .rgbio_metadata_for_record(metadata, record_id, i, n_records)
    feature_norm <- .rgbio_normalize_features(feature_i, record_id)
    metadata_norm <- .rgbio_normalize_metadata(metadata_i, record_id)

    if (validate) {
      if (!all(c("key", "location", "qualifiers") %in% names(feature_norm))) {
        stop("Normalized feature table is invalid.", call. = FALSE)
      }
      if (nrow(feature_norm) > 0) {
        .rgbio_ensure_qualifiers(feature_norm$qualifiers)
        for (loc in feature_norm$location) {
          .rgbio_validate_location(loc)
        }
      }
    }

    append_now <- append || i > 1
    res <- .rgbio_write_record(path, seq_norm$values[[i]], feature_norm, metadata_norm, append_now)
    if (is.character(res) && length(res) == 1) {
      stop(res, call. = FALSE)
    }
  }

  TRUE
}

#' Read a GenBank file
#'
#' Deprecated compatibility wrapper for `read_gbk()`.
#'
#' @param file Path to a GenBank file.
#' @return Parsed records in legacy list format.
#' @export
read_genbank <- function(file) {
  .Deprecated("read_gbk")
  if (!is.character(file) || length(file) != 1 || is.na(file) || !nzchar(file)) {
    stop("'file' must be a non-empty character path.", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }
  .rgbio_rust_records(normalizePath(file, mustWork = TRUE))
}

#' Write a GenBank file
#'
#' Deprecated compatibility wrapper for `write_gbk()`.
#'
#' @param file Path to output file.
#' @param sequence Sequence string.
#' @param features Feature table with columns key, location, qualifiers.
#' @param metadata Metadata list.
#' @return Logical TRUE on success.
#' @export
write_genbank <- function(file, sequence, features, metadata = list()) {
  .Deprecated("write_gbk")
  if (!is.character(sequence) || length(sequence) != 1 || is.na(sequence) || !nzchar(sequence)) {
    stop("'sequence' must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.data.frame(features) || !all(c("key", "location", "qualifiers") %in% names(features))) {
    stop("Legacy 'features' must be a data.frame with key, location, qualifiers columns.", call. = FALSE)
  }
  if (!is.list(metadata)) {
    stop("'metadata' must be a list.", call. = FALSE)
  }
  .rgbio_ensure_qualifiers(features$qualifiers)
  for (loc in features$location) {
    .rgbio_validate_location(loc)
  }

  path <- normalizePath(file, mustWork = FALSE)
  res <- .rgbio_write_record(path, toupper(sequence), features, metadata, FALSE)
  if (is.character(res) && length(res) == 1) {
    stop(res, call. = FALSE)
  }
  res
}
