read_fasta_single <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0 || !startsWith(lines[1], ">")) {
    stop("FASTA header not found.")
  }

  header <- sub("^>", "", lines[1])
  seq_lines <- lines[-1]
  seq_lines <- seq_lines[nzchar(seq_lines)]
  sequence <- toupper(gsub("\\s+", "", paste(seq_lines, collapse = "")))

  list(header = header, sequence = sequence)
}

parse_gff_attributes <- function(attr) {
  parts <- strsplit(attr, ";", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]

  out <- character(0)
  for (part in parts) {
    kv <- strsplit(part, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) {
      next
    }
    out[kv[1]] <- kv[2]
  }
  out
}

location_from_gff <- function(start, end, strand, start_range, end_range) {
  start_str <- as.character(start)
  end_str <- as.character(end)

  if (!is.na(start_range) && grepl("^\\.\\,", start_range)) {
    start_str <- paste0("<", start_str)
  }
  if (!is.na(end_range) && grepl("\\,\\.$", end_range)) {
    end_str <- paste0(">", end_str)
  }

  loc <- paste0(start_str, "..", end_str)
  if (strand == "-") {
    loc <- paste0("complement(", loc, ")")
  }
  loc
}

gff_to_features <- function(path, organism) {
  lines <- readLines(path)
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nzchar(lines)]

  cols <- strsplit(lines, "\t", fixed = TRUE)
  cols <- cols[vapply(cols, length, integer(1)) == 9]

  keys <- character(0)
  locs <- character(0)
  qualifiers <- list()

  for (row in cols) {
    type <- row[[3]]
    if (type == "exon") {
      next
    }

    attrs <- parse_gff_attributes(row[[9]])
    key <- if (type == "region") "source" else type

    start_range <- attrs["start_range"]
    end_range <- attrs["end_range"]
    start_val <- as.integer(row[[4]])
    end_val <- as.integer(row[[5]])
    strand <- row[[7]]

    loc <- location_from_gff(start_val, end_val, strand, start_range, end_range)

    q <- character(0)
    if (key == "source") {
      q["organism"] <- organism
    }

    if ("Dbxref" %in% names(attrs)) {
      q["db_xref"] <- attrs[["Dbxref"]]
    }
    if ("chromosome" %in% names(attrs)) {
      q["chromosome"] <- attrs[["chromosome"]]
    }
    if ("mol_type" %in% names(attrs)) {
      q["mol_type"] <- attrs[["mol_type"]]
    }

    copy_map <- c(
      gene = "gene",
      product = "product",
      protein_id = "protein_id",
      Note = "note"
    )

    for (name in names(copy_map)) {
      if (name %in% names(attrs)) {
        q[copy_map[[name]]] <- attrs[[name]]
      }
    }

    if (key == "CDS") {
      phase <- row[[8]]
      if (!is.na(phase) && phase != ".") {
        q["codon_start"] <- as.character(as.integer(phase) + 1)
      }
    }

    keys <- c(keys, key)
    locs <- c(locs, loc)
    qualifiers[[length(qualifiers) + 1]] <- q
  }

  features <- data.frame(
    key = keys,
    location = locs,
    stringsAsFactors = FALSE
  )
  features$qualifiers <- qualifiers
  features
}

compare_files_strict <- function(path_a, path_b) {
  lines_a <- readLines(path_a)
  lines_b <- readLines(path_b)
  identical(lines_a, lines_b)
}
