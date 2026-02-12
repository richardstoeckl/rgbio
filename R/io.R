#' Read a GenBank file
#'
#' Parses a GenBank file and returns a list of sequences, features, and metadata.
#'
#' @param file Path to the GenBank file.
#' @return A list of records, where each record contains:
#'   \item{metadata}{A list of metadata with supported fields:
#'     \itemize{
#'       \item \code{name} (Locus name)
#'       \item \code{definition}
#'       \item \code{accession}
#'       \item \code{version}
#'       \item \code{keywords} (character vector)
#'       \item \code{source}
#'       \item \code{organism}
#'       \item \code{molecule_type} (e.g., "DNA")
#'       \item \code{division}
#'       \item \code{topology} ("linear" or "circular")
#'       \item \code{date} (format: \code{DD-MON-YYYY})
#'       \item \code{references} (list of references; each reference may include
#'         \code{description}, \code{authors}, \code{consortium}, \code{title},
#'         \code{journal}, \code{pubmed}, \code{remark})
#'     }
#'   }
#'   \item{features}{A data frame of features (key, location, qualifiers)}
#'   \item{sequence}{The sequence as a string}
#' @export
#' @examples
#' \dontrun{
#'   records <- read_genbank("example.gb")
#' }
read_genbank <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  path <- normalizePath(file, mustWork = TRUE)
  res <- read_genbank_rs(path)
  if (is.character(res) && length(res) == 1) {
    stop(res)
  }
  if (is.list(res) && length(res) == 0) {
    stop("No records found in GenBank file.")
  }
  
  # Post-process to tidy data.frame
  lapply(res, function(rec) {
      # features is a list of vectors: list(key=..., location=..., qualifiers=...)
      # We convert this to a DataFrame manually to ensure list column is preserved key
      f <- rec$features
      # Ensure lengths match (they should from Rust)
      rec$features <- data.frame(
          key = f$key,
          location = f$location,
          qualifiers = I(f$qualifiers), # Protect list
          stringsAsFactors = FALSE
      )
      if (is.character(rec$sequence) && length(rec$sequence) == 1 && !is.na(rec$sequence)) {
        rec$sequence <- toupper(rec$sequence)
      }
      rec
  })
}

#' Write a GenBank file
#'
#' Writes sequences and features to a GenBank format file.
#'
#' @param file Path to the output file.
#' @param sequence Character string containing the sequence (DNA/RNA).
#' @param features Data frame of features. Must contain columns:
#'   \itemize{
#'     \item \code{key}: Feature type (e.g., "CDS", "gene")
#'     \item \code{location}: Location string (e.g., "1..100")
#'     \item \code{qualifiers}: List column of named character vectors
#'   }
#' @param metadata Named list of metadata (optional). Supported fields:
#'   \itemize{
#'     \item \code{name} (Locus Name)
#'     \item \code{definition}
#'     \item \code{accession}
#'     \item \code{version}
#'     \item \code{keywords} (character vector)
#'     \item \code{source}
#'     \item \code{organism}
#'     \item \code{molecule_type} (e.g., "DNA")
#'     \item \code{division}
#'     \item \code{topology} ("linear" or "circular")
#'     \item \code{date} (format: \code{DD-MON-YYYY})
#'     \item \code{references} (list of references; each reference may include
#'       \code{description}, \code{authors}, \code{consortium}, \code{title},
#'       \code{journal}, \code{pubmed}, \code{remark})
#'   }
#' @return Logical TRUE on success.
#' @export
#' @examples
#' \dontrun{
#'   meta <- list(definition = "Example Sequence", accession = "AB0001")
#'   feats <- data.frame(
#'     key = "source", 
#'     location = "1..100", 
#'     qualifiers = I(list(c(organism = "Homo sapiens")))
#'   )
#'   write_genbank("out.gb", "ATGC...", feats, meta)
#' }
write_genbank <- function(file, sequence, features, metadata = list()) {
  stopifnot(is.character(sequence), length(sequence) == 1)
  stopifnot(is.data.frame(features))
  stopifnot(is.list(metadata))

  if (is.na(sequence) || !nzchar(sequence)) {
    stop("Sequence must be a non-empty character string.")
  }

  sequence <- toupper(sequence)

  required_meta <- c("definition", "accession")
  missing_meta <- setdiff(required_meta, names(metadata))
  if (length(missing_meta) > 0) {
    stop("Metadata must include: definition, accession.")
  }
  for (field in required_meta) {
    value <- metadata[[field]]
    if (!is.character(value) || length(value) != 1 || is.na(value) || !nzchar(value)) {
      stop("Metadata must include: definition, accession.")
    }
  }
  
  if (!all(c("key", "location", "qualifiers") %in% colnames(features))) {
    stop("Features data frame must contain 'key', 'location', and 'qualifiers' columns.")
  }

  # Ensure qualifiers is a list column
  if (!is.list(features$qualifiers)) {
      stop("Features 'qualifiers' column must be a list.")
  }
    if (length(features$key) != nrow(features) ||
      length(features$location) != nrow(features) ||
      length(features$qualifiers) != nrow(features)) {
      stop("Features 'key', 'location', and 'qualifiers' must have the same length.")
    }
  for (i in seq_along(features$qualifiers)) {
      q <- features$qualifiers[[i]]
      if (!is.character(q)) {
          stop("Each element of 'qualifiers' must be a named character vector.")
      }
      q_names <- names(q)
      if (is.null(q_names) || any(is.na(q_names)) || any(q_names == "")) {
          stop("Each element of 'qualifiers' must have non-empty names.")
      }
  }

    validate_location <- function(loc) {
      if (!is.character(loc) || length(loc) != 1 || is.na(loc) || !nzchar(loc)) {
        stop("Each feature 'location' must be a non-empty character string.")
      }
      loc_trim <- gsub("\\s+", "", loc)
      loc_no_keywords <- gsub("(?i)(complement|join|order|one-of)", "", loc_trim, perl = TRUE)
      loc_no_allowed <- gsub("[0-9<>,.()^]", "", loc_no_keywords)
      if (grepl("[A-Za-z]", loc_no_allowed)) {
        stop("Location string contains unsupported tokens.")
      }
    }
    for (loc in features$location) {
      validate_location(loc)
    }

  path <- normalizePath(file, mustWork = FALSE)
  res <- write_genbank_rs(path, sequence, features, metadata)
  if (is.character(res) && length(res) == 1) {
      stop(res)
  }
  res
}
