# Megablast via NCBI BLAST API
# Requires: Biostrings, httr, xml2
# Uses NCBI's web BLAST - no local BLAST+ installation needed

# library(Biostrings)
# library(httr)
# library(xml2)

#' Parse BLAST XML2 result into a data frame
#'
#' @param text Character string of XML response
#' @param query_id Query sequence ID
#' @return Data frame of hits, or NULL if no hits
#'
.parse_blast_xml <- function(text, query_id) {
  if (!requireNamespace("xml2", quietly = TRUE)) return(NULL)
  doc <- tryCatch(xml2::read_xml(text), error = function(e) NULL)
  if (is.null(doc)) return(NULL)

  # Support both BlastOutput (classic) and BlastOutput2 (XML2) structures
  hits <- xml2::xml_find_all(doc, "//*[local-name()='Hit']")
  if (length(hits) == 0) return(NULL)

  rows <- lapply(hits, function(h) {
    hsp <- xml2::xml_find_first(h, ".//*[local-name()='Hsp']")
    if (inherits(hsp, "xml_missing")) return(NULL)
    id_el <- xml2::xml_find_first(h, ".//*[local-name()='Hit_id']")
    def_el <- xml2::xml_find_first(h, ".//*[local-name()='Hit_def']")
    ident_el <- xml2::xml_find_first(hsp, ".//*[local-name()='Hsp_identity']")
    len_el <- xml2::xml_find_first(hsp, ".//*[local-name()='Hsp_align-len']")
    eval_el <- xml2::xml_find_first(hsp, ".//*[local-name()='Hsp_evalue']")
    bit_el <- xml2::xml_find_first(hsp, ".//*[local-name()='Hsp_bit-score']")
    qf_el <- xml2::xml_find_first(hsp, ".//*[local-name()='Hsp_query-from']")
    qt_el <- xml2::xml_find_first(hsp, ".//*[local-name()='Hsp_query-to']")
    sf_el <- xml2::xml_find_first(hsp, ".//*[local-name()='Hsp_hit-from']")
    st_el <- xml2::xml_find_first(hsp, ".//*[local-name()='Hsp_hit-to']")
    gaps_el <- xml2::xml_find_first(hsp, ".//*[local-name()='Hsp_gaps']")

    if (inherits(ident_el, "xml_missing") || inherits(len_el, "xml_missing")) return(NULL)
    pident <- as.numeric(xml2::xml_text(ident_el)) / as.numeric(xml2::xml_text(len_el)) * 100

    data.frame(
      query_id = query_id,
      subject_id = if (!inherits(id_el, "xml_missing")) xml2::xml_text(id_el) else NA,
      subject_def = if (!inherits(def_el, "xml_missing")) xml2::xml_text(def_el) else NA,
      pident = pident,
      length = as.integer(xml2::xml_text(len_el)),
      mismatch = if (!inherits(gaps_el, "xml_missing")) as.integer(xml2::xml_text(gaps_el)) else NA_integer_,
      evalue = if (!inherits(eval_el, "xml_missing")) as.numeric(xml2::xml_text(eval_el)) else NA_real_,
      bitscore = if (!inherits(bit_el, "xml_missing")) as.numeric(xml2::xml_text(bit_el)) else NA_real_,
      qstart = if (!inherits(qf_el, "xml_missing")) as.integer(xml2::xml_text(qf_el)) else NA_integer_,
      qend = if (!inherits(qt_el, "xml_missing")) as.integer(xml2::xml_text(qt_el)) else NA_integer_,
      sstart = if (!inherits(sf_el, "xml_missing")) as.integer(xml2::xml_text(sf_el)) else NA_integer_,
      send = if (!inherits(st_el, "xml_missing")) as.integer(xml2::xml_text(st_el)) else NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows[!sapply(rows, is.null)])
  if (is.null(out) || nrow(out) == 0) NULL else out
}

#' Parse BLAST CSV/TSV result (fallback when XML not available)
#'
.parse_blast_tabular <- function(text, query_id) {
  # Extract content from <pre> if present (NCBI wraps in HTML)
  pre_m <- regexec("<pre[^>]*>([\\s\\S]*?)</pre>", text)
  pre_grp <- regmatches(text, pre_m)[[1]]
  if (length(pre_grp) >= 2) text <- pre_grp[2]
  else { text <- gsub("<[^>]+>", "\n", text); text <- gsub("&nbsp;|&quot;", " ", text) }

  lines <- strsplit(text, "\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  # NCBI Text tabular: "# Fields: query acc.ver, subject acc.ver, % identity, ..." (comma sep)
  # Data rows are TAB-separated - don't use comma!
  header_line <- NULL
  data_lines <- character(0)
  for (i in seq_along(lines)) {
    if (grepl("^# Fields:", lines[i], ignore.case = TRUE)) {
      header_line <- sub("^# Fields:\\s*", "", lines[i], ignore.case = TRUE)
    } else if (grepl("^#", lines[i]) || grepl("^Status=|^RID |QBlastInfo", lines[i])) {
      next
    } else if (grepl("\t", lines[i]) && (grepl("^[A-Za-z0-9_]+\\t", lines[i]) || grepl("^[A-Za-z0-9.]+\\t", lines[i]))) {
      data_lines <- c(data_lines, lines[i])
    }
  }

  if (is.null(header_line) || length(data_lines) == 0) return(NULL)

  # Parse tab-separated data (12 cols: query, subject, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)
  std_cols <- c("query_id", "subject_id", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  con <- textConnection(paste(data_lines, collapse = "\n"))
  tbl <- tryCatch(
    utils::read.delim(con, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                     col.names = std_cols, check.names = FALSE),
    error = function(e) NULL,
    finally = close(con)
  )
  if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
  if (ncol(tbl) > 12) tbl <- tbl[, seq_len(12), drop = FALSE]
  # Query id is already first column in BLAST tabular output.
  # Keep parsed value when present to avoid duplicate `query_id` names.
  if (!"query_id" %in% names(tbl)) tbl$query_id <- query_id
  tbl
}

#' Fetch taxonomy for BLAST subject accessions (optional)
#'
#' Uses the NCBI EUtils API via the `rentrez` package to look up:
#' - organism (scientific name) and taxid from `nuccore` summaries
#' - lineage from `taxonomy` summaries (v1.0 XML, if available)
#'
#' @param accessions Character vector of accession.version strings (e.g. \"MZ409288.1\")
#' @param verbose If TRUE, prints progress
#' @return data.frame(subject_id, organism, taxid, lineage) or NULL if `rentrez` not installed
#'
.fetch_taxonomy <- function(accessions, verbose = TRUE) {
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    if (verbose) message("  Install taxonomy helper: install.packages('rentrez')")
    return(NULL)
  }
  accessions <- unique(na.omit(accessions))
  if (length(accessions) == 0) return(NULL)

  tax_df <- data.frame(
    subject_id = accessions,
    organism = NA_character_,
    taxid = NA_character_,
    lineage = NA_character_,
    stringsAsFactors = FALSE
  )

  # Batch fetch nuccore summaries (fewer HTTP calls)
  batch_size <- 50
  for (start in seq(1, length(accessions), by = batch_size)) {
    batch <- accessions[seq(start, min(start + batch_size - 1, length(accessions)))]
    ok <- FALSE
    tryCatch({
      Sys.sleep(0.34) # be gentle to NCBI
      summ_list <- rentrez::entrez_summary(db = "nuccore", id = batch, always_return_list = TRUE)
      # Map by returned accessionversion/caption, not by list order.
      for (s in summ_list) {
        if (is.null(s)) next
        acc <- s[["accessionversion"]]
        if (is.null(acc)) acc <- s[["caption"]]
        if (is.null(acc)) next
        idx <- which(tax_df$subject_id == acc)
        if (!length(idx)) next
        org <- s[["organism"]]
        tid <- s[["taxid"]]
        tax_df$organism[idx] <- if (!is.null(org)) org else NA_character_
        tax_df$taxid[idx] <- if (!is.null(tid)) as.character(tid) else NA_character_
      }
      ok <- TRUE
    }, error = function(e) { ok <<- FALSE })

    # Fallback per accession for this batch if batch call failed/throttled.
    if (!ok) {
      for (acc in batch) {
        tryCatch({
          Sys.sleep(0.5)
          s <- rentrez::entrez_summary(db = "nuccore", id = acc)
          if (is.null(s)) next
          idx <- which(tax_df$subject_id == acc)
          if (!length(idx)) next
          org <- s[["organism"]]
          tid <- s[["taxid"]]
          tax_df$organism[idx] <- if (!is.null(org)) org else NA_character_
          tax_df$taxid[idx] <- if (!is.null(tid)) as.character(tid) else NA_character_
        }, error = function(e) NULL)
      }
    }
    if (verbose) message("  Taxonomy (nuccore): ", min(start + batch_size - 1, length(accessions)), "/", length(accessions))
  }

  # Second-pass fallback:
  # Retry unresolved accessions by searching accession -> UID, then summary by UID.
  unresolved <- which(is.na(tax_df$taxid) | is.na(tax_df$organism))
  if (length(unresolved) > 0) {
    max_rounds <- 3
    if (verbose) message("  Starting taxonomy fallback for ", length(unresolved), " unresolved accessions...")
    for (round_idx in seq_len(max_rounds)) {
      still_unresolved <- which(is.na(tax_df$taxid) | is.na(tax_df$organism))
      if (length(still_unresolved) == 0) break
      matched_this_round <- 0L
      for (idx in still_unresolved) {
        acc <- tax_df$subject_id[idx]
        tryCatch({
          Sys.sleep(0.5)
          srch <- rentrez::entrez_search(db = "nuccore", term = paste0(acc, "[ACCN]"), retmax = 1)
          if (length(srch$ids) == 0) next
          s <- rentrez::entrez_summary(db = "nuccore", id = srch$ids[[1]])
          if (is.null(s)) next
          org <- s[["organism"]]
          tid <- s[["taxid"]]
          had_missing <- is.na(tax_df$organism[idx]) || is.na(tax_df$taxid[idx])
          if (!is.null(org) && nzchar(as.character(org))) tax_df$organism[idx] <- as.character(org)
          if (!is.null(tid) && nzchar(as.character(tid))) tax_df$taxid[idx] <- as.character(tid)
          if (had_missing && !(is.na(tax_df$organism[idx]) || is.na(tax_df$taxid[idx]))) {
            matched_this_round <- matched_this_round + 1L
          }
        }, error = function(e) NULL)
      }
      if (verbose) {
        remain <- sum(is.na(tax_df$taxid) | is.na(tax_df$organism))
        message("  Taxonomy fallback round ", round_idx, ": matched ", matched_this_round, ", remaining ", remain)
      }
      if (matched_this_round == 0) break
    }
  }

  # Lineage lookup (best-effort)
  if (requireNamespace("xml2", quietly = TRUE)) {
    taxids <- unique(na.omit(tax_df$taxid))
    if (length(taxids) > 0) {
      lineage_map <- setNames(rep(NA_character_, length(taxids)), taxids)
      for (j in seq_along(taxids)) {
        tryCatch({
          Sys.sleep(0.34)
          ts_xml <- rentrez::entrez_summary(db = "taxonomy", id = taxids[j], version = "1.0", retmode = "xml")
          doc <- tryCatch(xml2::read_xml(as.character(ts_xml)), error = function(e) NULL)
          if (!is.null(doc)) {
            ln <- xml2::xml_text(xml2::xml_find_first(doc, "//*[local-name()='Lineage']"))
            if (!is.na(ln) && nzchar(ln)) lineage_map[[taxids[j]]] <- ln
          }
        }, error = function(e) NULL)
      }
      tax_df$lineage <- unname(lineage_map[tax_df$taxid])
    }
  } else if (verbose) {
    message("  (Optional) Install xml2 to include lineage: install.packages('xml2')")
  }

  tax_df
}

#' Run megablast on a FASTA file and export results
#'
#' Submits sequences to NCBI BLAST megablast and retrieves tabular results.
#' NCBI guidelines: poll no more than once per minute per RID.
#'
#' @param fasta_path Path to input FASTA file
#' @param output_path Path for output CSV (default: _megablast_results.csv in same dir as input)
#' @param database BLAST database (default: "nt" for nucleotide)
#' @param hitlist_size Max hits per query (default: 50)
#' @param email Your email - NCBI requests this for problem contact
#' @param tool Tool name for NCBI (default: "seaslug_barcoding")
#' @param max_sequences Max sequences to process (NULL = all; use e.g. 3 for testing)
#' @param verbose If TRUE, print status on each poll attempt
#' @param debug_save_raw If TRUE, save raw NCBI response to file when parsing fails
#' @param include_taxonomy If TRUE, add organism/taxid/lineage columns (requires `rentrez`)
#' @param temp_dir Directory to store per-query temporary CSV files (default: output folder/megablast_temp)
#' @return Data frame of BLAST hits, also written to output_path
#'
run_megablast <- function(fasta_path,
                          output_path = NULL,
                          database = "nt",
                          hitlist_size = 50,
                          email = "user@example.com",
                          tool = "seaslug_barcoding",
                          max_sequences = NULL,
                          verbose = TRUE,
                          debug_save_raw = FALSE,
                          include_taxonomy = TRUE,
                          temp_dir = NULL) {

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Please install Biostrings: BiocManager::install('Biostrings')")
  }
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Please install httr: install.packages('httr')")
  }

  # Read FASTA
  seqs <- Biostrings::readDNAStringSet(fasta_path)
  if (length(seqs) == 0) stop("No sequences found in FASTA file.")
  if (!is.null(max_sequences)) seqs <- seqs[seq_len(min(max_sequences, length(seqs)))]

  # Default output path
  if (is.null(output_path)) {
    base <- sub("\\.(fasta|fa|fna)$", "", fasta_path, ignore.case = TRUE)
    output_path <- paste0(base, "_megablast_results.csv")
  }
  if (is.null(temp_dir)) {
    out_dir <- dirname(output_path)
    temp_dir <- file.path(out_dir, "megablast_temp")
  }
  if (!dir.exists(temp_dir)) dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  blast_url <- "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
  temp_files <- character(0)

  for (i in seq_along(seqs)) {
    qname <- names(seqs)[i]
    qseq <- as.character(seqs[i])
    query_fasta <- paste0(">", qname, "\n", qseq)

    message("Submitting: ", qname, " (", i, "/", length(seqs), ")")

    # Submit job (megablast: PROGRAM=blastn, MEGABLAST=on)
    put_resp <- httr::POST(
      blast_url,
      body = list(
        CMD = "Put",
        PROGRAM = "blastn",
        MEGABLAST = "on",
        DATABASE = database,
        QUERY = query_fasta,
        HITLIST_SIZE = hitlist_size,
        FORMAT_TYPE = "Text",
        ALIGNMENT_VIEW = "Tabular",
        DESCRIPTIONS = hitlist_size,
        ALIGNMENTS = hitlist_size,
        email = email,
        tool = tool
      ),
      encode = "form"
    )

    put_text <- httr::content(put_resp, as = "text", encoding = "UTF-8")

    # Extract RID from response
    rid <- regmatches(put_text, regexpr("RID = ([A-Z0-9]+)", put_text, perl = TRUE))
    if (length(rid) == 0) {
      warning("No RID returned for ", qname, ". Skipping.")
      next
    }
    rid <- sub("RID = ", "", rid)

    # Extract RTOE (estimated time in seconds)
    rtoe <- regmatches(put_text, regexpr("RTOE = ([0-9]+)", put_text, perl = TRUE))
    wait_sec <- if (length(rtoe) > 0) as.numeric(sub("RTOE = ", "", rtoe)) else 30
    wait_sec <- min(max(wait_sec, 20), 120)

    message("  Waiting ", wait_sec, " sec (RID: ", rid, ")")
    Sys.sleep(wait_sec)

    # Poll until ready: use minimal GET for status check (NCBI: poll max once per minute)
    # When requesting CSV directly, the response format can differ and Status=READY may not appear
    max_attempts <- 20   # ~20 min max per sequence
    poll_interval <- 65
    for (attempt in seq_len(max_attempts)) {
      # Step 1: Status check only (minimal params - gets response with Status= in comment)
      status_resp <- tryCatch(
        httr::GET(
          blast_url,
          query = list(CMD = "Get", RID = rid),
          httr::timeout(30)
        ),
        error = function(e) NULL
      )
      if (is.null(status_resp) || httr::status_code(status_resp) != 200) {
        if (verbose) message("  Poll ", attempt, ": Request failed, retrying...")
        Sys.sleep(poll_interval)
        next
      }
      status_text <- httr::content(status_resp, as = "text", encoding = "UTF-8")

      # Extract status (format: "Status=READY" or "Status=WAITING" in QBlastInfo block)
      status <- if (grepl("Status\\s*=\\s*READY", status_text, ignore.case = TRUE)) {
        "READY"
      } else if (grepl("Status\\s*=\\s*FAILED|Status\\s*=\\s*UNKNOWN|Status\\s*=\\s*Error", status_text, ignore.case = TRUE)) {
        "FAILED"
      } else {
        "WAITING"
      }

      if (verbose) {
        message("  Poll ", attempt, "/", max_attempts, ": Status=", status)
      }

      if (status == "READY") {
        # Step 2: Fetch results - XML2 is most reliable, then CSV, then Text
        tbl <- NULL
        for (fmt in c("XML2", "CSV", "Text")) {
          qry <- list(CMD = "Get", RID = rid, FORMAT_TYPE = fmt,
                     DESCRIPTIONS = hitlist_size, ALIGNMENTS = hitlist_size)
          if (fmt != "XML2") qry$ALIGNMENT_VIEW <- "Tabular"
          result_resp <- httr::GET(blast_url, query = qry, httr::timeout(60))
          result_text <- httr::content(result_resp, as = "text", encoding = "UTF-8")
          tbl <- if (fmt == "XML2") {
            .parse_blast_xml(result_text, qname)
          } else {
            .parse_blast_tabular(result_text, qname)
          }
          if (!is.null(tbl) && nrow(tbl) > 0) {
            # Save each query result first, then merge all files later.
            tmp_file <- file.path(temp_dir, paste0(sprintf("%03d", i), "_", qname, ".csv"))
            utils::write.csv(tbl, tmp_file, row.names = FALSE)
            temp_files <- c(temp_files, tmp_file)
            if (verbose && fmt != "XML2") message("  Parsed results (", fmt, " format)")
            break
          }
        }
        if (is.null(tbl) || nrow(tbl) == 0) {
          if (grepl("No hits found|no hits found|0 hits", result_text, ignore.case = TRUE)) {
            if (verbose) message("  No BLAST hits found for ", qname)
          } else {
            if (debug_save_raw) {
              raw_file <- paste0("blast_raw_", qname, "_", rid, ".txt")
              writeLines(result_text, raw_file)
              warning("Could not parse results for ", qname, ". Raw response saved to ", raw_file)
            } else {
              warning("Could not parse BLAST results for ", qname, ". RID ", rid,
                      " - use debug_save_raw=TRUE to save raw response, or view at https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=", rid)
            }
          }
        }
        break
      }
      if (status == "FAILED") {
        warning("BLAST job failed for ", qname, ". Check RID ", rid, " at https://blast.ncbi.nlm.nih.gov/")
        break
      }

      Sys.sleep(poll_interval)
    }
    if (attempt >= max_attempts && status == "WAITING") {
      warning("Timeout (", max_attempts, " polls) for ", qname, ". RID ", rid,
              " - check manually at https://blast.ncbi.nlm.nih.gov/")
    }
  }

  # Merge all per-query temp files
  temp_files <- unique(temp_files[file.exists(temp_files)])
  if (length(temp_files) == 0) {
    message("No results to export.")
    return(NULL)
  }

  out_list <- lapply(temp_files, function(fp) {
    tryCatch(utils::read.csv(fp, stringsAsFactors = FALSE), error = function(e) NULL)
  })
  out_list <- out_list[!sapply(out_list, is.null)]
  if (length(out_list) == 0) {
    message("No results to export.")
    return(NULL)
  }
  out <- do.call(rbind, out_list)

  if (include_taxonomy && "subject_id" %in% names(out)) {
    if (verbose) message("Fetching taxonomy for ", length(unique(out$subject_id)), " unique accessions...")
    tax_df <- .fetch_taxonomy(out$subject_id, verbose = verbose)
    if (!is.null(tax_df) && nrow(tax_df) > 0) {
      out$organism <- tax_df$organism[match(out$subject_id, tax_df$subject_id)]
      out$taxid <- tax_df$taxid[match(out$subject_id, tax_df$subject_id)]
      out$lineage <- tax_df$lineage[match(out$subject_id, tax_df$subject_id)]
    }
  }

  utils::write.csv(out, output_path, row.names = FALSE)
  message("Results saved to: ", output_path)
  if (verbose) message("Per-query temporary files saved in: ", temp_dir)
  return(out)
}

#' Test parser with sample NCBI BLAST XML (no network required)
test_megablast_parser <- function() {
  sample_xml <- '<?xml version="1.0"?>
<BlastOutput><BlastOutput_iterations>
<Iteration><Iteration_hits>
<Hit><Hit_num>1</Hit_num><Hit_id>NR_123456.1</Hit_id>
<Hit_def>Nudibranch 16S ribosomal RNA</Hit_def>
<Hit_hsps><Hsp><Hsp_identity>400</Hsp_identity><Hsp_align-len>410</Hsp_align-len>
<Hsp_evalue>1e-150</Hsp_evalue><Hsp_bit-score>550</Hsp_bit-score>
<Hsp_query-from>1</Hsp_query-from><Hsp_query-to>400</Hsp_query-to>
<Hsp_hit-from>50</Hsp_hit-from><Hsp_hit-to>459</Hsp_hit-to><Hsp_gaps>2</Hsp_gaps></Hsp></Hit_hsps></Hit>
</Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>'
  tbl <- .parse_blast_xml(sample_xml, "NUD003_16S")
  stopifnot(!is.null(tbl), nrow(tbl) == 1, tbl$subject_id[1] == "NR_123456.1",
            tbl$pident[1] > 97, tbl$bitscore[1] == 550)
  message("Parser test OK")
  invisible(tbl)
}

# --- Test / Example ---
# Run: test_megablast_parser() to verify parser (no network)
# Run full test: Sys.setenv(RUN_MEGABLAST_TEST="1"); source("5_megablast_fasta.R")
if (Sys.getenv("RUN_MEGABLAST_TEST") == "1") {
  message("Running parser unit test...")
  test_megablast_parser()
  message("Running megablast with example FASTA (1 seq)...")
  setwd("D:/MSL Miscellaneous/20240408 BGI Pipeline Github")
  results <- run_megablast(
    fasta_path = "example_data/16S_trimmed_seq_full.fasta",
    output_path = "example_data/16S_megablast_results_full.csv",
    email = "user@example.com",
    max_sequences = 5,
    hitlist_size = 20,
    verbose = TRUE,
    debug_save_raw = TRUE
  )
  if (!is.null(results)) {
    message("Test OK. First few rows:")
    print(head(results, 3))
  }
}
