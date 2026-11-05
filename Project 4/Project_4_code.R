# Project 4: TCGA-LUAD NKX2-1 expression analysis
# Base R only: no add-on R packages are required.

options(stringsAsFactors = FALSE)

sample_sheet <- "gdc_sample_sheet.tsv"
output_plot  <- "NKX2-1_expression.png"
output_table <- "NKX2-1_expression_values.tsv"

# ---------- 1. Read and check the sample sheet ----------
if (!file.exists(sample_sheet)) {
  stop("Cannot find gdc_sample_sheet.tsv. Run this script from the folder containing the sample sheet.")
}

samples <- read.delim(
  sample_sheet,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

# Normalize column names so both 'Sample Type' and 'Sample.Type' are recognized.
normal_names <- tolower(gsub("[^[:alnum:]]", "", names(samples)))

find_column <- function(choices, description, required = TRUE) {
  for (choice in choices) {
    hit <- which(normal_names == choice)
    if (length(hit) > 0) return(names(samples)[hit[1]])
  }
  if (required) {
    stop(paste0(
      "Could not identify the ", description, " column in gdc_sample_sheet.tsv. Columns found: ",
      paste(names(samples), collapse = ", ")
    ))
  }
  NA_character_
}

sample_type_column <- find_column(c("sampletype"), "sample type")
file_name_column   <- find_column(c("filename"), "file name")
file_id_column     <- find_column(c("fileid", "fileuuid"), "file ID", required = FALSE)

sample_types <- trimws(as.character(samples[[sample_type_column]]))
file_names   <- trimws(as.character(samples[[file_name_column]]))

# Keep only the two sample groups and the STAR gene-count files for this analysis.
keep <- sample_types %in% c("Solid Tissue Normal", "Primary Tumor") &
        grepl("rna_seq\\.augmented_star_gene_counts\\.tsv$", file_names)
selected <- samples[keep, , drop = FALSE]

if (nrow(selected) == 0) {
  stop("No Primary Tumor or Solid Tissue Normal STAR gene-count rows were found in the sample sheet.")
}

selected_types <- trimws(as.character(selected[[sample_type_column]]))
selected_names <- trimws(as.character(selected[[file_name_column]]))
if (!is.na(file_id_column)) {
  selected_ids <- trimws(as.character(selected[[file_id_column]]))
} else {
  selected_ids <- rep(NA_character_, nrow(selected))
}

cat("Sample sheet loaded successfully.\n")
cat("Tumor/normal STAR-count rows to check:", nrow(selected), "\n\n")

# ---------- 2. Index the extracted STAR gene-count files ----------
count_files <- list.files(
  path = ".",
  pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE
)

if (length(count_files) == 0) {
  stop(paste0(
    "No STAR gene-count TSV files were found under the current folder.\n",
    "Extract tcga_data.tar.gz first with: tar -xzf tcga_data.tar.gz"
  ))
}

# Try the GDC file ID first (often present in the extracted folder path),
# then fall back to an exact filename match.
locate_count_file <- function(file_name, file_id, available_files) {
  if (!is.na(file_id) && nzchar(file_id)) {
    id_hits <- available_files[grepl(file_id, available_files, fixed = TRUE)]
    if (length(id_hits) > 0) {
      name_hits <- id_hits[basename(id_hits) == file_name]
      if (length(name_hits) > 0) return(name_hits[1])
      return(id_hits[1])
    }
  }

  exact_hits <- available_files[basename(available_files) == file_name]
  if (length(exact_hits) > 0) return(exact_hits[1])

  NA_character_
}

# ---------- 3. Read the NKX2-1 TPM from one gene-count file ----------
read_nkx2_tpm <- function(path) {
  text_lines <- readLines(path, warn = FALSE)

  # GDC STAR-count files can contain comment/metadata lines before the table.
  header_hits <- which(startsWith(text_lines, paste0("gene_id", "\t")))
  if (length(header_hits) == 0) {
    header_hits <- which(startsWith(text_lines, paste0("\ufeff", "gene_id", "\t")))
  }
  if (length(header_hits) == 0) {
    stop("Could not locate the gene_id header row.")
  }
  header_position <- header_hits[1]

  column_names <- strsplit(text_lines[header_position], "\t", fixed = TRUE)[[1]]
  column_names[1] <- sub("^\ufeff", "", column_names[1])
  gene_name_position <- match("gene_name", column_names)
  tpm_position <- match("tpm_unstranded", column_names)

  if (is.na(gene_name_position)) stop("The gene_name column is missing.")
  if (is.na(tpm_position)) stop("The tpm_unstranded column is missing.")

  # Match the gene symbol as its own tab-separated field, not as a substring.
  if (header_position >= length(text_lines)) {
    stop("No gene rows follow the header.")
  }
  gene_rows <- text_lines[(header_position + 1L):length(text_lines)]
  split_rows <- strsplit(gene_rows, "\t", fixed = TRUE)
  is_target_gene <- vapply(
    split_rows,
    function(fields) length(fields) >= gene_name_position && fields[gene_name_position] == "NKX2-1",
    logical(1)
  )
  gene_rows <- gene_rows[is_target_gene]
  if (length(gene_rows) == 0) stop("NKX2-1 was not found in this file.")

  fields <- strsplit(gene_rows[1], "\t", fixed = TRUE)[[1]]
  if (length(fields) < max(gene_name_position, tpm_position)) {
    stop("The NKX2-1 row has fewer fields than expected.")
  }
  if (fields[gene_name_position] != "NKX2-1") {
    stop("The matched gene row did not have the exact NKX2-1 symbol.")
  }

  value <- suppressWarnings(as.numeric(fields[tpm_position]))
  if (length(value) != 1 || !is.finite(value) || value < 0) {
    stop("The NKX2-1 TPM value is missing or invalid.")
  }
  value
}

# ---------- 4. Match each sample to its expression file ----------
found_paths <- rep(NA_character_, nrow(selected))
tpm_values  <- rep(NA_real_, nrow(selected))
status      <- rep("Not processed", nrow(selected))

for (i in seq_len(nrow(selected))) {
  path <- locate_count_file(selected_names[i], selected_ids[i], count_files)
  found_paths[i] <- path

  if (is.na(path)) {
    status[i] <- "Expression file not found"
    next
  }

  extraction <- tryCatch(
    list(value = read_nkx2_tpm(path), message = "OK"),
    error = function(e) list(
      value = NA_real_,
      message = paste("Could not read NKX2-1 TPM:", conditionMessage(e))
    )
  )
  tpm_values[i] <- extraction$value
  status[i] <- extraction$message

  if (i %% 50 == 0 || i == nrow(selected)) {
    cat("Checked", i, "of", nrow(selected), "samples\r")
  }
}
cat("\n")

expression_data <- data.frame(
  Sample_Type = selected_types,
  File_ID = selected_ids,
  File_Name = selected_names,
  TPM = tpm_values,
  Log2_TPM = ifelse(is.finite(tpm_values), log2(tpm_values + 1), NA_real_),
  File_Path = found_paths,
  Status = status,
  stringsAsFactors = FALSE
)

usable <- is.finite(expression_data$TPM) & is.finite(expression_data$Log2_TPM)
analysis_data <- expression_data[usable, , drop = FALSE]
analysis_data$Sample_Type <- factor(
  analysis_data$Sample_Type,
  levels = c("Solid Tissue Normal", "Primary Tumor")
)

# Keep a record of successfully extracted expression values and any failed rows.
write.table(
  expression_data,
  file = output_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

cat("\nExpression extraction summary:\n")
cat("  Selected samples:", nrow(selected), "\n")
cat("  Samples with usable NKX2-1 TPM:", nrow(analysis_data), "\n")
cat("  Samples skipped:", sum(!usable), "\n")
if (any(!usable)) {
  cat("  Detailed file/status information is saved in", output_table, "\n")
}

if (nrow(analysis_data) == 0) {
  stop("No usable NKX2-1 expression values were extracted. Check the file paths and the Status column in the TSV output.")
}

sample_counts <- table(analysis_data$Sample_Type)
print(sample_counts)

normal_n <- sum(analysis_data$Sample_Type == "Solid Tissue Normal")
tumor_n  <- sum(analysis_data$Sample_Type == "Primary Tumor")
if (normal_n < 2 || tumor_n < 2) {
  stop("At least two valid samples in each group are required for the two-sample t-test and boxplot comparison.")
}

# ---------- 5. Statistical comparison ----------
# Welch's t-test compares group means without assuming equal variances.
test_result <- stats::t.test(
  Log2_TPM ~ Sample_Type,
  data = analysis_data,
  var.equal = FALSE
)

cat("\nWelch two-sample t-test on log2(TPM + 1):\n")
print(test_result)

# ---------- 6. Save the requested boxplot ----------
p_value_label <- format(test_result$p.value, scientific = TRUE, digits = 3)

grDevices::png(filename = output_plot, width = 900, height = 700)
tryCatch(
  {
    graphics::par(las = 1)
    graphics::boxplot(
      Log2_TPM ~ Sample_Type,
      data = analysis_data,
      main = "TCGA-LUAD: NKX2-1 expression",
      xlab = "Sample type",
      ylab = "log2(TPM + 1)",
      col = c("lightblue", "pink")
    )
    graphics::mtext(paste("p =", p_value_label), side = 3, line = 0.2)
  },
  finally = grDevices::dev.off()
)

cat("\nPlot saved to:", output_plot, "\n")
cat("Expression table saved to:", output_table, "\n")
cat("Project 4 analysis finished.\n")
