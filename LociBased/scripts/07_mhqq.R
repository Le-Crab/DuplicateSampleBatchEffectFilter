#!/usr/bin/env Rscript

# Script: 07_mhqq.R
# Purpose: Turn all PLINK .assoc files into Manhattan and QQ plots.
# Inputs:
#   - WD/outputs/gwas/*.assoc  (produced by 05_bcf2pca.sh)
# Outputs:
#   - WD/outputs/plots/*.jpg

cat("SUBSTEP 0: START\n")

#   CUSTOMIZATION
#   1. qqplot Y-lim - set every time, default == 16
# qqplot_ylim <- 16
#   2. Upper threshold for QQplot p-values phred saled
qqplot_ylim_switch <- FALSE
qqplot_phred_threshold <- 10
gwas_threshold <- 1e-8
suggestive_threshold <- 5e-8


#  STRUCTURE
#   1. Define functions for plotting and some stats calculation
#   2. Execution

cat("SUBSTEP 1: LOADING R PACKAGES\n")
suppressPackageStartupMessages({
    if (!require(CMplot, quietly = TRUE)) {
        stop("CMplot package is required. Install with: install.packages('CMplot')")
    }
    if (!require(dplyr, quietly = TRUE)) {
        stop("dplyr package is required. Install with: install.packages('dplyr')")
    }
})

# Determine repo root (WD) from script location.
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path) == 0) {
    script_dir <- getwd()
} else {
    script_dir <- dirname(normalizePath(script_path))
}
WD <- normalizePath(file.path(script_dir, ".."))

DIR_GWAS <- file.path(WD, "outputs", "gwas")
DIR_PLOTS <- file.path(WD, "outputs", "plots")
dir.create(DIR_PLOTS, recursive = TRUE, showWarnings = FALSE)

cat(paste0("WD: ", WD, "\n"))
cat(paste0("GWAS dir: ", DIR_GWAS, "\n"))
cat(paste0("Plots dir: ", DIR_PLOTS, "\n"))

# Function to read and validate GWAS data
read_gwas_data <- function(file_path) {
    if (!file.exists(file_path)) {
        stop(paste("File not found:", file_path))
    }
    
    cat("Reading GWAS data from:", file_path, "\n")
    data <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
    
    # Validate required columns
    required_cols <- c("SNP", "CHR", "BP", "P")
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
        stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    }
    
    # Clean and validate data
    data_clean <- data[, required_cols]
    data_clean <- data_clean[complete.cases(data_clean), ]
    data_clean <- data_clean[data_clean$P >= 0 & data_clean$P <= 1, ]
    return(data_clean)
}

# Function to calculate genome-wide significance statistics
calculate_gwas_stats <- function(data, label = "Dataset") {

    
    total_snps <- nrow(data)
    significant_snps <- sum(data$P < gwas_threshold, na.rm = TRUE)
    suggestive_snps <- sum(data$P < suggestive_threshold, na.rm = TRUE)
    
    cat("\n", label, "Statistics:\n")
    cat("  Total SNPs:", total_snps, "\n")
    cat("  Genome-wide significant for RARE variants (P < 1e-8):", significant_snps, "\n")
    cat("  Genome-wide significant for COMMON variants(P < 5e-8):", suggestive_snps, "\n")
    cat("  Significance rate:", round(significant_snps/total_snps * 100, 4), "%\n")
    
    return(list(
        total = total_snps,
        significant = significant_snps,
        suggestive = suggestive_snps
    ))
}

# Function to create Manhattan plot
create_manhattan_plot <- function(data, output_name, plot_title = "Manhattan Plot") {
    cat("Creating Manhattan plot:", output_name, "\n")
    
    CMplot(data, 
           plot.type = "m",
           col = c("grey30", "grey60"),
           cex = 0.3,
           band = 1,
           LOG10 = TRUE,
           verbose = FALSE,
           width = 15,
           height = 6,
           main = plot_title,
           dpi = 300,
           ylim = c(0, 16),
           signal.col = c("#ff7300", "#6c3d88"),
           signal.cex = c(0.3, 0.3),
           amplify = TRUE,
           threshold = c(gwas_threshold, suggestive_threshold),
           threshold.lwd = c(1, 1),
           threshold.col = c("black", "grey"),
           threshold.lty = c(1, 2),
           file.output = TRUE,
           file = "jpg",
           file.name = output_name)
}

# Function to create QQ plot
create_qq_plot <- function(data, output_name, plot_title = "QQ Plot", ylim_max = qqplot_phred_threshold) {
    cat("Creating QQ plot:", output_name, "\n")
    # adjust the behavor of QQplotting according to Y lim handling stated on top of script
    if (qqplot_ylim_switch == TRUE) {
        clip_p <- 10^(-ylim_max)
        data_ylim_considered <- data
        if (is.data.frame(data_ylim_considered) && "P" %in% names(data_ylim_considered)) {
            data_ylim_considered$P <- pmax(data_ylim_considered$P, clip_p)
        } else if (is.numeric(data_ylim_considered)) {
        data_ylim_considered <- pmax(data_ylim_considered, clip_p)
        }
    } else {
        data_ylim_considered <- data
    }

    CMplot(data_ylim_considered,
           plot.type = "q",
           box = FALSE,
           file = "jpg",
           cex = 0.35,
           file.name = output_name,
           dpi = 300,
           #ylim = c(0, ylim_max),
           conf.int = TRUE,
           conf.int.col = NULL,
           threshold.col = "red",
           threshold.lty = 2,
           file.output = TRUE,
           verbose = FALSE,
           width = 5,
           height = 5,
           main = plot_title)
}

# Function to find Removed SNPs between datasets
find_removed_snps <- function(control_data, filtered_data) {
    if (is.null(control_data) || is.null(filtered_data)) {
        return(NULL)
    }
    
    removed <- setdiff(control_data$SNP, filtered_data$SNP)
    cat("removed SNPs found:", length(removed), "\n")
    return(removed)
}


# ========================== MAIN ========================== #
cat("\nSUBSTEP 2: DISCOVERING .assoc FILES\n")

if (!dir.exists(DIR_GWAS)) {
    stop(paste0("[e] GWAS directory not found: ", DIR_GWAS))
}

assoc_files <- list.files(path = DIR_GWAS, pattern = "\\.assoc$", full.names = TRUE)

if (length(assoc_files) == 0) {
    stop(paste0("[e] No .assoc files found in: ", DIR_GWAS))
}

cat("Found .assoc files:\n")
for (f in assoc_files) cat(paste0("  - ", basename(f), "\n"))

# Initialize variables
control_data <- NULL
filtered_data <- NULL

# Process files based on naming convention
for (assoc_file in assoc_files) {
    if (grepl("control", basename(assoc_file), ignore.case = TRUE)) {
        control_data <- read_gwas_data(assoc_file)
        control_file_name <- sub("\\.assoc$", "", basename(assoc_file))
        calculate_gwas_stats(control_data, "Control")
    } else if (grepl("Filtered|LociBased", basename(assoc_file), ignore.case = TRUE)) {
        filtered_data <- read_gwas_data(assoc_file)
        filtered_file_name <- sub("\\.assoc$", "", basename(assoc_file))
        calculate_gwas_stats(filtered_data, "Filtered")
    } else {
            # If naming convention is unclear, use first file as control, second as prediction
            if (is.null(control_data)) {
                control_data <- read_gwas_data(assoc_file)
                control_file_name <- sub("\\.assoc$", "", basename(assoc_file))
                calculate_gwas_stats(control_data, "Control")
            } else if (is.null(filtered_data)) {
                filtered_data <- read_gwas_data(assoc_file)
                filtered_file_name <- sub("\\.assoc$", "", basename(assoc_file))
                calculate_gwas_stats(filtered_data, "Filtered")
            }
        }
}
cat("\nSUBSTEP 3: DATASETS\n")
if (!is.null(control_data))  cat("  - Control dataset:", control_file_name , "\n")
if (!is.null(filtered_data)) cat("  - Filtered dataset:", filtered_file_name , "\n")
# Find removed SNPs if both datasets available
removed_snps <- find_removed_snps(control_data, filtered_data)

cat("\nSUBSTEP 4: GENERATING PLOTS\n")

# Generate plots
if (!is.null(control_data)) {
    # Control plots
    oldwd <- getwd()
    setwd(DIR_PLOTS)
    create_manhattan_plot(control_data, control_file_name, paste0(control_file_name))
    create_qq_plot(control_data, control_file_name, paste0(control_file_name))
    setwd(oldwd)
}

if (!is.null(filtered_data)) {
    # Prediction plots
    oldwd <- getwd()
    setwd(DIR_PLOTS)
    create_manhattan_plot(filtered_data, filtered_file_name, paste0(filtered_file_name))
    create_qq_plot(filtered_data, filtered_file_name, paste0(filtered_file_name))
    setwd(oldwd)
}

# Summary comparison if both datasets available
if (!is.null(control_data) && !is.null(filtered_data)) {
    control_stats <- calculate_gwas_stats(control_data, "")
    filtered_stats <- calculate_gwas_stats(filtered_data, "")
    
    if (control_stats$significant > 0 && filtered_stats$significant > 0) {
        reduction_ratio <- round(control_stats$significant / filtered_stats$significant, 1)
        cat("\nSignificant SNP reduction ratio:", reduction_ratio, "times\n")
    }
}

cat("\nScript finished. Plotting completed successfully!\n")