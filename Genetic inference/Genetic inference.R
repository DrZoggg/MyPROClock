# =====================================================================
# ADVANCED MENDELIAN RANDOMIZATION ANALYSIS PIPELINE
# =====================================================================
# Comprehensive MR Analysis Framework for Causal Inference Studies
# Incorporates: 
#   - Instrumental Variable Selection
#   - Harmonization Procedures
#   - Multiple MR Methods
#   - Sensitivity Analyses
#   - Visualization Suite
# =====================================================================

library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(viridis)
library(data.table)
library(tidyr)
library(qqman)
library(forestplot)

# =====================================================================
# DATA PREPARATION MODULE
# =====================================================================

#' Process eQTL Data from GTEx Portal
#' 
#' @param file_path Path to GTEx Portal CSV file
#' @return Formatted eQTL data for MR analysis
process_gtex_eqtl_data <- function(file_path) {
    # Read and validate input data
    if (!file.exists(file_path)) {
        stop("GTEx data file not found at specified path")
    }
    
    eqtl_data <- fread(file_path) %>%
        mutate(
            chr = gsub("_.*", "", Variant.Id),
            pos = as.integer(gsub(".*_(\\d+)_.*", "\\1", Variant.Id)),
            other_allele = gsub(".*_(\\w)_(\\w)_b38", "\\1", Variant.Id),
            effect_allele = gsub(".*_(\\w)_(\\w)_b38", "\\2", Variant.Id),
            se = abs(beta / qnorm(1 - pval / 2))
        ) %>%
        select(
            SNP = SNP.Id, 
            chr, 
            pos, 
            other_allele, 
            beta = NES,
            effect_allele,
            pval = P.Value,
            gene_symbol = Gene.Symbol
        ) %>%
        mutate(
            id.exposure = "GTEx_eQTL",
            chr = as.character(15)  # Force chromosome to 15 for MFAP4
        )
    
    return(eqtl_data)
}

#' Process GWAS VCF Data
#' 
#' @param file_path Path to GWAS VCF file
#' @return Formatted GWAS data for MR analysis
process_gwas_vcf_data <- function(file_path) {
    # Read and validate input data
    if (!file.exists(file_path)) {
        stop("GWAS VCF file not found at specified path")
    }
    
    gwas_data <- fread(file_path, skip = "#CHROM") %>%
        separate(
            `ukb-e-I21_CSA`, 
            into = c("beta", "se", "pval", "af"), 
            sep = ":", 
            convert = TRUE
        ) %>%
        rename(
            SNP = ID,
            effect_allele.outcome = ALT,
            other_allele.outcome = REF,
            beta.outcome = beta,
            se.outcome = se,
            pval.outcome = pval
        ) %>%
        mutate(
            outcome = "AcuteMyocardialInfarction",
            id.outcome = "I21",
            eaf.outcome = NA_real_
        ) %>%
        select(
            SNP, 
            effect_allele.outcome, 
            other_allele.outcome,
            beta.outcome, 
            se.outcome, 
            pval.outcome, 
            outcome,
            id.outcome,
            eaf.outcome
        )
    
    return(gwas_data)
}

# =====================================================================
# MR ANALYSIS CORE FUNCTIONS
# =====================================================================

#' Execute Full MR Workflow
#' 
#' @param exposure_id ID for exposure dataset
#' @param outcome_id ID for outcome dataset
#' @param pval_threshold Instrument selection p-value threshold
#' @return Comprehensive MR results object
execute_mr_workflow <- function(exposure_id, outcome_id, pval_threshold = 5e-8) {
    # Instrument extraction
    instruments <- extract_instruments(
        outcomes = exposure_id,
        p1 = pval_threshold
    )
    
    # Outcome data extraction
    outcome_data <- extract_outcome_data(
        snps = instruments$SNP,
        outcomes = outcome_id
    )
    
    # Data harmonization
    harmonized_data <- harmonise_data(
        exposure_dat = instruments,
        outcome_dat = outcome_data,
        action = 2
    )
    
    # MR analysis with multiple methods
    mr_results <- mr(
        harmonized_data,
        method_list = c(
            "mr_ivw",
            "mr_egger_regression",
            "mr_weighted_median",
            "mr_wald_ratio",
            "mr_simple_mode_nome",
            "mr_weighted_mode_nome"
        )
    ) %>%
        generate_odds_ratios()
    
    # Sensitivity analyses
    heterogeneity <- mr_heterogeneity(harmonized_data)
    pleiotropy <- mr_pleiotropy_test(harmonized_data)
    leaveoneout <- mr_leaveoneout(harmonized_data)
    scatter_plot <- mr_scatter_plot(mr_results, harmonized_data)
    
    # MR-PRESSO analysis
    presso_input <- harmonized_data %>%
        select(
            beta.exposure,
            beta.outcome,
            se.exposure,
            se.outcome
        )
    
    presso_results <- mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = presso_input,
        NbDistribution = 1000,
        SignifThreshold = 0.05
    )
    
    # Compile comprehensive results
    results <- list(
        instruments = instruments,
        harmonized_data = harmonized_data,
        mr_results = mr_results,
        heterogeneity = heterogeneity,
        pleiotropy = pleiotropy,
        leaveoneout = leaveoneout,
        scatter_plot = scatter_plot,
        presso_results = presso_results
    )
    
    return(results)
}

#' Perform Reverse MR Analysis
#' 
#' @param exposure_id ID for original exposure
#' @param outcome_id ID for original outcome
#' @return Reverse MR results
perform_reverse_mr <- function(exposure_id, outcome_id) {
    rev_instruments <- extract_instruments(outcomes = outcome_id)
    rev_outcome <- extract_outcome_data(
        snps = rev_instruments$SNP,
        outcomes = exposure_id
    )
    
    rev_harmonized <- harmonise_data(
        exposure_dat = rev_instruments,
        outcome_dat = rev_outcome,
        action = 2
    )
    
    rev_results <- mr(
        rev_harmonized,
        method_list = c(
            "mr_ivw",
            "mr_weighted_median",
            "mr_wald_ratio",
            "mr_simple_mode_nome",
            "mr_weighted_mode_nome"
        )
    )
    
    return(list(
        harmonized_data = rev_harmonized,
        mr_results = rev_results
    ))
}

# =====================================================================
# VISUALIZATION MODULE
# =====================================================================

#' Create F-Statistic Visualization
#' 
#' @param instrument_data Instrumental variable dataset
#' @return ggplot object of F-statistics
visualize_f_statistics <- function(instrument_data) {
    instrument_data <- instrument_data %>%
        mutate(
            F_stat = (beta.exposure^2) / (se.exposure^2),
            SNP_order = factor(SNP, levels = SNP[order(F_stat)])
    
    ggplot(instrument_data, aes(x = SNP_order, y = F_stat, fill = F_stat)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_c(option = "C", name = "F Value") +
        theme_minimal(base_size = 13) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank()
        ) +
        labs(
            title = "Instrument Strength Assessment (F-Statistics)",
            x = "Genetic Variant",
            y = "F Statistic Value"
        )
}

#' Create Directional Comparison Forest Plot
#' 
#' @param directional_results Main MR results
#' @param reverse_results Reverse MR results
#' @param analysis_name Name of analysis for title
#' @return forestplot object
create_direction_comparison_plot <- function(directional_results, reverse_results, analysis_name) {
    # Extract relevant estimates
    directional_est <- directional_results %>%
        filter(method == "Inverse variance weighted") %>%
        select(b, lo_ci, up_ci)
    
    reverse_est <- reverse_results %>%
        filter(method == "Inverse variance weighted") %>%
        select(b, lo_ci, up_ci)
    
    # Prepare data frame
    comparison_df <- data.frame(
        Analysis = c("Directional MR", "Reverse MR"),
        Estimate = c(directional_est$b, reverse_est$b),
        Lower_CI = c(directional_est$lo_ci, reverse_est$lo_ci),
        Upper_CI = c(directional_est$up_ci, reverse_est$up_ci)
    )
    
    # Prepare text labels
    table_text <- list(
        c("Analysis", "Estimate (95% CI)"),
        c("Directional MR", sprintf("%.4f (%.4f-%.4f)", 
            comparison_df$Estimate[1], 
            comparison_df$Lower_CI[1], 
            comparison_df$Upper_CI[1])),
        c("Reverse MR", sprintf("%.4f (%.4f-%.4f)", 
            comparison_df$Estimate[2], 
            comparison_df$Lower_CI[2], 
            comparison_df$Upper_CI[2]))
    )
    
    # Create plot
    forestplot(
        labeltext = table_text,
        mean = c(NA, comparison_df$Estimate),
        lower = c(NA, comparison_df$Lower_CI),
        upper = c(NA, comparison_df$Upper_CI),
        zero = 0,
        boxsize = 0.15,
        lwd.ci = 2.5,
        xlab = "Causal Effect Estimate",
        title = paste("Directional Comparison:", analysis_name),
        col = fpColors(
            box = c("#1f78b4", "#e31a1c"),
            line = c("#1f78b4", "#e31a1c"),
            zero = "gray40"
        ),
        txt_gp = fpTxtGp(
            label = gpar(fontsize = 11, fontface = "plain"),
            ticks = gpar(fontsize = 12),
            xlab = gpar(fontsize = 14, fontface = "bold"),
            title = gpar(fontsize = 15, fontface = "bold")
        ),
        graph.pos = 3,
        is.summary = c(TRUE, FALSE, FALSE),
        clip = c(min(comparison_df$Lower_CI) - 0.1, 
                 max(comparison_df$Upper_CI) + 0.1),
        xticks = pretty(c(comparison_df$Lower_CI, comparison_df$Upper_CI), n = 8)
    )
}

#' Create QTL Manhattan Plot
#' 
#' @param qtl_data Formatted QTL data
#' @param chromosome Target chromosome
#' @return ggplot object
visualize_qtl_manhattan <- function(qtl_data, chromosome) {
    plot_data <- qtl_data %>%
        mutate(
            logp = -log10(pval.exposure),
            chr = as.numeric(chr.exposure),
            pos = as.numeric(pos.exposure)
        ) %>%
        arrange(pos)
    
    ggplot(plot_data, aes(x = pos, y = logp)) +
        geom_point(aes(color = logp), alpha = 0.8, size = 2.5) +
        scale_color_viridis_c(option = "B", direction = -1, name = "-log10(p)") +
        geom_hline(
            yintercept = -log10(5e-8), 
            linetype = "dashed", 
            color = "red",
            size = 0.8
        ) +
        labs(
            x = paste("Genomic Position on Chromosome", chromosome),
            y = "-log10(p-value)",
            title = paste("QTL Association Mapping: Chromosome", chromosome),
            subtitle = "Red line indicates genome-wide significance threshold (5e-8)"
        ) +
        theme_minimal() +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linewidth = 0.2),
            legend.position = "right",
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 10)
        )
}

# =====================================================================
# SENSITIVITY ANALYSIS MODULE
# =====================================================================

#' Perform MR-PRESSO Analysis
#' 
#' @param harmonized_data Harmonized exposure-outcome dataset
#' @return MR-PRESSO results
perform_mr_presso_analysis <- function(harmonized_data) {
    presso_input <- data.frame(
        beta.exposure = harmonized_data$beta.exposure,
        beta.outcome = harmonized_data$beta.outcome,
        se.exposure = harmonized_data$se.exposure,
        se.outcome = harmonized_data$se.outcome
    )
    
    mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = presso_input,
        NbDistribution = 1000,
        SignifThreshold = 0.05
    )
}

#' Perform Outlier Detection and Correction
#' 
#' @param harmonized_data Harmonized exposure-outcome dataset
#' @param outlier_threshold Residual quantile threshold (default = 0.95)
#' @return Filtered dataset with outliers removed
detect_and_remove_outliers <- function(harmonized_data, outlier_threshold = 0.95) {
    egger_model <- lm(
        beta.outcome ~ beta.exposure, 
        weights = 1/se.outcome^2, 
        data = harmonized_data
    )
    
    harmonized_data$residuals <- abs(resid(egger_model))
    threshold <- quantile(harmonized_data$residuals, outlier_threshold)
    outlier_snps <- harmonized_data$SNP[harmonized_data$residuals > threshold]
    
    filtered_data <- harmonized_data %>%
        filter(!SNP %in% outlier_snps)
    
    return(filtered_data)
}

# =====================================================================
# ANALYSIS EXECUTION MODULE
# =====================================================================

# --------------------------
# Data Preparation
# --------------------------
gtex_data <- process_gtex_eqtl_data("D:/12109/Documents/GTEx Portal4.csv")
gwas_data <- process_gwas_vcf_data("C:/Users/12109/Desktop/ukb-e-I21_CSA.vcf.gz")

# Format GTEx data for MR
mfap4_data <- format_data(
    gtex_data,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    gene_col = "gene_symbol",
    id_col = "id.exposure"
) %>%
    slice(-c(56, 58))  # Remove specific problematic rows

# --------------------------
# Analysis 1: Exposure 'ukb-a-504' -> Outcome 'ieu-a-798'
# --------------------------
analysis1 <- execute_mr_workflow(
    exposure_id = "ukb-a-504",
    outcome_id = "ieu-a-798",
    pval_threshold = 4.7e-6
)

# Visualization
f_stat_plot1 <- visualize_f_statistics(analysis1$instruments)
reverse_analysis1 <- perform_reverse_mr("ukb-a-504", "ieu-a-798")
direction_plot1 <- create_direction_comparison_plot(
    analysis1$mr_results,
    reverse_analysis1$mr_results,
    "Body Composition → Type 2 Diabetes"
)

# --------------------------
# Analysis 2: Exposure 'ukb-e-3426_CSA' -> Outcome 'prot-a-1885'
# --------------------------
analysis2 <- execute_mr_workflow(
    exposure_id = "ukb-e-3426_CSA",
    outcome_id = "prot-a-1885",
    pval_threshold = 3e-6
)

# Visualization
f_stat_plot2 <- visualize_f_statistics(analysis2$instruments)
reverse_analysis2 <- perform_reverse_mr("ukb-e-3426_CSA", "prot-a-1885")
direction_plot2 <- create_direction_comparison_plot(
    analysis2$mr_results,
    reverse_analysis2$mr_results,
    "Metabolic Factor → Protein Biomarker"
)

# --------------------------
# Analysis 3: Exposure 'MFAP4' -> Outcome 'Acute Myocardial Infarction'
# --------------------------
analysis3_data <- harmonise_data(
    exposure_dat = mfap4_data,
    outcome_dat = gwas_data,
    action = 2
)

analysis3_results <- mr(
    analysis3_data,
    method_list = c(
        "mr_ivw",
        "mr_weighted_median",
        "mr_wald_ratio",
        "mr_simple_mode_nome",
        "mr_weighted_mode_nome"
    )
) %>%
    generate_odds_ratios()

# Visualization
manhattan_plot <- visualize_qtl_manhattan(mfap4_data, 15)
f_stat_plot3 <- visualize_f_statistics(mfap4_data)

# --------------------------
# Sensitivity Analyses for Analysis 3
# --------------------------
heterogeneity3 <- mr_heterogeneity(analysis3_data)
pleiotropy3 <- mr_pleiotropy_test(analysis3_data)
leaveoneout3 <- mr_leaveoneout(analysis3_data)
presso3 <- perform_mr_presso_analysis(analysis3_data)

# Outlier detection and correction
filtered_data <- detect_and_remove_outliers(analysis3_data, 0.95)
filtered_results <- mr(filtered_data)

# =====================================================================
# SUPPLEMENTARY ANALYSES
# =====================================================================

#' Generate Manhattan Plot for GWAS Data
#' 
#' @param gwas_id GWAS study ID
#' @param pval_threshold P-value threshold
#' @param plot_title Plot title
visualize_gwas_manhattan <- function(gwas_id, pval_threshold = 1e-2, plot_title) {
    gwas_instruments <- extract_instruments(
        outcomes = gwas_id,
        p1 = pval_threshold
    ) %>%
        mutate(chr.exposure = as.numeric(chr.exposure))
    
    manhattan(
        gwas_instruments,
        chr = "chr.exposure",
        bp = "pos.exposure",
        p = "pval.exposure",
        snp = "SNP",
        genomewideline = -log10(5e-8),
        suggestiveline = -log10(1e-5),
        col = c("#1f77b4", "#ff7f0e"),
        ylim = c(2, 12),
        main = plot_title,
        cex = 0.7
    )
}

# Generate supplementary Manhattan plots
visualize_gwas_manhattan("ukb-a-504", plot_title = "Body Composition GWAS")
visualize_gwas_manhattan("ieu-a-798", plot_title = "Type 2 Diabetes GWAS")
visualize_gwas_manhattan("ukb-e-3426_CSA", plot_title = "Metabolic Factor GWAS")
visualize_gwas_manhattan("prot-a-1885", plot_title = "Protein Biomarker GWAS")
visualize_gwas_manhattan("ukb-b-453", plot_title = "Cardiovascular Disease GWAS")

# =====================================================================
# RESULTS EXPORT MODULE
# =====================================================================

#' Export Analysis Results to CSV
#' 
#' @param results_list List of analysis results
#' @param output_dir Output directory
export_results <- function(results_list, output_dir = "results") {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }
    
    # Export MR results
    write.csv(
        results_list$mr_results, 
        file.path(output_dir, "mr_primary_results.csv"),
        row.names = FALSE
    )
    
    # Export sensitivity analyses
    write.csv(
        results_list$heterogeneity, 
        file.path(output_dir, "heterogeneity_results.csv"),
        row.names = FALSE
    )
    
    write.csv(
        results_list$pleiotropy, 
        file.path(output_dir, "pleiotropy_results.csv"),
        row.names = FALSE
    )
    
    # Export MR-PRESSO results
    presso_df <- as.data.frame(results_list$presso_results$`Main MR results`)
    write.csv(
        presso_df, 
        file.path(output_dir, "mr_presso_results.csv"),
        row.names = FALSE
    )
}

# Export all results
export_results(analysis1, "analysis1_results")
export_results(analysis2, "analysis2_results")

# =====================================================================
# SESSION INFORMATION
# =====================================================================
session_info <- sessionInfo()
capture.output(session_info, file = "analysis_session_info.txt")