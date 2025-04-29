#!/usr/bin/env Rscript

#getwd()

suppressPackageStartupMessages({
  library(Biostrings)
  library(ggplot2)
  library(pwalign)
  library(gridExtra)
})

# ---------- user paths ----------
seq_dir  <- "fasta"
ref_file <- "rCRS.fasta"
out_pdf  <- "sanger_damage_report.pdf"
boot_iter <- 10000
term_len  <- 10
# --------------------------------

clean_dna <- function(x) {
  DNAString(gsub("[^ACGTacgt]", "", toupper(as.character(x))))
}
# helper to extract ends of an alignment and record CT/GA events
score_ends <- function(aln, term = term_len) {
  ref <- as.character(pattern(aln))
  qry <- as.character(subject(aln))
  ref_vec <- strsplit(ref, "")[[1]]
  qry_vec <- strsplit(qry, "")[[1]]
  # keep only positions that are not gaps in reference
  keep <- ref_vec != "-"
  ref_vec <- ref_vec[keep]
  qry_vec <- qry_vec[keep]
  n <- length(ref_vec)
  left  <- 1:term
  right <- (n - term + 1):n
  # CT at 5' means ref=C, qry=T at left end
  ct_5  <- sum(ref_vec[left] == "C" & qry_vec[left] == "T")
  ga_5  <- sum(ref_vec[left] == "G" & qry_vec[left] == "A")
  ct_3  <- sum(ref_vec[right] == "C" & qry_vec[right] == "T")
  ga_3  <- sum(ref_vec[right] == "G" & qry_vec[right] == "A")
  len_5 <- sum(ref_vec[left] %in% c("A", "C", "G", "T"))
  len_3 <- sum(ref_vec[right] %in% c("A", "C", "G", "T"))
  data.frame(ct_5, ga_5, ct_3, ga_3, len_5, len_3)
}
ref <- clean_dna(readDNAStringSet(ref_file)[[1]])
message("Reading Sanger sequences …")
# define substitution matrix for nucleotide alignment
subsMat <- nucleotideSubstitutionMatrix(match = 1,
                                        mismatch = -1,
                                        baseOnly = FALSE)
files <- list.files(seq_dir, pattern = "\\.fasta$", full.names = TRUE)
# compute per-file damage stats and total letters
scores <- lapply(files, function(f) {
  seq <- readDNAStringSet(f)[[1]]
  aln <- pwalign::pairwiseAlignment(
    ref,
    seq,
    substitutionMatrix = subsMat,
    gapOpening = -5,
    gapExtension = -2,
    type = "global-local"
  )
  # extract aligned sequences, drop ref gaps
  r <- as.character(pattern(aln))
  q <- as.character(subject(aln))
  rv <- strsplit(r, "")[[1]]
  qv <- strsplit(q, "")[[1]]
  keep <- rv != "-"
  rv <- rv[keep]
  qv <- qv[keep]
  n <- length(rv)
  # define terminal windows
  left_i  <- seq_len(term_len)
  right_i <- seq.int(n - term_len + 1, n)
  # count C->T and G->A
  ct5 <- rv[left_i] == "C" & qv[left_i] == "T"
  ga5 <- rv[left_i] == "G" & qv[left_i] == "A"
  ct3 <- rv[right_i] == "C" & qv[right_i] == "T"
  ga3 <- rv[right_i] == "G" & qv[right_i] == "A"
  # valid reference positions
  len5 <- sum(rv[left_i] %in% c("A", "C", "G", "T"))
  len3 <- sum(rv[right_i] %in% c("A", "C", "G", "T"))
  # total canonical bases in sequence
  total_bases <- sum(qv %in% c("A", "C", "G", "T"))
  data.frame(
    file = basename(f),
    ct_5 = sum(ct5),
    ga_5 = sum(ga5),
    len_5 = len5,
    ct_3 = sum(ct3),
    ga_3 = sum(ga3),
    len_3 = len3,
    total_bases = total_bases,
    stringsAsFactors = FALSE
  )
})
scores <- do.call(rbind, scores)
scores$dmg_5 <- (scores$ct_5 + scores$ga_5) / scores$len_5
scores$dmg_3 <- (scores$ct_3 + scores$ga_3) / scores$len_3
scores$p5 <- pbinom(scores$ct_5 + scores$ga_5,
                    scores$len_5, 0.01, lower.tail = FALSE)
scores$p3 <- pbinom(scores$ct_3 + scores$ga_3,
                    scores$len_3, 0.01, lower.tail = FALSE)
scores$sample <- sub("-(F|R)\\.fasta$", "", scores$file)
scores$strand <- sub(".*-([FR])\\.fasta$", "\\1", scores$file)
samples <- unique(scores$sample)
summary_list <- lapply(samples, function(s) {
  df <- scores[scores$sample == s,]
  make_row <- function(sdf, st) {
    if (st %in% sdf$strand) {
      r <- sdf[sdf$strand == st,]
      data.frame(sample = s, strand = st,
                 total_bases = r$total_bases,
                 dmg_5 = round(r$dmg_5,3),
                 dmg_3 = round(r$dmg_3,3),
                 p5 = round(r$p5,4), p3 = round(r$p3,4),
                 stringsAsFactors = FALSE)
    } else {
      data.frame(sample = s, strand = st,
                 total_bases = NA, dmg_5 = NA, dmg_3 = NA,
                 p5 = NA, p3 = NA, stringsAsFactors = FALSE)
    }
  }
  row_F <- make_row(df, "F")
  row_R <- make_row(df, "R")
  if (all(c("F","R") %in% df$strand)) {
    sum_ct5 <- sum(df$ct_5); sum_ga5 <- sum(df$ga_5)
    sum_ct3 <- sum(df$ct_3); sum_ga3 <- sum(df$ga_3)
    sum_len5 <- sum(df$len_5); sum_len3 <- sum(df$len_3)
    sum_letters <- sum(df$total_bases)
    row_FR <- data.frame(sample = s, strand = "FR",
                          total_bases = sum_letters,
                          dmg_5 = round((sum_ct5 + sum_ga5)/sum_len5, 3),
                          dmg_3 = round((sum_ct3 + sum_ga3)/sum_len3, 3),
                          p5 = round(pbinom(sum_ct5 + sum_ga5, sum_len5, 0.01,
                                            lower.tail = FALSE), 4),
                          p3 = round(pbinom(sum_ct3 + sum_ga3, sum_len3, 0.01,
                                            lower.tail = FALSE), 4))
  } else {
    row_FR <- data.frame(sample = s, strand = "FR",
                 total_bases = NA, dmg_5 = NA, dmg_3 = NA,
                 p5 = NA, p3 = NA, stringsAsFactors = FALSE)
  }
  rbind(row_F, row_R, row_FR)
})
summary_df <- do.call(rbind, summary_list)

qc_dmg_thresh <- 0.02
qc_p_thresh <- 0.05

summary_df$QC_auth <- ""
idx <- with(summary_df,
            strand == "FR" & !is.na(dmg_5) & !is.na(dmg_3))
summary_df$QC_auth[idx] <-
  ifelse(summary_df$dmg_5[idx] > qc_dmg_thresh &
           summary_df$dmg_3[idx] > qc_dmg_thresh &
           summary_df$p5[idx]   < qc_p_thresh &
           summary_df$p3[idx]   < qc_p_thresh,
         "OK", "FAIL")

obs <- within(scores, {
  dmg_5 <- (ct_5 + ga_5) / len_5
  dmg_3 <- (ct_3 + ga_3) / len_3
})
obs_overall <- with(obs, c(n        = nrow(obs),
                           dmg_5    = mean(dmg_5, na.rm = TRUE),
                           dmg_3    = mean(dmg_3, na.rm = TRUE)))

# bootstrap sampling distribution
set.seed(123)
boots <- replicate(boot_iter, {
  idx <- sample(seq_len(nrow(obs)), replace = TRUE)
  c(mean(obs$dmg_5[idx], na.rm = TRUE), mean(obs$dmg_3[idx], na.rm = TRUE))
})
boots <- t(boots)
colnames(boots) <- c("boot_dmg_5", "boot_dmg_3")

# one-tailed p-value: how often is boot damage ≤ expected modern error (0.01)?
p_5 <- mean(boots[,"boot_dmg_5"] <= 0.01)
p_3 <- mean(boots[,"boot_dmg_3"] <= 0.01)

# PDF report
pdf(out_pdf, width = 8, height = 6)
# per-sequence summary table: split into pages with 20 samples each
chunk_size <- 20
n_chunks <- ceiling(nrow(summary_df) / chunk_size)
for (i in seq_len(n_chunks)) {
  start <- (i - 1) * chunk_size + 1
  end <- min(i * chunk_size, nrow(summary_df))
  df_chunk <- summary_df[start:end, ]
  if (i > 1) grid::grid.newpage()
  grid.table(df_chunk)
}
# new page for damage plots
grid::grid.newpage()
par(mfrow = c(1,2), mar = c(4,4,2,1))
# misincorporation barplot
barplot(height = c(obs_overall["dmg_5"], obs_overall["dmg_3"]),
        names.arg = c("5' C->T + G->A", "3' C->T + G->A"),
        ylim = c(0, max(0.05, obs_overall["dmg_5"], obs_overall["dmg_3"])*1.2),
        ylab = "Proportion of transitions",
        col = c("steelblue3","tomato"))
abline(h = 0.01, lty = 2)  # modern background threshold
mtext(sprintf("n = %d sequences", obs_overall["n"]), side = 3, line = 0.3, cex = 0.8)

# bootstrap histogram
df <- data.frame(boot = boots[,"boot_dmg_5"])
  ggplot(df, aes(boot)) +
    geom_histogram(binwidth = 0.002) +
    geom_vline(xintercept = obs_overall["dmg_5"], colour = "red") +
    geom_vline(xintercept = 0.01, lty = 2) +
    labs(x = "Bootstrapped mean 5' C->T + G->A", y = "Count") +
    theme_minimal()
  # smile plots: per-sample damage profile across terminal positions
  grid::grid.newpage()
  for (s in samples) {
    # initialize counts and valid base trackers
    left_counts <- rep(0, term_len); left_valid <- rep(0, term_len)
    right_counts <- rep(0, term_len); right_valid <- rep(0, term_len)
    # aggregate over strands F and R
    for (st in c("F", "R")) {
      fpath <- file.path(seq_dir, paste0(s, "-", st, ".fasta"))
      if (file.exists(fpath)) {
        seq <- clean_dna(readDNAStringSet(fpath)[[1]])
        aln <- pwalign::pairwiseAlignment(ref, seq,
                         substitutionMatrix = subsMat,
                         gapOpening = -5, gapExtension = -2,
                         type = "global-local")
        rchar <- as.character(pattern(aln)); qchar <- as.character(subject(aln))
        rvec <- strsplit(rchar, "")[[1]]; qvec <- strsplit(qchar, "")[[1]]
        keep_idx <- rvec != "-"; rvec <- rvec[keep_idx]; qvec <- qvec[keep_idx]
        npos <- length(rvec)
        left_idx <- seq_len(term_len)
        right_idx <- seq.int(npos - term_len + 1, npos)
        for (i in seq_len(term_len)) {
          # 5' end
          if (rvec[left_idx[i]] %in% c("A","C","G","T")) {
            left_valid[i] <- left_valid[i] + 1
            if ((rvec[left_idx[i]] == "C" & qvec[left_idx[i]] == "T") |
                (rvec[left_idx[i]] == "G" & qvec[left_idx[i]] == "A")) {
              left_counts[i] <- left_counts[i] + 1
            }
          }
          # 3' end
          if (rvec[right_idx[i]] %in% c("A","C","G","T")) {
            right_valid[i] <- right_valid[i] + 1
            if ((rvec[right_idx[i]] == "C" & qvec[right_idx[i]] == "T") |
                (rvec[right_idx[i]] == "G" & qvec[right_idx[i]] == "A")) {
              right_counts[i] <- right_counts[i] + 1
            }
          }
        }
      }
    }
    # compute frequencies, handle zero-valid cases
    left_freq <- ifelse(left_valid > 0, left_counts / left_valid, NA)
    right_freq <- ifelse(right_valid > 0, right_counts / right_valid, NA)
    # prepare data frame for plotting
    pos_vals <- c(-(term_len:1), 1:term_len)
    freq_vals <- c(left_freq, right_freq)
    df_smile <- data.frame(pos = pos_vals,
                           freq = freq_vals,
                           end = factor(rep(c("5'","3'"), each = term_len),
                                        levels = c("5'","3'")))
    # plot smile profile
    p_smile <- ggplot(df_smile, aes(x = pos, y = freq, color = end)) +
      geom_line() + geom_point() +
      scale_x_continuous(breaks = pos_vals) +
      labs(title = paste("Damage plot for sample", s),
           x = "Distance from read terminus (positions)",
           y = "Misincorporation rate") +
      theme_minimal()
    print(p_smile)
  }
  # final consensus smile plots for merged sequences in 'final' directory
  grid::grid.newpage()
  final_files <- list.files("final", pattern = "_merged\\.fasta$", full.names = TRUE)
  if (length(final_files) > 0) {
    for (f in final_files) {
      sample <- sub("_merged\\.fasta$", "", basename(f))
      seq <- clean_dna(readDNAStringSet(f)[[1]])
      aln <- pwalign::pairwiseAlignment(ref, seq,
                       substitutionMatrix = subsMat,
                       gapOpening = -5, gapExtension = -2,
                       type = "global-local")
      rchar <- as.character(pattern(aln)); qchar <- as.character(subject(aln))
      rvec <- strsplit(rchar, "")[[1]]; qvec <- strsplit(qchar, "")[[1]]
      keep_idx <- rvec != "-"
      rvec <- rvec[keep_idx]; qvec <- qvec[keep_idx]
      npos <- length(rvec)
      left_idx <- seq_len(term_len)
      right_idx <- seq.int(npos - term_len + 1, npos)
      left_counts <- rep(0, term_len); left_valid <- rep(0, term_len)
      right_counts <- rep(0, term_len); right_valid <- rep(0, term_len)
      for (i in seq_len(term_len)) {
        if (rvec[left_idx[i]] %in% c("A","C","G","T")) {
          left_valid[i] <- left_valid[i] + 1
          if ((rvec[left_idx[i]] == "C" & qvec[left_idx[i]] == "T") |
              (rvec[left_idx[i]] == "G" & qvec[left_idx[i]] == "A")) {
            left_counts[i] <- left_counts[i] + 1
          }
        }
        if (rvec[right_idx[i]] %in% c("A","C","G","T")) {
          right_valid[i] <- right_valid[i] + 1
          if ((rvec[right_idx[i]] == "C" & qvec[right_idx[i]] == "T") |
              (rvec[right_idx[i]] == "G" & qvec[right_idx[i]] == "A")) {
            right_counts[i] <- right_counts[i] + 1
          }
        }
      }
      left_freq <- ifelse(left_valid > 0, left_counts / left_valid, NA)
      right_freq <- ifelse(right_valid > 0, right_counts / right_valid, NA)
      pos_vals <- c(-(term_len:1), 1:term_len)
      freq_vals <- c(left_freq, right_freq)
      df_final <- data.frame(pos = pos_vals,
                             freq = freq_vals,
                             end = factor(rep(c("5'","3'"), each = term_len),
                                          levels = c("5'","3'")))
      p_final <- ggplot(df_final, aes(x = pos, y = freq, color = end)) +
        geom_line() + geom_point() +
        scale_x_continuous(breaks = pos_vals) +
        labs(title = paste("Final damage plot for sample", sample),
             x = "Distance from read terminus (positions)",
             y = "Misincorporation rate") +
        theme_minimal()
      print(p_final)
    }
  }
  invisible(dev.off())

# textual report
cat("\n----- Ancient DNA authentication summary -----\n")
cat("Total sequences analysed:       ", obs_overall["n"], "\n")
cat("Mean 5' C->T + G->A proportion:     ", sprintf("%.3f", obs_overall["dmg_5"]), "\n")
cat("Mean 3' C->T + G->A proportion:     ", sprintf("%.3f", obs_overall["dmg_3"]), "\n")
cat("Bootstrap P(5' damage ≤0.01):  ", sprintf("%.4f", p_5), "\n")
cat("Bootstrap P(3' damage ≤0.01):  ", sprintf("%.4f", p_3), "\n\n")
if (obs_overall["dmg_5"] > 0.02 && p_5 < 0.05 &&
    obs_overall["dmg_3"] > 0.02 && p_3 < 0.05) {
  cat("Interpretation: the damage signal at both termini comfortably exceeds modern\n",
      "background levels. Together with independent extraction and lab controls,\n",
      "these data support authenticity of the mitochondrial fragments.\n")
} else {
  cat("Interpretation: the terminal C->T / G->A excess is modest; further evidence such as\n",
      "replicate extractions, library-wide fragment length and contamination testing\n",
      "with next-generation sequencing is recommended before making strong claims.\n")
}
