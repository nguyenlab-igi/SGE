library(dplyr)
library(tidyr)
library(readr)

clean_Mutation_code <- function(x) {
  # Remove whitespace and make uppercase 
  toupper(trimws(x))
}


# ---- 1. Data Loading and Preprocessing ----

file1 = "ekansh_annotated/annotated_21g_1Libp_S21_L001_VS_LibPlasmid_S1_L001_interface_final_14days_ABSelected.csv"
file2 = "ekansh_annotated/annotated_5g_1Libp_S7_L001_VS_LibTemplate_S2_L001_interface_final_3Day_ABSelected.csv"
loess_span = 0.15

donor1 <-
  read_csv(file1, show_col_types = FALSE)  %>%
    drop_na( AminoAcid_Position, Variant_AminoAcid, Variant_Codon) %>%
    mutate(Mutation_code = clean_Mutation_code(paste(WildType_AminoAcid, AminoAcid_Position, Variant_AminoAcid, "|", Variant_Codon, sep = ""))) %>%
    drop_na(cDNA_Position, Variant_Score, Mutation_Type, Mutation_code)


donor2 <-   read_csv(file2, show_col_types = FALSE) %>%
  drop_na( AminoAcid_Position, Variant_AminoAcid, Variant_Codon) %>%
  mutate(Mutation_code = clean_Mutation_code(paste(WildType_AminoAcid, AminoAcid_Position, Variant_AminoAcid, "|", Variant_Codon, sep = ""))) %>%
  drop_na(cDNA_Position, Variant_Score, Mutation_Type, Mutation_code)

# ---- 2. LOESS Normalization, Join, and Averaging ----
df1 <- donor1
df2 <- donor2

cDNA_all <- sort(unique(c(df1$cDNA_Position, df2$cDNA_Position)))
cDNA_range <- min(cDNA_all):max(cDNA_all)
# --- DONOR 1: LOESS FITTING ---
syn_min1 <- min(df1$Variant_Score[df1$Mutation_Type == "Syn"], na.rm = TRUE)
df1_loess <- df1 %>% filter(Variant_Score >= syn_min1)
loess1 <- loess(Variant_Score ~ cDNA_Position, data = df1_loess, span = loess_span)
loess_pred_donor1 <- data.frame(
  `cDNA_pos` = cDNA_range,
  loess_donor1 = predict(loess1, newdata = data.frame(cDNA_Position = cDNA_range))
)

# --- DONOR 2: LOESS FITTING ---
syn_min2 <- min(df2$Variant_Score[df2$Mutation_Type == "Syn"], na.rm = TRUE)
df2_loess <- df2 %>% filter(Variant_Score >= syn_min2)
loess2 <- loess(Variant_Score ~ cDNA_Position, data = df2_loess, span = loess_span)
loess_pred_donor2 <- data.frame(
  `cDNA_pos` = cDNA_range,
  loess_donor2 = predict(loess2, newdata = data.frame(cDNA_Position = cDNA_range))
)

# Merge both LOESS predictions onto the joined table
shared_variants <- full_join(df1, df2, by = "Mutation_code", suffix = c("_d1", "_d2"))
#shared_variants <- left_join(shared_variants, loess_pred_donor1, by = "cDNA..")
#shared_variants <- left_join(shared_variants, loess_pred_donor2, by = "cDNA..")

# Calculate adjusted values using LOESS predicted from each donor
shared_variants <- shared_variants %>%
  mutate(
    `cDNA_pos` = as.integer(cDNA_Position_d1),
    ClinVar_Classifications = ClinVar_Classification_d1
  ) %>%
  left_join(loess_pred_donor1, by = "cDNA_pos") %>%
  left_join(loess_pred_donor2, by = "cDNA_pos") %>%
  mutate(
    adjusted_d1 = Variant_Score_d1 - loess_donor1,
    adjusted_d2 = Variant_Score_d2 - loess_donor2,
    Mutation.type = Mutation_Type_d1
  )
# Linear median adjustment scaling
shared_variants_syn = shared_variants[which(shared_variants$Mutation.type == "Syn"),]
shared_variants_non = shared_variants[which(shared_variants$Mutation.type == "Non"),]
syn_scale_value = mean(c(median(shared_variants_syn$adjusted_d1), median(shared_variants_syn$adjusted_d2)))
non_scale_value = mean(c(median(shared_variants_non$adjusted_d1), median(shared_variants_non$adjusted_d2)))

syn_median_d1 <- median(shared_variants$adjusted_d1[shared_variants$Mutation.type == "Syn"])
syn_median_d2 <- median(shared_variants$adjusted_d2[shared_variants$Mutation.type == "Syn"])
non_median_d1 <- median(shared_variants$adjusted_d1[shared_variants$Mutation.type == "Non"])
non_median_d2 <- median(shared_variants$adjusted_d2[shared_variants$Mutation.type == "Non"])

# Linear transformation parameters
scale_factor_d1 <- (non_scale_value - syn_scale_value) / (non_median_d1 - syn_median_d1)
shift_factor_d1 <- syn_scale_value - scale_factor_d1 * syn_median_d1
scale_factor_d2 <- (non_scale_value - syn_scale_value) / (non_median_d2 - syn_median_d2)
shift_factor_d2 <- syn_scale_value - scale_factor_d2 * syn_median_d2

# Apply the transformation
shared_variants$ScaledFoldChange_d1 <- shared_variants$adjusted_d1 * scale_factor_d1 + shift_factor_d1
shared_variants$ScaledFoldChange_d2 <- shared_variants$adjusted_d2 * scale_factor_d2 + shift_factor_d2
shared_variants$ScaledFoldChange_mean <- rowMeans(cbind(shared_variants$ScaledFoldChange_d1, shared_variants$ScaledFoldChange_d2))
shared_variants$adj_avg_score <- shared_variants$ScaledFoldChange_mean

# ---- 3. EM Gaussian Mixture Model & Posterior Assignment ----

estimate_extrapolated_cutoffs <- function(df, n_trials, n_points, 
                                          min_val = input$min_val, 
                                          max_val = input$max_val, 
                                          cutoff_mode = "random", 
                                          trim = 0.05, mu_tol = 0.005, 
                                          return_all = FALSE) {
  data_thresholdnon_syn <- df %>%
    filter(Mutation.type %in% c("Non", "Syn")) %>%
    select(adj_avg_score)
  
  # === Mode: Direct EM fit to original Non/Syn data ===
  
    mixmdl <- tryCatch(
      mixtools::normalmixEM(data_thresholdnon_syn$adj_avg_score, k = 2, maxit = 1000, epsilon = 1e-8, verb = FALSE),
      error = function(e) NULL
    )
    
    if (is.null(mixmdl)) {
      return(list(
        cutoffs = data.frame(Est_Cutoff_Low = NA, Est_Cutoff_High = NA),
        random_points = data.frame()
      ))
    }
    
    # Sort components by mean (so lower mean = nonfunc if that's convention)
    ordering <- order(mixmdl$mu)
    mu1 <- mixmdl$mu[ordering[1]]
    mu2 <- mixmdl$mu[ordering[2]]
    sigma1 <- mixmdl$sigma[ordering[1]]
    sigma2 <- mixmdl$sigma[ordering[2]]
    lambda1 <- mixmdl$lambda[ordering[1]]
    lambda2 <- mixmdl$lambda[ordering[2]]
    
    # Posterior for "non-functional"
    ## probablitliy that any given data point comes from component1
    posterior1_fun <- function(x) {
      num <- lambda1 * dnorm(x, mean = mu1, sd = sigma1)
      denom <- num + lambda2 * dnorm(x, mean = mu2, sd = sigma2)
      return(num / denom)
    }
    
    # --- Find cutoffs for intermediate boundaries ---
    score_range <- range(df$adj_avg_score, na.rm = TRUE)
    # Lower: posterior = 0.99 (lower boundary of intermediate)
    lower_cut <- tryCatch(
      uniroot(function(x) posterior1_fun(x) - 0.99, interval = score_range, extendInt = "yes")$root,
      error = function(e) NA_real_
    )
    # Upper: posterior = 0.01 (upper boundary of intermediate)
    upper_cut <- tryCatch(
      uniroot(function(x) posterior1_fun(x) - 0.01, interval = score_range, extendInt = "yes")$root,
      error = function(e) NA_real_
    )
    
    # Apply posterior calculation to all data for output
    df_full <- shared_variants %>%
      mutate(
        posterior_nonfunc = posterior1_fun(adj_avg_score),
        label = case_when(
          posterior_nonfunc > 0.99 ~ "non-functional",
          posterior_nonfunc < 0.01 ~ "functional",
          TRUE ~ "intermediate"
        )
      )
    
    out <- data.frame(
      Est_Cutoff_Low = lower_cut,
      Est_Cutoff_High = upper_cut
    )
    
    return(list(
      cutoffs = out,
      random_points = df_full,
      model = mixmdl  
    ))
  
  
  
  
  
}

df <- shared_variants %>% filter(Mutation.type %in% c("Non", "Syn"))
num_trials = 1000
n_points = 25
min_val = -2
max_val = -0.5
output = estimate_extrapolated_cutoffs(shared_variants,n_trials = num_trials,n_points = num_random,min_val = min_val,max_val = max_val,cutoff_mode = "direct",  return_all = TRUE)
cutoff_low = output$cutoffs$Est_Cutoff_Low
cutoff_high = output$cutoffs$Est_Cutoff_High



full = shared_variants

est_low = output$cutoffs$Est_Cutoff_Low
est_high = output$cutoffs$Est_Cutoff_High
mixmdl = output$model
if (!is.na(est_low) && !is.na(est_high)) {
  full$auto_classification <- dplyr::case_when(
    full$adj_avg_score < est_low  ~ "Non-functional",
    full$adj_avg_score > est_high ~ "Functional",
    TRUE                          ~ "Intermediate"
  )
} else {
  full$auto_classification <- NA_character_
}
if (!is.null(mixmdl)) {
  ord <- order(mixmdl$mu)
  mu1 <- mixmdl$mu[ord[1]];  mu2 <- mixmdl$mu[ord[2]]
  s1  <- mixmdl$sigma[ord[1]]; s2 <- mixmdl$sigma[ord[2]]
  l1  <- mixmdl$lambda[ord[1]]; l2 <- mixmdl$lambda[ord[2]]
  
  posterior1_fun <- function(x) {
    num <- l1 * dnorm(x, mean = mu1, sd = s1)
    den <- num + l2 * dnorm(x, mean = mu2, sd = s2)
    num / den
  }
  full$posterior1 <- posterior1_fun(full$adj_avg_score)
} else {
  full$posterior1 <- NA_real_
}
needed_nums <- c(
  "loess_donor1","loess_donor2","loess_donor3",
  "adjusted_d1","adjusted_d2","adjusted_d3",
  "ScaledFoldChange_d1","ScaledFoldChange_d2","ScaledFoldChange_d3",
  "ScaledFoldChange_mean",
  "Variant_Score_d1","Variant_Score_d2","Variant_Score_d3",
  "cDNA_pos"
)
for (nm in needed_nums) if (!nm %in% names(full)) full[[nm]] <- NA_real_
if ("Variant_Score" %in% names(full) && all(is.na(full$Variant_Score_d3))) {
  full$Variant_Score_d3 <- full$Variant_Score
}

needed_chars <- c(
  "cDNA_Position_d1","cDNA_Position_d2","cDNA_Position",
  "AminoAcid_Position_d1","AminoAcid_Position_d2","AminoAcid_Position",
  "ClinVar_Classification_d1","ClinVar_Classification_d2","ClinVar_Classification",
  "Mutation.type","Mutation_code"
)
for (nm in needed_chars) if (!nm %in% names(full)) full[[nm]] <- NA
full <- full %>%
  mutate(
    cDNA_Position = coalesce(
      cDNA_Position_d1,
      cDNA_Position_d2,
      cDNA_Position,
      `cDNA_pos`
    )
  )
full = full %>%
  mutate(across(where(is.list), ~as.character(.))) %>%
  select(-adj_avg_score) %>%
  rename(
    LOESS_Adjustment_d1      = loess_donor1,
    LOESS_Adjustment_d2      = loess_donor2,
    LOESS_Adjustment_d3      = loess_donor3,
    LOESS_Adjusted_Score_d1  = adjusted_d1,
    LOESS_Adjusted_Score_d2  = adjusted_d2,
    LOESS_Adjusted_Score_d3  = adjusted_d3,
    Median_Adjusted_Score_d1 = ScaledFoldChange_d1,
    Median_Adjusted_Score_d2 = ScaledFoldChange_d2,
    Median_Adjusted_Score_d3 = ScaledFoldChange_d3,
    New_Variant_Score        = ScaledFoldChange_mean,
    Variant_Classification   = auto_classification,
    `Non-functional_Probability` = posterior1
  )

renames <- list(
  ScaledFoldChange_d1   = "Median_Adjusted_Score_d1",
  ScaledFoldChange_d2   = "Median_Adjusted_Score_d2",
  ScaledFoldChange_d3   = "Median_Adjusted_Score_d3",
  ScaledFoldChange_mean = "New_Variant_Score",
  auto_classification   = "Variant_Classification",
  posterior1            = "Non-functional_Probability"
)
for (old in names(renames)) {
  new <- renames[[old]]
  if (old %in% names(full) && !(new %in% names(full))) {
    names(full)[names(full) == old] <- new
  }
}



n <- nrow(full)
as_chr_vec <- function(x, n) {
  if (is.null(x)) return(rep(NA_character_, n))
  if (is.list(x)) {
    v <- vapply(x, function(el) if (length(el)) as.character(el[[1]]) else NA_character_, character(1))
    return(v)
  }
  as.character(x)
}
as_num_vec <- function(x, n) {
  if (is.null(x)) return(rep(NA_real_, n))
  if (is.list(x)) {
    v <- vapply(x, function(el) if (length(el)) suppressWarnings(as.numeric(el[[1]])) else NA_real_, numeric(1))
    return(v)
  }
  suppressWarnings(as.numeric(x))
}
getv_chr <- function(nm) if (nm %in% names(full)) as_chr_vec(full[[nm]], n) else rep(NA_character_, n)
getv_num <- function(nm) if (nm %in% names(full)) as_num_vec(full[[nm]], n) else rep(NA_real_, n)
cDNA_Position <- dplyr::coalesce(
  suppressWarnings(as.integer(getv_chr("cDNA_Position_d1"))),
  suppressWarnings(as.integer(getv_chr("cDNA_Position_d2"))),
  suppressWarnings(as.integer(getv_chr("cDNA_Position"))),
  suppressWarnings(as.integer(getv_chr("cDNA_pos")))
)

AminoAcid_Position <- dplyr::coalesce(
  suppressWarnings(as.integer(getv_chr("AminoAcid_Position_d1"))),
  suppressWarnings(as.integer(getv_chr("AminoAcid_Position_d2"))),
  suppressWarnings(as.integer(getv_chr("AminoAcid_Position")))
)

ClinVar_Classifications <- dplyr::coalesce(
  getv_chr("ClinVar_Classification_d1"),
  getv_chr("ClinVar_Classification_d2"),
  getv_chr("ClinVar_Classification")
)

Variant_Codon <- dplyr::coalesce(
  getv_chr("Variant_Codon_d2"),
  getv_chr("Variant_Codon_d1"),
  getv_chr("Variant_Codon")
)

WT_AA <- dplyr::coalesce(
  getv_chr("WildType_AminoAcid_d2"),
  getv_chr("WildType_AminoAcid_d1"),
  getv_chr("WildType_AminoAcid")
)

mas1 <- getv_num("Median_Adjusted_Score_d1")
mas2 <- getv_num("Median_Adjusted_Score_d2")
mas3 <- getv_num("Median_Adjusted_Score_d3")
Variant_SD <- apply(cbind(mas1, mas2, mas3), 1, function(v) if (all(is.na(v))) NA_real_ else sd(v, na.rm = TRUE))

.genetic_code <- c(
  "TTT"="F","TTC"="F","TTA"="L","TTG"="L",
  "TCT"="S","TCC"="S","TCA"="S","TCG"="S",
  "TAT"="Y","TAC"="Y","TAA"="*","TAG"="*",
  "TGT"="C","TGC"="C","TGA"="*","TGG"="W",
  "CTT"="L","CTC"="L","CTA"="L","CTG"="L",
  "CCT"="P","CCC"="P","CCA"="P","CCG"="P",
  "CAT"="H","CAC"="H","CAA"="Q","CAG"="Q",
  "CGT"="R","CGC"="R","CGA"="R","CGG"="R",
  "ATT"="I","ATC"="I","ATA"="I","ATG"="M",
  "ACT"="T","ACC"="T","ACA"="T","ACG"="T",
  "AAT"="N","AAC"="N","AAA"="K","AAG"="K",
  "AGT"="S","AGC"="S","AGA"="R","AGG"="R",
  "GTT"="V","GTC"="V","GTA"="V","GTG"="V",
  "GCT"="A","GCC"="A","GCA"="A","GCG"="A",
  "GAT"="D","GAC"="D","GAA"="E","GAG"="E",
  "GGT"="G","GGC"="G","GGA"="G","GGG"="G"
)

codon_to_aa <- function(codon) {
  codon <- toupper(trimws(codon))
  if (is.na(codon) || nchar(codon) != 3 || grepl("[^ACGT]", codon)) return(NA_character_)
  aa <- .genetic_code[[codon]]
  if (is.null(aa)) NA_character_ else aa
}

infer_wt_base <- function(variant_codon, codon_pos, wt_aa) {
  if (is.na(variant_codon) || is.na(codon_pos) || is.na(wt_aa)) return(NA_character_)
  variant_codon <- toupper(variant_codon)
  if (nchar(variant_codon) != 3 || !(codon_pos %in% 1:3)) return(NA_character_)
  wt_aa <- toupper(wt_aa)
  for (b in c("A","C","G","T")) {
    cand <- variant_codon
    substr(cand, codon_pos, codon_pos) <- b
    aa <- codon_to_aa(cand)
    if (!is.na(aa) && aa == wt_aa) return(b)
  }
  NA_character_
}

aa1_to3 <- c(
  A="Ala", R="Arg", N="Asn", D="Asp", C="Cys",
  Q="Gln", E="Glu", G="Gly", H="His", I="Ile",
  L="Leu", K="Lys", M="Met", F="Phe", P="Pro",
  S="Ser", T="Thr", W="Trp", Y="Tyr", V="Val",
  "*"="Ter", X="Ter"
)

to_aa3 <- function(aa1) {
  aa1 <- toupper(aa1)
  val <- aa1_to3[[aa1]]
  if (is.null(val)) aa1 else val
}

pc_1L_to_3L <- function(pc) {
  if (is.na(pc) || pc == "") return(NA_character_)
  s <- sub("^p\\.", "", pc)
  m <- stringr::str_match(s, "^([A-Z\\*])([0-9]+)([A-Z\\*])$")
  if (is.na(m[1,1])) return(pc)  
  
  from1 <- m[2]; pos <- m[3]; to1 <- m[4]
  from3 <- to_aa3(from1); to3 <- to_aa3(to1)
  
  if (!is.na(from1) && !is.na(to1) && from1 == to1) {
    # synonymous: p.His242=
    paste0("p.", from3, pos, "=")
  } else {
    # missense/nonsense: p.His242Tyr, p.Gln100Ter, etc.
    paste0("p.", from3, pos, to3)
  }
}

# Protein change from Mutation_code left side -> "p.Ser225Arg"
mc_left     <- sub("\\|.*$", "", getv_chr("Mutation_code"))
prot_change <- ifelse(is.na(mc_left) | mc_left == "", NA_character_, paste0("p.", mc_left))
prot_change_3L <- vapply(prot_change, pc_1L_to_3L, character(1))

min_pos <- suppressWarnings(min(cDNA_Position, na.rm = TRUE))
codon_pos <- if (is.finite(min_pos)) ((as.integer(cDNA_Position) - min_pos) %% 3L) + 1L else rep(NA_integer_, length(cDNA_Position))
alt_base <- ifelse(!is.na(Variant_Codon) & !is.na(codon_pos) & nchar(Variant_Codon) == 3,
                   substring(Variant_Codon, codon_pos, codon_pos),
                   NA_character_)
wt_base  <- mapply(infer_wt_base, Variant_Codon, codon_pos, WT_AA, USE.NAMES = FALSE)
#cdna_offset <- input$cdna_offset  # adjust if gene changes
cdna_offset = -38
cdot_pos <- suppressWarnings(as.integer(cDNA_Position) + cdna_offset)
ClinVar_ID <- ifelse(
  is.na(cdot_pos) | is.na(wt_base) | is.na(alt_base) | is.na(prot_change_3L) |
    nchar(wt_base) != 1 | nchar(alt_base) != 1,
  NA_character_,
  paste0("c.", cdot_pos, wt_base, ">", alt_base, " (", prot_change_3L, ")")
)
vs1 <- getv_num("Variant_Score_d1")
vs2 <- getv_num("Variant_Score_d2")
vs3 <- dplyr::coalesce(getv_num("Variant_Score_d3"), getv_num("Variant_Score"))
output_clean = tibble::tibble(
  cDNA_Position                = as.integer(cDNA_Position),
  Mutation_Type                = getv_chr("Mutation.type"),
  AminoAcid_Position           = as.integer(AminoAcid_Position),
  ClinVar_Classification       = ClinVar_Classifications,
  Mutation_Code                = getv_chr("Mutation_code"),
  Variant_Score_d1             = vs1,
  Variant_Score_d2             = vs2,
  Variant_Score_d3             = vs3,
  Variant_SD                   = Variant_SD,
  New_Variant_Score            = getv_num("New_Variant_Score"),
  Variant_Classification       = getv_chr("Variant_Classification"),
  `Non-functional_Probability` = getv_num("Non-functional_Probability"),
  ClinVar_ID                   = ClinVar_ID
)

## looking at the sd analysis

### sd only
variant_sd <- apply(
  cbind(df$ScaledFoldChange_d1, df$ScaledFoldChange_d2, df$ScaledFoldChange_d3),
  1,
  function(v) { v <- as.numeric(v); if (sum(!is.na(v)) < 2) NA_real_ else stats::sd(v, na.rm = TRUE) }
)

both_avail <- !is.na(cutoff_low) && !is.na(cutoff_high)
int_min    <- if (both_avail) min(cutoff_low, cutoff_high) else NA_real_
int_max    <- if (both_avail) max(cutoff_low, cutoff_high) else NA_real_
eps=0
sd_data = df %>%
  mutate(
    Variant_SD = variant_sd,
    score_mean = ScaledFoldChange_mean,
    Functionality = dplyr::case_when(
      !is.na(score_mean) & !is.na(cutoff_low)  & score_mean < cutoff_low  ~ "Non-functional",
      !is.na(score_mean) & !is.na(cutoff_high) & score_mean > cutoff_high ~ "Functional",
      TRUE ~ "Intermediate"
    )
  ) %>%
  filter(!is.na(Variant_SD), !is.na(score_mean)) %>%
  mutate(
    bar_min    = score_mean - Variant_SD,
    bar_max    = score_mean + Variant_SD,
    spans_both = both_avail &
      (bar_min <= (int_min + eps)) &
      (bar_max >= (int_max - eps)),
    Overlap    = ifelse(spans_both, "Spans both thresholds", "Does not span both thresholds"),
    .hover = paste0(
      "cDNA Pos: ", cDNA..,
      "<br>Mutation Code: ", Mutation_code,
      "<br>New Variant Score: ", round(score_mean, 3),
      "<br>SD: ", round(Variant_SD, 3),
      "<br>Donor1: ", round(ScaledFoldChange_d1, 3),
      "<br>Donor2: ", round(ScaledFoldChange_d2, 3),
      "<br>Mutation Type: ", Mutation.type,
      "<br>ClinVar: ", ClinVar_Classifications,
      "<br>Classification: ", Functionality
    )
  ) %>%
  arrange(Overlap) # non-overlappers first; overlappers on top

