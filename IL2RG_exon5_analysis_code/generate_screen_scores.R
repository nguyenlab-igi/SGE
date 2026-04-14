library(dplyr)
library(tidyverse)
library(openxlsx)
#______________________________________________________________________
#Donor 1 gDNA 3 Day
annotated_5g_1Libp <- read.delim("diff_exp_results_annotated_lib_with_remaining_variant/annotated_5g_1Libp_S7_L001_VS_LibTemplate_S2_L001.csv")%>%drop_na(`WTaa`)
Mutcode<- paste(annotated_5g_1Libp$`WTaa`, annotated_5g_1Libp$`pro..`, annotated_5g_1Libp$mutaminoacid, "|", annotated_5g_1Libp$`mutcodon`, sep = "")
annotated_5g_1Libp$Mutcode <- Mutcode
#Donor 2 gDNA 3Day
annotated_13g_2Libp <- read.delim("diff_exp_results_annotated_lib_with_remaining_variant/annotated_13g_2Libp_S15_L001_VS_LibTemplate_S2_L001.csv")%>%drop_na(`WTaa`)
Mutcode<- paste(annotated_13g_2Libp$`WTaa`, annotated_13g_2Libp $`pro..`, annotated_13g_2Libp $mutaminoacid,"|", annotated_13g_2Libp$`mutcodon`,sep = "")
annotated_13g_2Libp$Mutcode <- Mutcode
#creating a combined column
data_5g_13g_combined <- data.frame(
  LogFC_5g = annotated_5g_1Libp$fold_change,
  LogFC_13g = annotated_13g_2Libp$fold_change
)

# Calculate additional columns
data_5g_13g_combined <- data_5g_13g_combined %>%
  mutate(
    avgFC = (LogFC_5g + LogFC_13g) / 2,
    varFC = (LogFC_5g - avgFC)^2,
    sdFC = sqrt(varFC),
    `cDNA..` = annotated_5g_1Libp$`cDNA..`,
    `pro..` = annotated_5g_1Libp$`pro..`,
    Mutation.type = annotated_5g_1Libp$Mutation.type,
    ClinVar = annotated_5g_1Libp$ClinVar, 
    Mutcode = annotated_5g_1Libp$Mutcode, 
    ClinVar = coalesce(ClinVar, "Unclassified"))

data_5g_Loess <- data_5g_13g_combined%>%filter(LogFC_5g >=	
                                                 -1.8081)
loess_15_5g <- loess(LogFC_5g ~ `cDNA..`, data = data_5g_Loess, span = 0.15)
data_5g_Loess$smoothed15 <- predict(loess_15_5g) 
data_13g_Loess <- data_5g_13g_combined%>%filter(LogFC_13g >=	
                                                  -1.7643)
## fixed in july 2025
loess_15_13g <- loess(LogFC_13g ~ `cDNA..`, data = data_13g_Loess, span = 0.15)
data_13g_Loess$smoothed15 <- predict(loess_15_13g) 

loess_predictions_5g <- data.frame(predict(loess_15_5g , newdata = data.frame('cDNA..' = 633:795)))
colnames(loess_predictions_5g) = "LOESS_5g"
loess_predictions_5g <- loess_predictions_5g%>%mutate('cDNA..' = 633:795)

loess_predictions_13g <- data.frame(predict(loess_15_13g , newdata = data.frame(`cDNA..` = 633:795)))
colnames(loess_predictions_13g) = "LOESS_13g"

loess_predictions_13g <- loess_predictions_13g%>%mutate(`cDNA..` = 633:795)

data_5g_13g_combined<- left_join(data_5g_13g_combined, loess_predictions_5g , by = "cDNA..")
data_5g_13g_combined<- left_join(data_5g_13g_combined, loess_predictions_13g , by = "cDNA..")
data_5g_13g_combined <-data_5g_13g_combined%>%mutate(logFC_5g_adj = (LogFC_5g - LOESS_5g))%>%mutate(logFC_13g_adj = (LogFC_13g - LOESS_13g))
data_5g_13g_combined<-data_5g_13g_combined%>%mutate(logFC_adj_avg = (logFC_13g_adj+logFC_5g_adj)/2)

data_5g_13g_combined_syn = data_5g_13g_combined[which(data_5g_13g_combined$Mutation.type=="Syn"),]
data_5g_13g_combined_non = data_5g_13g_combined[which(data_5g_13g_combined$Mutation.type=="Non"),]
syn_scale_value = mean(c(median(data_5g_13g_combined_syn$logFC_5g_adj),median(data_5g_13g_combined_syn$logFC_13g_adj)))
non_scale_value = mean(c(median(data_5g_13g_combined_non$logFC_5g_adj),median(data_5g_13g_combined_non$logFC_13g_adj)))

syn_median_d1 <- median(data_5g_13g_combined$logFC_5g_adj[data_5g_13g_combined$Mutation.type == "Syn"])
syn_median_d2 <- median(data_5g_13g_combined$logFC_13g_adj[data_5g_13g_combined$Mutation.type == "Syn"])

non_median_d1 <- median(data_5g_13g_combined$logFC_5g_adj[data_5g_13g_combined$Mutation.type == "Non"])
non_median_d2 <- median(data_5g_13g_combined$logFC_13g_adj[data_5g_13g_combined$Mutation.type == "Non"])

# Define the linear transformation parameters
scale_factor_d1 <- (non_scale_value - syn_scale_value) / (non_median_d1 - syn_median_d1)
shift_factor_d1 <- syn_scale_value - scale_factor_d1 * syn_median_d1


scale_factor_d2 <- (non_scale_value - syn_scale_value) / (non_median_d2 - syn_median_d2)
shift_factor_d2 <- syn_scale_value - scale_factor_d2 * syn_median_d2

# Apply the transformation to the fold change values
data_5g_13g_combined$ScaledFoldChange_d1 <- data_5g_13g_combined$logFC_5g_adj * scale_factor_d1 + shift_factor_d1
data_5g_13g_combined$ScaledFoldChange_d2 <- data_5g_13g_combined$logFC_13g_adj * scale_factor_d2 + shift_factor_d2
data_5g_13g_combined$ScaledFoldChange_mean = rowMeans(cbind(data_5g_13g_combined$ScaledFoldChange_d1,data_5g_13g_combined$ScaledFoldChange_d2))


### redo cdna 
annotated_3m_1Libp<-read.delim("diff_exp_results_annotated_lib_with_remaining_variant/annotated_3m_1Libp_S34_L001_VS_LibTemplate_S2_L001.csv")%>%drop_na(`WTaa`)
Mutcode<- paste(annotated_3m_1Libp$`WTaa`, annotated_3m_1Libp$`pro..`, annotated_3m_1Libp$mutaminoacid, "|", annotated_3m_1Libp$`mutcodon`,sep = "")
annotated_3m_1Libp$Mutcode <- Mutcode
annotated_11m_2Libp <- read.delim("diff_exp_results_annotated_lib_with_remaining_variant/annotated_11m_2Libp_S42_L001_VS_LibTemplate_S2_L001.csv")%>%drop_na(`WTaa`)
Mutcode<- paste(annotated_11m_2Libp $`WTaa`, annotated_11m_2Libp $`pro..`, annotated_11m_2Libp $mutaminoacid,"|", annotated_11m_2Libp$`mutcodon`,sep = "")
annotated_11m_2Libp $Mutcode <- Mutcode

data_3m_11m_combined <- data.frame(
  LogFC_3m = annotated_3m_1Libp$fold_change,
  LogFC_11m = annotated_11m_2Libp$fold_change
)
data_3m_11m_combined <- data_3m_11m_combined %>%
  mutate(
    avgFC = (LogFC_3m + LogFC_11m) / 2,
    varFC = (LogFC_3m - avgFC)^2,
    sdFC = sqrt(varFC),
    `cDNA..` = annotated_3m_1Libp$`cDNA..`,
    `pro..` = annotated_3m_1Libp$`pro..`,
    Mutation.type = annotated_3m_1Libp$Mutation.type,
    ClinVar = annotated_3m_1Libp$ClinVar, 
    Mutcode = annotated_3m_1Libp$Mutcode, 
    ClinVar = coalesce(ClinVar, "Unclassified"))
# minimum of the synonymous
data_3m_Loess <- data_3m_11m_combined%>%filter(LogFC_3m >=	
                                                 -2.7342)
loess_15_3m <- loess(LogFC_3m ~ `cDNA..`, data = data_3m_Loess, span = 0.15)
data_3m_Loess$smoothed15 <- predict(loess_15_3m) 


## the new variant becomes missing aftewrds, will set to 0



#########LOESSfit for 11m (Donor2)
data_11m_Loess <- data_3m_11m_combined%>%filter(LogFC_11m >=	
                                                  -5.1448)
## fixed in july 2025
loess_15_11m <- loess(LogFC_11m ~ `cDNA..`, data = data_11m_Loess, span = 0.15)
data_11m_Loess$smoothed15 <- predict(loess_15_11m) 


#adding LOESS predictions to data set 
loess_predictions_3m <- data.frame(predict(loess_15_3m , newdata = data.frame(`cDNA..` = 633:795)))
colnames(loess_predictions_3m) = "LOESS_3m"
loess_predictions_3m <- loess_predictions_3m%>%mutate(`cDNA..` = 633:795)
loess_predictions_3m$LOESS_3m[is.na(loess_predictions_3m$LOESS_3m)] = 0
loess_predictions_11m <- data.frame(predict(loess_15_11m , newdata = data.frame(`cDNA..` = 633:795)))
colnames(loess_predictions_11m) = "LOESS_11m"
loess_predictions_11m <- loess_predictions_11m%>%mutate(`cDNA..` = 633:795)

data_3m_11m_combined<- left_join(data_3m_11m_combined, loess_predictions_3m , by = "cDNA..")
data_3m_11m_combined<- left_join(data_3m_11m_combined, loess_predictions_11m , by = "cDNA..")
data_3m_11m_combined <-data_3m_11m_combined%>%mutate(logFC_3m_adj = (LogFC_3m - LOESS_3m))%>%mutate(logFC_11m_adj = (LogFC_11m - LOESS_11m))
data_3m_11m_combined<-data_3m_11m_combined%>%mutate(logFC_adj_avg = (logFC_11m_adj+logFC_3m_adj)/2)
data_3m_11m_combined_syn = data_3m_11m_combined[which(data_3m_11m_combined$Mutation.type=="Syn"),]
data_3m_11m_combined_non = data_3m_11m_combined[which(data_3m_11m_combined$Mutation.type=="Non"),]
syn_scale_value = mean(c(median(data_3m_11m_combined_syn$logFC_3m_adj),median(data_3m_11m_combined_syn$logFC_11m_adj)))
non_scale_value = mean(c(median(data_3m_11m_combined_non$logFC_3m_adj),median(data_3m_11m_combined_non$logFC_11m_adj)))

syn_median_d1 <- median(data_3m_11m_combined$logFC_3m_adj[data_3m_11m_combined$Mutation.type == "Syn"])
syn_median_d2 <- median(data_3m_11m_combined$logFC_11m_adj[data_3m_11m_combined$Mutation.type == "Syn"])

non_median_d1 <- median(data_3m_11m_combined$logFC_3m_adj[data_3m_11m_combined$Mutation.type == "Non"])
non_median_d2 <- median(data_3m_11m_combined$logFC_11m_adj[data_3m_11m_combined$Mutation.type == "Non"])

# Define the linear transformation parameters
scale_factor_d1 <- (non_scale_value - syn_scale_value) / (non_median_d1 - syn_median_d1)
shift_factor_d1 <- syn_scale_value - scale_factor_d1 * syn_median_d1


scale_factor_d2 <- (non_scale_value - syn_scale_value) / (non_median_d2 - syn_median_d2)
shift_factor_d2 <- syn_scale_value - scale_factor_d2 * syn_median_d2

# Apply the transformation to the fold change values
data_3m_11m_combined$ScaledFoldChange_d1 <- data_3m_11m_combined$logFC_3m_adj * scale_factor_d1 + shift_factor_d1
data_3m_11m_combined$ScaledFoldChange_d2 <- data_3m_11m_combined$logFC_11m_adj * scale_factor_d2 + shift_factor_d2
data_3m_11m_combined$ScaledFoldChange_mean = rowMeans(cbind(data_3m_11m_combined$ScaledFoldChange_d1,data_3m_11m_combined$ScaledFoldChange_d2))

master_table = read.xlsx("IL2RG_EXON5_table_92824.xlsx")

## remove silent mutations until they work...
data_3m_11m_combined = data_3m_11m_combined[-which(is.na(data_3m_11m_combined$pro..)),]
data_5g_13g_combined = data_5g_13g_combined[-which(is.na(data_5g_13g_combined$pro..)),]


old_master_table = read.xlsx("28m_2Linp_S56_VS_LibTemplate_S2_1_coords-lookupLogFCs_by_mut_type-update9-2021-validate_with_last_variant_aug2025.xlsx")
old_master_table$Mutcode=Mutcode<- paste(old_master_table$`WTaa`, old_master_table$`pro.#`, old_master_table$mutaminoacid,"|", old_master_table$`mutcodon`,sep = "")
old_master_table_wt = old_master_table[which(old_master_table$Mutation.type=="WT_"),]
old_master_table_non_wt = old_master_table[-which(old_master_table$Mutation.type=="WT_"),]
#old_master_table_non_wt = old_master_table_non_wt[-which(is.na(old_master_table_non_wt$mutcodon)),]

match_to_master_2_gdna_d1 = c()
match_to_master_2_gdna_d2 = c()
match_to_master_2_gdna_avg = c()

match_to_master_2_cdna_d1 = c()
match_to_master_2_cdna_d2 = c()
match_to_master_2_cdna_avg = c()
match_to_master_annotations = c()
for (m in old_master_table_non_wt$Mutcode) { 
  match_to_master_2_gdna_d1 = c(match_to_master_2_gdna_d1,data_5g_13g_combined$ScaledFoldChange_d1[which(data_5g_13g_combined$Mutcode==m)])
  match_to_master_2_gdna_d2 = c(match_to_master_2_gdna_d2,data_5g_13g_combined$ScaledFoldChange_d2[which(data_5g_13g_combined$Mutcode==m)])
  match_to_master_2_gdna_avg = c(match_to_master_2_gdna_avg, data_5g_13g_combined$ScaledFoldChange_mean[which(data_5g_13g_combined$Mutcode==m)])
  match_to_master_2_cdna_d1 = c(match_to_master_2_cdna_d1,data_3m_11m_combined$ScaledFoldChange_d1[which(data_3m_11m_combined$Mutcode==m)])
  match_to_master_2_cdna_d2 = c(match_to_master_2_cdna_d2,data_3m_11m_combined$ScaledFoldChange_d2[which(data_3m_11m_combined$Mutcode==m)])
  match_to_master_2_cdna_avg = c(match_to_master_2_cdna_avg,data_3m_11m_combined$ScaledFoldChange_mean[which(data_3m_11m_combined$Mutcode==m)])
  
  match_to_master_annotations = c(match_to_master_annotations,old_master_table$ClinVar[which(old_master_table$Mutcode==m)])
}
old_master_table_non_wt$new_score_gdna_d1 = match_to_master_2_gdna_d1
old_master_table_non_wt$new_score_gdna_d2 = match_to_master_2_gdna_d2
old_master_table_non_wt$new_score_gdna_avg = match_to_master_2_gdna_avg
old_master_table_non_wt$new_score_cdna_d1 = match_to_master_2_cdna_d1
old_master_table_non_wt$new_score_cdna_d2 = match_to_master_2_cdna_d2
old_master_table_non_wt$new_score_cdna_avg = match_to_master_2_cdna_avg
old_master_table_non_wt$ClinVar = match_to_master_annotations
old_master_table_wt$new_score_gdna_d1 = rep(NA,nrow(old_master_table_wt))
old_master_table_wt$new_score_gdna_d2 = rep(NA,nrow(old_master_table_wt))
old_master_table_wt$new_score_gdna_avg = rep(NA,nrow(old_master_table_wt))
old_master_table_wt$new_score_cdna_d1 = rep(NA,nrow(old_master_table_wt))
old_master_table_wt$new_score_cdna_d2 = rep(NA,nrow(old_master_table_wt))
old_master_table_wt$new_score_cdna_avg= rep(NA,nrow(old_master_table_wt))
edited_master_table = rbind(old_master_table_non_wt,old_master_table_wt)

edited_master_table = edited_master_table[order(edited_master_table$`cDNA.#`,edited_master_table$Y),]
#write.xlsx(edited_master_table,"master_table_with_new_scores_after_excluding_11624.xlsx")

## let's try the EM scores
#mixmodel with nonsense/ synonomus 
data_5g_13g_nonsyn <- data_5g_13g_combined%>%filter(Mutation.type %in% c("Non", "Syn"))


# Assuming SCID_data3 is a data frame containing the 'logFC' column
#data_thresholdnon_syn <- data_5g_13g_nonsyn$ScaledFoldChange_mean
library(mixtools)
# Fit a Gaussian mixture model with 2 components
# set.seed(70)
# mixmdl_5g_13g_nonsyn <- normalmixEM(LOGfc_5g_13g_nonsyn , k = 2)
# mixmdl_5g_13g_nonsyn$Mutation.type <- data_5g_13g_nonsyn$Mutation.type
# ## number of iterations= 4
# posterior_5g_13g_nonsyn <- mixmdl_5g_13g_nonsyn$posterior
# # Define thresholds for classifying functional SNVs
# thresholds_5g_13g_nonsyn<- c(.99, .01)
# functional_classification_5g_13g_nonsyn <- data.frame(logFC_adj_avg =data_5g_13g_nonsyn$ScaledFoldChange_mean, 
#                                                       posterior = mixmdl_5g_13g_nonsyn$posterior,
#                                                       label = ifelse(
#                                                         posterior_5g_13g_nonsyn[, 1] > thresholds_5g_13g_nonsyn[1], 'non-functional',
#                                                         ifelse(
#                                                           posterior_5g_13g_nonsyn[, 1] < thresholds_5g_13g_nonsyn[2], 'functional', 'intermediate'
#                                                         )
#                                                       )
# )
# data_thresholdnon_syn <- data_5g_13g_nonsyn%>%select(ScaledFoldChange_mean)
# set.seed(72) # Ensure reproducibility
# 
# # Define a function to fit the GMM, classify data, and repeat the process
# run_gmm_classification <- function( num_iterations = 1000) {
#   # Initialize a data frame to store results
#   new_point_labels <- data.frame(
#     ScaledFoldChange_mean = numeric(),
#     posterior.comp.1 = numeric(),
#     posterior.comp.2 = numeric(),
#     label = character(),
#     iteration = integer(),
#     stringsAsFactors = FALSE
#   )
#   for (iteration in 1:num_iterations) {
#     # Add a random data point within the specified range
#     new_point <- data.frame(ScaledFoldChange_mean = runif(1, min = -2, max = -0.5))
#     data <- rbind(data_thresholdnon_syn , new_point)
#     
#     # Fit the Gaussian mixture model with 2 components
#     mixmdl <- normalmixEM(data$ScaledFoldChange_mean, k = 2)
#     posterior <- mixmdl$posterior
#     
#     # Define thresholds for classifying functional SNVs
#     thresholds <- c(0.99, 0.01)
#     
#     # Classify SNVs based on posterior probabilities for the new point
#     new_point_posterior <- posterior[nrow(posterior), ]
#     
#     # Determine the label for the new point
#     new_point_label <- ifelse(
#       new_point_posterior[1] > thresholds[1], 'non-functional',
#       ifelse(new_point_posterior[1] < thresholds[2], 'functional', 'intermediate')
#     )
#     
#     # Store only the new point's classification
#     new_point_labels <- rbind(new_point_labels, data.frame(
#       ScaledFoldChange_mean = new_point$ScaledFoldChange_mean,
#       posterior.comp.1 = new_point_posterior[1],
#       posterior.comp.2 = new_point_posterior[2],
#       label = new_point_label,
#       iteration = iteration
#     ))
#     
#   }
#   return(new_point_labels)
# }
# 
# results <- run_gmm_classification()
# run_gmm_classification_and_extract_ranges <- function() {
#   # Initialize an empty list to store range results
#   range_results_list <- vector("list", 25)
#   
#   # Run the classification 100 times
#   for (i in 1:25) {
#     # Run the classification
#     result <- run_gmm_classification()
#     
#     # Filter the results to get only the rows labeled "intermediate"
#     intermediate_results <- result[result$label == "intermediate", ]
#     
#     # Check if there are any intermediate results
#     if (nrow(intermediate_results) > 0) {
#       # Calculate the range of values in the "intermediate" label column
#       value_range <- range(intermediate_results$ScaledFoldChange_mean, na.rm = TRUE)
#     } else {
#       # If no intermediate results, set the range to NA
#       value_range <- c(NA, NA)
#     }
#     
#     # Store the range in the list
#     range_results_list[[i]] <- data.frame(
#       Run = i,
#       Min = value_range[1],
#       Max = value_range[2]
#     )
#   }
#   
#   # Combine all ranges into one data frame
#   combined_ranges <- do.call(rbind, range_results_list)
#   
#   return(combined_ranges)
# }
# combined_ranges <- run_gmm_classification_and_extract_ranges()
# summary_threshold <- combined_ranges%>%summarize(Min = mean(Min), Max = mean(Max))
# #min -1.46 max -0.827
# filtered_results <- results %>%
#   filter(
#     !(label == 'functional' & ScaledFoldChange_mean < -1.46) &
#       !(label == 'non-functional' & ScaledFoldChange_mean > -0.827) &
#       !(
#         label == 'intermediate' &
#           (posterior.comp.1 < 0.25 & ScaledFoldChange_mean <= -1.25)
#       ) &
#       !(
#         label == 'intermediate' &
#           (ScaledFoldChange_mean >= -1.25 & posterior.comp.1 > 0.75)
#       )
#   )

mixmdl =  mixtools::normalmixEM(data_5g_13g_nonsyn$ScaledFoldChange_mean, k = 2, maxit = 1000, epsilon = 1e-8, verb = FALSE)
ordering <- order(mixmdl$mu)
mu1 <- mixmdl$mu[ordering[1]]
mu2 <- mixmdl$mu[ordering[2]]
sigma1 <- mixmdl$sigma[ordering[1]]
sigma2 <- mixmdl$sigma[ordering[2]]
lambda1 <- mixmdl$lambda[ordering[1]]
lambda2 <- mixmdl$lambda[ordering[2]]

posterior1_fun <- function(x) {
  num <- lambda1 * dnorm(x, mean = mu1, sd = sigma1)
  denom <- num + lambda2 * dnorm(x, mean = mu2, sd = sigma2)
  return(num / denom)
}

score_range <- range(data_5g_13g_nonsyn$ScaledFoldChange_mean, na.rm = TRUE)
lower_cut =  uniroot(function(x) posterior1_fun(x) - 0.99, interval = score_range, extendInt = "yes")$root
upper_cut = uniroot(function(x) posterior1_fun(x) - 0.01, interval = score_range, extendInt = "yes")$root


###
classification = c()
for (i in 1:nrow(edited_master_table)) {
  score = edited_master_table$new_score_gdna_avg[i]
  if (is.na(score)) {classification = c(classification, NA); next}
  if (score > upper_cut) {classification = c(classification,"functional")}
  else if (score < lower_cut) {classification = c(classification,"non-functional")}
  else(classification = c(classification,"intermediate"))
}

edited_master_table$classification = classification
#write.xlsx(edited_master_table,"master_table_with_new_scores_and_classifications_71425.xlsx")
up_to_date_annos = read.xlsx("master_table_with_new_scores_and_classifications_and_clinvar_31025.xlsx")
match_to_master_annotations = c()
for (m in edited_master_table$Mutcode) { 
 
  match_to_master_annotations = c(match_to_master_annotations,up_to_date_annos$ClinVar[which(up_to_date_annos$Mutcode==m)[1]])
}
edited_master_table_fix_annos = edited_master_table

edited_master_table_fix_annos$ClinVar = match_to_master_annotations
write.xlsx(edited_master_table_fix_annos,"master_table_with_new_scores_and_classifications_and_clinvar_82625.xlsx")

### per david - classification algorithm with 0.95/0.05?
thresholds_5g_13g_nonsyn_2<- c(.95, .05)
functional_classification_5g_13g_nonsyn_2 <- data.frame(logFC_adj_avg =data_5g_13g_nonsyn$ScaledFoldChange_mean, 
                                                      posterior = mixmdl_5g_13g_nonsyn$posterior,
                                                      label = ifelse(
                                                        posterior_5g_13g_nonsyn[, 1] > thresholds_5g_13g_nonsyn_2[1], 'non-functional',
                                                        ifelse(
                                                          posterior_5g_13g_nonsyn[, 1] < thresholds_5g_13g_nonsyn_2[2], 'functional', 'intermediate'
                                                        )
                                                      )
)
data_thresholdnon_syn <- data_5g_13g_nonsyn%>%select(ScaledFoldChange_mean)
set.seed(72) # Ensure reproducibility

# Define a function to fit the GMM, classify data, and repeat the process
run_gmm_classification <- function( num_iterations = 1000) {
  # Initialize a data frame to store results
  new_point_labels <- data.frame(
    ScaledFoldChange_mean = numeric(),
    posterior.comp.1 = numeric(),
    posterior.comp.2 = numeric(),
    label = character(),
    iteration = integer(),
    stringsAsFactors = FALSE
  )
  for (iteration in 1:num_iterations) {
    # Add a random data point within the specified range
    new_point <- data.frame(ScaledFoldChange_mean = runif(1, min = -2, max = -0.5))
    data <- rbind(data_thresholdnon_syn , new_point)
    
    # Fit the Gaussian mixture model with 2 components
    mixmdl <- normalmixEM(data$ScaledFoldChange_mean, k = 2)
    posterior <- mixmdl$posterior
    
    # Define thresholds for classifying functional SNVs
    thresholds <- c(0.95, 0.05)
    
    # Classify SNVs based on posterior probabilities for the new point
    new_point_posterior <- posterior[nrow(posterior), ]
    
    # Determine the label for the new point
    new_point_label <- ifelse(
      new_point_posterior[1] > thresholds[1], 'non-functional',
      ifelse(new_point_posterior[1] < thresholds[2], 'functional', 'intermediate')
    )
    
    # Store only the new point's classification
    new_point_labels <- rbind(new_point_labels, data.frame(
      ScaledFoldChange_mean = new_point$ScaledFoldChange_mean,
      posterior.comp.1 = new_point_posterior[1],
      posterior.comp.2 = new_point_posterior[2],
      label = new_point_label,
      iteration = iteration
    ))
    
  }
  return(new_point_labels)
}

results <- run_gmm_classification()
run_gmm_classification_and_extract_ranges <- function() {
  # Initialize an empty list to store range results
  range_results_list <- vector("list", 25)
  
  # Run the classification 100 times
  for (i in 1:25) {
    # Run the classification
    result <- run_gmm_classification()
    
    # Filter the results to get only the rows labeled "intermediate"
    intermediate_results <- result[result$label == "intermediate", ]
    
    # Check if there are any intermediate results
    if (nrow(intermediate_results) > 0) {
      # Calculate the range of values in the "intermediate" label column
      value_range <- range(intermediate_results$ScaledFoldChange_mean, na.rm = TRUE)
    } else {
      # If no intermediate results, set the range to NA
      value_range <- c(NA, NA)
    }
    
    # Store the range in the list
    range_results_list[[i]] <- data.frame(
      Run = i,
      Min = value_range[1],
      Max = value_range[2]
    )
  }
  
  # Combine all ranges into one data frame
  combined_ranges <- do.call(rbind, range_results_list)
  
  return(combined_ranges)
}
combined_ranges <- run_gmm_classification_and_extract_ranges()
summary_threshold <- combined_ranges%>%summarize(Min = mean(Min), Max = mean(Max))
#min -1.39 max -1.02
filtered_results_2 <- results %>%
  filter(
    !(label == 'functional' & ScaledFoldChange_mean < -1.39) &
      !(label == 'non-functional' & ScaledFoldChange_mean > -1.02) &
      !(
        label == 'intermediate' &
          (posterior.comp.1 < 0.25 & ScaledFoldChange_mean <= -1.39)
      ) &
      !(
        label == 'intermediate' &
          (ScaledFoldChange_mean >= -1.39 & posterior.comp.1 > 0.75)
      )
  )

classification = c()
for (i in 1:nrow(master_with_variants_to_test)) {
  score = master_with_variants_to_test$new_score_gdna_avg[i]
  if (is.na(score)) {classification = c(classification, NA); next}
  if (score > -1.02) {classification = c(classification,"functional")}
  else if (score < -1.39) {classification = c(classification,"non-functional")}
  else(classification = c(classification,"intermediate"))
}

master_with_variants_to_test$classification_2_TEST = classification

