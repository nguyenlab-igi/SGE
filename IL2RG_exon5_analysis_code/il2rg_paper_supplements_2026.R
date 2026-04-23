library(openxlsx)
library(dplyr)
og_annotation = read.xlsx("master_table_with_new_scores_and_classifications_and_clinvar_42226.xlsx")
og_annotation$clinical_numbering = og_annotation$`cDNA.#`-38
og_annotation$base2 = toupper(sapply(strsplit(og_annotation$base, " "), "[", 2))
##get hgvs format
og_annotation$clinical_numbering = og_annotation$`cDNA.#`-38
og_annotation$base2 = toupper(sapply(strsplit(og_annotation$base, " "), "[", 2))
mutation_coords <- og_annotation %>%
  group_by(clinical_numbering) %>%
  mutate(
    wt_base = base2[Mutation.type == "WT_"][1]
  ) %>%
  ungroup() %>%
  #filter(Mutation.type != "WT_") %>%
  mutate(
    coord = paste0("c.", clinical_numbering, wt_base, ">", base2)
  ) %>%
  select(clinical_numbering, wt_base, new_base = base2, coord)

og_annotation$hgvs = mutation_coords$coord
nwt = og_annotation[-which(og_annotation$Mutation.type == "WT_"),]
cleaned_table = nwt[,c(32,42,9,33,34,35,36,37,38,39,22,31)]
colnames(cleaned_table)=c("AA_Change|Codon","HGVS_cDNA","Type","Score_Donor1_gDNA","Score_Donor2_gDNA","Score_Average_gDNA",
                          "Score_Donor1_cDNA","Score_Donor2_cDNA","Score_Average_cDNA","Screen_Classification","DB_Annotation","Template")
#write.xlsx(cleaned_table,"screen_data_gdna_and_cdna_formatted_12026.xlsx")

mutation_colors <- c("Syn" = "#BAE0AD", "Mis" = "#6D91BA", "Non" = "#B3776F")

summary_stats <- function(y) {
  r <- quantile(y, probs = c(0.25, 0.5, 0.75))
  names(r) <- c("ymin", "y", "ymax")
  r
}
my_comparisons <- list( c("Syn", "Non"), c("Syn", "Mis"), c("Mis", "Non") )
library(ggpubr)
for_means = compare_means(new_score_cdna_avg ~ Mutation.type,  data = non_wild_types,method = "t.test")
stat_df <- compare_means(
  new_score_cdna_avg ~ Mutation.type,
  data = non_wild_types,
  comparisons = my_comparisons
)

# Create formatted label
stat_df$label <- paste0("p=", signif(stat_df$p, 2))
stat_df$y.position <- c(6, 5, 4)
#

# ekansh_loess_2 = read.csv("data_5g_13g_combined_with_cons.csv")
# ekansh_loess_2$Mutcode <- gsub("\\|.*-> ", "|", ekansh_loess_2$Mutcode2)
# ekansh_loess_2 = ekansh_loess_2[-which(ekansh_loess_2$Mutcode %in% excluded_snplist$Mutcode),]
# ekansh_loess_2$Mutation.type <- factor(ekansh_loess_2$Mutation.type, levels = c("Syn", "Mis", "Non"))
ggplot(non_wild_types, aes(x = Mutation.type, y = new_score_cdna_avg)) +
  geom_jitter(aes(color = Mutation.type), width = 0.2, size = 3, alpha = 0.6) +  # Add transparency to jittered points
  stat_summary(fun.data = summary_stats, geom = "errorbar",width=0.1, color = "black") +  # Increase size of error bars
  stat_summary(fun = median, geom = "point", size = 1, color = "black") +  # Increase size of median points
  #stat_compare_means(comparisons = my_comparisons) +
  scale_color_manual(values = mutation_colors) +
  labs(x = "", y = "Score", title = "cDNA sequencing") +  # Add title
  theme_minimal(base_size = 7.5) + stat_pvalue_manual(
    stat_df,
    label = "label",
    tip.length = 0.01
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 15),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +  # Add horizontal line at y = 0
  scale_y_continuous(breaks = seq(-12, 7, by = 1), limits = c(-12, 7))  # Adjust y-axis limits and breaks


ggsave("figures3_donor_average_scores_cdna_42226.pdf",device = "pdf")

dp_gdna_d1 = read.delim("diff_exp_results_annotated_lib_with_remaining_variant/annotated_21g_1Libp_S21_L001_VS_LibTemplate_S2_L001.csv")
match_to_master = c()
for (i in 1:nrow(non_wild_types)) {
  seq = non_wild_types$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG[i]
  pos = which(toupper(dp_gdna_d1$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG)==toupper(seq))
  if (length(pos)==0) { match_to_master = c(match_to_master,NA); next}
  match_to_master = c(match_to_master,dp_gdna_d1$fold_change[pos])
}
non_wild_types$dp_gdna_d1 = match_to_master
dp_gdna_d2 = read.delim("diff_exp_results_annotated_lib_with_remaining_variant/annotated_29g_2Libp_S28_L001_VS_LibTemplate_S2_L001.csv")
match_to_master = c()
for (i in 1:nrow(non_wild_types)) {
  seq = non_wild_types$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG[i]
  pos = which(toupper(dp_gdna_d2$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG)==toupper(seq))
  if (length(pos)==0) { match_to_master = c(match_to_master,NA); next}
  match_to_master = c(match_to_master,dp_gdna_d2$fold_change[pos])
}
non_wild_types$dp_gdna_d2 = match_to_master
non_wild_types$dp_gdna_mean = rowMeans(cbind(non_wild_types$dp_gdna_d1,non_wild_types$dp_gdna_d2))
data_21g_29g_combined <- data.frame(
  LogFC_21g = non_wild_types$dp_gdna_d1,
  LogFC_29g = non_wild_types$dp_gdna_d2
)
data_21g_29g_combined <- data_21g_29g_combined %>%
  mutate(
    avgFC = (LogFC_21g + LogFC_29g) / 2,
    #varFC = (LogFC_3m - avgFC)^2,
    #sdFC = sqrt(varFC),
    `cDNA..` = non_wild_types$`cDNA.#`,
    `pro..` = non_wild_types$`pro.#`,
    Mutation.type = non_wild_types$Mutation.type,
    ClinVar = non_wild_types$ClinVar, 
    Mutcode = non_wild_types$Mutcode, 
    ClinVar = coalesce(ClinVar, "Unclassified"))
data_21g_Loess <- data_21g_29g_combined%>%filter(LogFC_21g >= min(data_21g_29g_combined$LogFC_21g[data_21g_29g_combined$Mutation.type=="Syn"]), !(Mutcode %in% exclusion_list))
loess_15_21g <- loess(LogFC_21g ~ `cDNA..`, data = data_21g_Loess, span = 0.15)
data_21g_Loess$smoothed15 <- predict(loess_15_21g) 


#########LOESSfit for 11m (Donor2)
data_29g_Loess <- data_21g_29g_combined%>%filter(LogFC_29g >= min(data_21g_29g_combined$LogFC_29g[data_21g_29g_combined$Mutation.type=="Syn"]), !(Mutcode %in% exclusion_list))
loess_15_29g <- loess(LogFC_29g ~ `cDNA..`, data = data_29g_Loess, span = 0.15)
data_29g_Loess$smoothed15 <- predict(loess_15_29g) 


#adding LOESS predictions to data set 
loess_predictions_21g <- data.frame(predict(loess_15_21g , newdata = data.frame(`cDNA..` = 633:795)))
colnames(loess_predictions_21g) = "LOESS_21g"
loess_predictions_21g <- loess_predictions_21g%>%mutate(`cDNA..` = 633:795)

loess_predictions_29g <- data.frame(predict(loess_15_29g , newdata = data.frame(`cDNA..` = 633:795)))
colnames(loess_predictions_29g) = "LOESS_29g"
loess_predictions_29g <- loess_predictions_29g%>%mutate(`cDNA..` = 633:795)

data_21g_29g_combined<- left_join(data_21g_29g_combined, loess_predictions_21g , by = "cDNA..")
data_21g_29g_combined<- left_join(data_21g_29g_combined, loess_predictions_29g , by = "cDNA..")
data_21g_29g_combined <-data_21g_29g_combined%>%mutate(logFC_21g_adj = (LogFC_21g - LOESS_21g))%>%mutate(logFC_29g_adj = (LogFC_29g - LOESS_29g))
data_21g_29g_combined<-data_21g_29g_combined%>%mutate(logFC_adj_avg = (logFC_21g_adj+logFC_29g_adj)/2)
data_21g_29g_combined_syn = data_21g_29g_combined[which(data_21g_29g_combined$Mutation.type=="Syn"),]
data_21g_29g_combined_non = data_21g_29g_combined[which(data_21g_29g_combined$Mutation.type=="Non"),]
syn_scale_value = mean(c(median(data_21g_29g_combined_syn$logFC_21g_adj,na.rm=T),median(data_21g_29g_combined_syn$logFC_29g_adj,na.rm=T)))
non_scale_value = mean(c(median(data_21g_29g_combined_non$logFC_21g_adj,na.rm=T),median(data_21g_29g_combined_non$logFC_29g_adj,na.rm=T)))

syn_median_d1 <- median(data_21g_29g_combined$logFC_21g_adj[data_21g_29g_combined$Mutation.type == "Syn"],na.rm=T)
syn_median_d2 <- median(data_21g_29g_combined$logFC_29g_adj[data_21g_29g_combined$Mutation.type == "Syn"],na.rm=T)

non_median_d1 <- median(data_21g_29g_combined$logFC_21g_adj[data_21g_29g_combined$Mutation.type == "Non"],na.rm=T)
non_median_d2 <- median(data_21g_29g_combined$logFC_29g_adj[data_21g_29g_combined$Mutation.type == "Non"],na.rm=T)

scale_factor_d1 <- (non_scale_value - syn_scale_value) / (non_median_d1 - syn_median_d1)
shift_factor_d1 <- syn_scale_value - scale_factor_d1 * syn_median_d1


scale_factor_d2 <- (non_scale_value - syn_scale_value) / (non_median_d2 - syn_median_d2)
shift_factor_d2 <- syn_scale_value - scale_factor_d2 * syn_median_d2

# Apply the transformation to the fold change values
data_21g_29g_combined$ScaledFoldChange_d1 <- data_21g_29g_combined$logFC_21g_adj * scale_factor_d1 + shift_factor_d1
data_21g_29g_combined$ScaledFoldChange_d2 <- data_21g_29g_combined$logFC_29g_adj * scale_factor_d2 + shift_factor_d2
data_21g_29g_combined$ScaledFoldChange_mean = rowMeans(cbind(data_21g_29g_combined$ScaledFoldChange_d1,data_21g_29g_combined$ScaledFoldChange_d2))
#non_wild_types = data_21g_29g_combined[-which(data_21g_29g_combined$Mutation.type=="WT_"),]
mutation_colors <- c("Syn" = "#BAE0AD", "Mis" = "#6D91BA", "Non" = "#B3776F")

stat_df <- compare_means(
  ScaledFoldChange_mean ~ Mutation.type,
  data = data_21g_29g_combined,
  comparisons = my_comparisons
)

stat_df$label <- paste0("p=", signif(stat_df$p, 2))
stat_df$y.position <- c(6, 5, 4)
#

ggplot(data_21g_29g_combined, aes(x = Mutation.type, y = ScaledFoldChange_mean)) +
  geom_jitter(aes(color = Mutation.type), width = 0.2, size = 3, alpha = 0.6) +  # Add transparency to jittered points
  stat_summary(fun.data = summary_stats, geom = "errorbar",width=0.25, color = "black") +  # Increase size of error bars
  stat_summary(fun = median, geom = "point", size = 1, color = "black") +  # Increase size of median points
  #stat_compare_means(comparisons = my_comparisons,method = "t.test") +
  scale_color_manual(values = mutation_colors) +
  labs(x = "", y = "Score", title = "Late Timepoint") +  # Add title
  theme_minimal(base_size = 7.5) + stat_pvalue_manual(
    stat_df,
    label = "label",
    tip.length = 0.01
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 15),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +  # Add horizontal line at y = 0
  scale_y_continuous(breaks = seq(-10, 6, by = 1), limits = c(-10, 6))  # Adjust y-axis limits and breaks
ggsave("figures3_donor_average_scores_late_timepoint_12126.pdf",device = "pdf")

cleaned_table$Score_Day14_D1 <- data_21g_29g_combined$ScaledFoldChange_d1[match(cleaned_table$`AA_Change|Codon`, data_21g_29g_combined$Mutcode)]
cleaned_table$Score_Day14_D2 <- data_21g_29g_combined$ScaledFoldChange_d2[match(cleaned_table$`AA_Change|Codon`, data_21g_29g_combined$Mutcode)]
cleaned_table$Score_Day14_Average <- data_21g_29g_combined$ScaledFoldChange_mean[match(cleaned_table$`AA_Change|Codon`, data_21g_29g_combined$Mutcode)]

early_late = merge(non_wild_types,data_21g_29g_combined,by="Mutcode")
ggplot(early_late,aes(x = new_score_gdna_avg, y = ScaledFoldChange_mean)) +
  geom_point(aes(color = Mutation.type.x),size=3) +
  labs(x = "Day 3 Score (Two Donors)", y = "Day 14 Score (Two Donors)",caption="r=0.74") +
  scale_color_manual(values = c("Mis" = "#1f77b4", "Non" = "#e63f29", "Syn" = "#42bd0d"))+
  theme_minimal(base_size = 15) + xlim(-12,4)+ geom_smooth(method='lm', formula= y~x,se=F,fullrange=T)+
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid = element_blank(), legend.text = element_text(size=20)# Remove gridlines
  )+ theme(axis.text = element_text(size = 15),  # Adjust axis text size
           axis.ticks = element_line(size = 0.5),
           panel.grid = element_blank(), axis.title.y.right = element_blank(),
           axis.title.x = element_text(size = 20, , margin = margin(t = 15)),
           axis.title.y = element_text(size = 20, margin = margin(r = 15)),
           aspect.ratio = .7, plot.title = element_text(hjust = 0.5),
           plot.caption = element_text(hjust = 0.5, size = 15, face = "italic", margin = margin(t = 10))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")

ggsave("supp_compare_earlylate_42226.pdf",device = "pdf")


screen_18g = read.delim("corrected_annotated_DN_18g_14_500U_S18_L001_VS_DN_Lib_HDRT_S20_L001.csv")
screen_18g$Mutcode = paste(paste0(screen_18g$WTaa,screen_18g$pro..,screen_18g$mutaminoacid),screen_18g$mutcodon,sep="|")
data_18g_Loess <- screen_18g %>%filter(fold_change >= min(screen_18g$fold_change[screen_18g$Mutation.type=="Syn"]), !(Mutcode %in% exclusion_list))


loess_15_18g <- loess(fold_change ~ `cDNA..`, data = data_18g_Loess, span = 0.15)
data_18g_Loess$smoothed15 <- predict(loess_15_18g) 



#adding LOESS predictions to data set 
loess_predictions_18g <- data.frame(predict(loess_15_18g , newdata = data.frame(`cDNA..` = 633:795)))
colnames(loess_predictions_18g) = "LOESS_18g"
loess_predictions_18g <- loess_predictions_18g%>%mutate(`cDNA..` = 633:795)

screen_18g<- left_join(screen_18g, loess_predictions_18g , by = "cDNA..")
screen_18g <-screen_18g%>%mutate(logFC_18g = (fold_change - LOESS_18g))

screen_18g_nwt = screen_18g[which(screen_18g$Mutation.type %in% c("Mis","Syn","Non")),]

stat_df <- compare_means(
  logFC_18g ~ Mutation.type,
  data = screen_18g_nwt,
  comparisons = my_comparisons
)

stat_df$label <- paste0("p=", signif(stat_df$p, 2))
stat_df$y.position <- c(6, 5, 4)
#
cleaned_table$Score_Donor3 <- screen_18g_nwt$logFC_18g[match(cleaned_table$`AA_Change|Codon`, screen_18g_nwt$Mutcode)]

ggplot(screen_18g_nwt, aes(x = Mutation.type, y = logFC_18g)) +
  geom_jitter(aes(color = Mutation.type), width = 0.2, size = 3, alpha = 0.6) +  # Add transparency to jittered points
  stat_summary(fun.data = summary_stats, geom = "errorbar",width=0.25, color = "black") +  # Increase size of error bars
  #stat_summary(fun = median, geom = "point", size = 1, color = "black") +  # Increase size of median points
  #stat_compare_means(comparisons = my_comparisons) +
  scale_color_manual(values = mutation_colors) +
  labs(x = "", y = "Score", title = "Third Donor (no IL2RG+ sort)") +  # Add title
  theme_minimal(base_size = 7.5) +
  stat_pvalue_manual(
    stat_df,
    label = "label",
    tip.length = 0.01
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 15),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +  # Add horizontal line at y = 0
  scale_y_continuous(breaks = seq(-10, 6, by = 1), limits = c(-10, 6))  # Adjust y-axis limits and breaks


ggsave("figures3_donor_average_scores_third_donor_42226.pdf",device = "pdf")

screen_18g_nwt$first_donors = non_wild_types$new_score_gdna_avg
non_wild_types$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG = toupper(non_wild_types$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG)
three_donors = merge(screen_18g_nwt,non_wild_types,by="GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG")
ggplot(three_donors,aes(x = new_score_gdna_avg, y = logFC_18g)) +
  geom_point(aes(color = Mutation.type.x),size=3) +
  labs(x = "Donors 1 and 2 Average Score", y = "Donor 3 Score",caption="r=0.85") +
  scale_color_manual(values = c("Mis" = "#1f77b4", "Non" = "#e63f29", "Syn" = "#42bd0d"))+
  theme_minimal(base_size = 15) + xlim(-12,4)+ geom_smooth(method='lm', formula= y~x,se=F,fullrange=T)+
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid = element_blank(), legend.text = element_text(size=20)# Remove gridlines
  )+ theme(axis.text = element_text(size = 15),  # Adjust axis text size
           axis.ticks = element_line(size = 0.5),
           panel.grid = element_blank(), axis.title.y.right = element_blank(),
           axis.title.x = element_text(size = 20, , margin = margin(t = 15)),
           axis.title.y = element_text(size = 20, margin = margin(r = 15)),
           aspect.ratio = .7, plot.title = element_text(hjust = 0.5),
           plot.caption = element_text(hjust = 0.5, size = 15, face = "italic", margin = margin(t = 10))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")

ggsave("figures3_compare_to_third_donor_42226.pdf",device = "pdf")


supp_dotplot = ggplot(screen_18g_nwt, aes(x = Mutation.type, y = logFC_18g)) +
  geom_jitter(aes(color = Mutation.type), width = 0.2, size = 3, alpha = 0.6) +  # Add transparency to jittered points
  stat_summary(fun.data = summary_stats, geom = "errorbar",width=0.25, color = "black") +  # Increase size of error bars
  #stat_summary(fun = median, geom = "point", size = 1, color = "black") +  # Increase size of median points
  #stat_compare_means(comparisons = my_comparisons) +
  scale_color_manual(values = mutation_colors) +
  labs(x = "", y = "Score", title = "Third Donor (no IL2RG+ sort)",color="") +  # Add title
  theme_minimal(base_size = 7.5) +
  stat_pvalue_manual(
    stat_df,
    label = "label",
    tip.length = 0.01
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +  # Add horizontal line at y = 0
  scale_y_continuous(breaks = seq(-10, 6, by = 1), limits = c(-10, 6))  # Adjust y-axis limits and breaks

supp_dotplot_legend = get_legend(supp_dotplot)
as_ggplot(supp_dotplot_legend)
ggsave("figs3_legend.pdf",device = "pdf")

mutation_colors <- c(
  "Benign" = "#00ff00", 
  "Likely Benign" = "lightgreen", 
  "VUS" = "#C11C84", 
  "Likely Pathogenic" = "#de8bf9",  # A lighter red for Likely Pathogenic
  "X-SCID" = "#f94040", 
  "Unclassified" = "lightgrey") 


ggplot(screen_w_cadd,aes(x=new_score_gdna_avg,y=PolyPhenVal)) + geom_point(aes(color = ClinVar), size = 4, alpha = 0.8) +
  scale_color_manual(values = mutation_colors) +
  labs(x = "Screen Score", y = "PolyPhen Score", title = "") +  # Add title
  theme_minimal(base_size = 7.5) +
  #geom_smooth(method = "loess", formula = y ~ x, color = "blue", se = FALSE,size=0.5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 15,face="bold"), legend.text = element_text(size=10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),legend.title=element_blank()
  ) 

ggsave("figs4_compare_to_polyphen_42226.pdf",device = "pdf")

ggplot(screen_w_cadd,aes(x=new_score_gdna_avg,y=SIFTval)) + geom_point(aes(color = ClinVar), size = 4, alpha = 0.8) +
  scale_color_manual(values = mutation_colors) +
  labs(x = "Screen Score", y = "SIFT Score", title = "") +  # Add title
  theme_minimal(base_size = 7.5) +
  #geom_smooth(method = "loess", formula = y ~ x, color = "blue", se = FALSE,size=0.5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 15,face="bold"), legend.text = element_text(size=10),legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),legend.title=element_blank()
  ) 
ggsave("figs4_compare_to_sift_42226.pdf",device = "pdf")

ratios = read.csv("~/Downloads/count files/ratios_test_benign_32_mutations.csv")
variants = ratios$aa_change
#G247E and G247A did not have a working hdrt value so I will exclude them
big_table_t= c()
for (v in variants){
  if (v %in% c("H236H","G247E","G247A","G247G")) {next}
  subset = ratios[which(ratios$aa_change==v),-c(1,2,3)]
  values = c(subset$ratio_hdrt,
             mean(c(subset$ratio_day3_d1,subset$ratio_day3_d2),na.rm=T),
             mean(c(subset$ratio_day5_no_il2_d1,subset$ratio_day5_no_il2_d2),na.rm=T),
             mean(c(subset$ratio_day5_il2_d1,subset$ratio_day5_il2_d2),na.rm=T),
             mean(c(subset$ratio_day8_no_il2_d1,subset$ratio_day8_no_il2_d2),na.rm=T),
             mean(c(subset$ratio_day8_il2_d1,subset$ratio_day8_il2_d2),na.rm=T))
  values = values/values[1]
  category = c("HDRT","Day 3","Day 5 IL2-","Day 5 IL2+","Day 8 IL2-","Day 8 IL2+")
  variant = rep(v,length(values)) 
  variant_out = data.frame(category,values,variant)
  big_table_t = rbind(big_table_t,variant_out)
}

big_table_t$category = factor(big_table_t$category,levels = c("HDRT","Day 3","Day 5 IL2-","Day 5 IL2+","Day 8 IL2-","Day 8 IL2+"))
big_table_t$variant = factor(big_table_t$variant,levels = rev(c("F227C","R222P","R226C","R222G","R226H","S225R","W237*",
                                                            "W237R","W237S","A234P","E239*","R226P","R224G",
                                                            "R224P","W246S","R226G","A234A","G215E","H245D",
                                                            "H245L","H245H","H245N","I244I","I244V","K252R","S233T","G247R","G232R")))
ggplot(big_table_t, aes(x = category, y = variant, fill = values)) +
  geom_tile() + 
  #scale_fill_gradient(low = "white", high = "black") +  # Grayscale gradient
  scale_fill_gradientn(
    colors = c("white", "black"),  # White to black gradient
    values = c(0, 1.5),       # Define stops in data range
    limits = c(0, 1.5),            # Clamp values at 1.5
    oob = scales::squish           # Anything above 1.5 stays black
  ) +
  
  theme_minimal() +
  labs(title = "",
       x = "Condition",
       y = "Variant",
       fill = "Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=25,face="bold"),
        axis.text.y = element_text(angle = 45, hjust = 1,size=25,face="bold"),
        axis.title.x = element_text(size = 25,face="bold"),
        axis.title.y = element_text(size = 25,face="bold"),
        legend.text=element_text(size=20,face="bold"),legend.title=element_text(size=20,face="bold"))
ggsave("figures5_heatmap.pdf",dev="pdf",height=15,width=10)
#ggsave("t_cell_differentiation_31025.eps",dev="eps",width=15)

#og_annotation = read.xlsx("master_table_with_new_scores_and_classifications_after_excluding_11724.xlsx")
#og_annotation$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG = toupper(og_annotation$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG)
#technically 69 is probably D1 and 87 D2
dp_gdna_d1 = read.csv("diff_exp_output_hdrt_add_1_to_0_counts/69g_A1_TN_TG_S71_L004_VS_H2_HDRT-ssLib_PG_S159_L004.csv")

match_to_master = c()
for (i in 1:nrow(non_wild_types)) {
  seq = non_wild_types$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG[i]
  pos = which(rownames(dp_gdna_d1)==toupper(seq))
  if (length(pos)==0) { match_to_master = c(match_to_master,NA); next}
  match_to_master = c(match_to_master,dp_gdna_d1$logFC[pos])
}
non_wild_types$dp_gdna_d1 = match_to_master

dp_gdna_d2 = read.csv("diff_exp_output_hdrt_add_1_to_0_counts/87g_B1_TN_TG_S89_L004_VS_H2_HDRT-ssLib_PG_S159_L004.csv")
match_to_master = c()
for (i in 1:nrow(non_wild_types)) {
  seq = non_wild_types$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG[i]
  pos = which(rownames(dp_gdna_d2)==toupper(seq))
  if (length(pos)==0) { match_to_master = c(match_to_master,NA); next}
  match_to_master = c(match_to_master,dp_gdna_d2$logFC[pos])
}
non_wild_types$dp_gdna_d2 = match_to_master
non_wild_types$dp_gdna_mean = rowMeans(cbind(non_wild_types$dp_gdna_d1,non_wild_types$dp_gdna_d2))

data_69g_87g_combined <- data.frame(
  LogFC_69g = non_wild_types$dp_gdna_d1,
  LogFC_87g = non_wild_types$dp_gdna_d2
)
data_69g_87g_combined <- data_69g_87g_combined %>%
  mutate(
    avgFC = (LogFC_69g + LogFC_87g) / 2,
    #varFC = (LogFC_3m - avgFC)^2,
    #sdFC = sqrt(varFC),
    `cDNA..` = non_wild_types$`cDNA.#`,
    `pro..` = non_wild_types$`pro.#`,
    Mutation.type = non_wild_types$Mutation.type,
    ClinVar = non_wild_types$ClinVar, 
    Mutcode = non_wild_types$Mutcode, 
    ClinVar = coalesce(ClinVar, "Unclassified"))
# minimum of the synonymous
data_69g_Loess <- data_69g_87g_combined%>%filter(LogFC_69g >=	
                                                   min(data_69g_87g_combined$LogFC_69g[which(data_69g_87g_combined$Mutation.type=="Syn")]),!(Mutcode %in% exclusion_list)) 
loess_15_69g <- loess(LogFC_69g ~ `cDNA..`, data = data_69g_Loess, span = 0.15)
data_69g_Loess$smoothed15 <- predict(loess_15_69g) 


#########LOESSfit for 11m (Donor2)
data_87g_Loess <- data_69g_87g_combined%>%filter(LogFC_87g >=	
                                                   min(data_69g_87g_combined$LogFC_87g[which(data_69g_87g_combined$Mutation.type=="Syn")]),!(Mutcode %in% exclusion_list)) 
loess_15_87g <- loess(LogFC_87g ~ `cDNA..`, data = data_87g_Loess, span = 0.15)
data_87g_Loess$smoothed15 <- predict(loess_15_87g) 


#adding LOESS predictions to data set 
loess_predictions_69g <- data.frame(predict(loess_15_69g , newdata = data.frame(`cDNA..` = 633:795)))
colnames(loess_predictions_69g) = "LOESS_69g"
loess_predictions_69g <- loess_predictions_69g%>%mutate(`cDNA..` = 633:795)

loess_predictions_87g <- data.frame(predict(loess_15_87g , newdata = data.frame(`cDNA..` = 633:795)))
colnames(loess_predictions_87g) = "LOESS_87g"
loess_predictions_87g <- loess_predictions_87g%>%mutate(`cDNA..` = 633:795)

data_69g_87g_combined<- left_join(data_69g_87g_combined, loess_predictions_69g , by = "cDNA..")
data_69g_87g_combined<- left_join(data_69g_87g_combined, loess_predictions_87g , by = "cDNA..")
data_69g_87g_combined <-data_69g_87g_combined%>%mutate(logFC_69g_adj = (LogFC_69g - LOESS_69g))%>%mutate(logFC_87g_adj = (LogFC_87g - LOESS_87g))
data_69g_87g_combined<-data_69g_87g_combined%>%mutate(logFC_adj_avg = (logFC_69g_adj+logFC_87g_adj)/2)
data_69g_87g_combined_syn = data_69g_87g_combined[which(data_69g_87g_combined$Mutation.type=="Syn"),]
data_69g_87g_combined_non = data_69g_87g_combined[which(data_69g_87g_combined$Mutation.type=="Non"),]
syn_scale_value = mean(c(median(data_69g_87g_combined_syn$logFC_69g_adj,na.rm=T),median(data_69g_87g_combined$logFC_87g_adj,na.rm=T)))
non_scale_value = mean(c(median(data_69g_87g_combined_non$logFC_69g_adj,na.rm=T),median(data_69g_87g_combined$logFC_87g_adj,na.rm=T)))

syn_median_d1 <- median(data_69g_87g_combined$logFC_69g_adj[data_69g_87g_combined$Mutation.type == "Syn"],na.rm=T)
syn_median_d2 <- median(data_69g_87g_combined$logFC_87g_adj[data_69g_87g_combined$Mutation.type == "Syn"],na.rm=T)

non_median_d1 <- median(data_69g_87g_combined$logFC_69g_adj[data_69g_87g_combined$Mutation.type == "Non"],na.rm=T)
non_median_d2 <- median(data_69g_87g_combined$logFC_87g_adj[data_69g_87g_combined$Mutation.type == "Non"],na.rm=T)

# Define the linear transformation parameters
scale_factor_d1 <- (non_scale_value - syn_scale_value) / (non_median_d1 - syn_median_d1)
shift_factor_d1 <- syn_scale_value - scale_factor_d1 * syn_median_d1


scale_factor_d2 <- (non_scale_value - syn_scale_value) / (non_median_d2 - syn_median_d2)
shift_factor_d2 <- syn_scale_value - scale_factor_d2 * syn_median_d2

# Apply the transformation to the fold change values
data_69g_87g_combined$ScaledFoldChange_d1 <- data_69g_87g_combined$logFC_69g_adj * scale_factor_d1 + shift_factor_d1
data_69g_87g_combined$ScaledFoldChange_d2 <- data_69g_87g_combined$logFC_87g_adj * scale_factor_d2 + shift_factor_d2
data_69g_87g_combined$ScaledFoldChange_mean = rowMeans(cbind(data_69g_87g_combined$ScaledFoldChange_d1,data_69g_87g_combined$ScaledFoldChange_d2))
#non_wild_types = data_69g_87g_combined[-which(data_69g_87g_combined$Mutation.type=="WT_"),]
mutation_colors <- c("Syn" = "#BAE0AD", "Mis" = "#6D91BA", "Non" = "#B3776F")

stat_df <- compare_means(
  ScaledFoldChange_mean ~ Mutation.type,
  data = data_69g_87g_combined,
  comparisons = my_comparisons
)

stat_df$label <- paste0("p=", signif(stat_df$p, 2))
stat_df$y.position <- c(8, 7, 6)

ggplot(data_69g_87g_combined, aes(x = Mutation.type, y = ScaledFoldChange_mean)) +
  geom_jitter(aes(color = Mutation.type), width = 0.2, size = 3, alpha = 0.6) +  # Add transparency to jittered points
  #stat_summary(fun.data = summary_stats, geom = "errorbar",width=0.25, color = "black") +  # Increase size of error bars
  #stat_summary(fun = median, geom = "point", size = 1, color = "black") +  # Increase size of median points
  #stat_compare_means(comparisons = my_comparisons,method = "t.test") +
  scale_color_manual(values = mutation_colors) +
  labs(x = "", y = "Score", title = "") +  # Add title
  theme_minimal(base_size = 7.5) + stat_pvalue_manual(
    stat_df,
    label = "label",
    tip.length = 0.01
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 15),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black")
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +  # Add horizontal line at y = 0
  scale_y_continuous(breaks = seq(-6, 9, by = 1), limits = c(-6, 9))  # Adjust y-axis limits and breaks

ggsave("figures8_pcam02_42326.pdf",device = "pdf")



hsc_validation_old_to_use_wk3 = rbind(hsc_validation_round0[grep("HDRT",rownames(hsc_validation_round0)),],
                                      hsc_validation_round1[grep("Wk3",rownames(hsc_validation_round1)),])
variants_round1 = colnames(hsc_validation_old_to_use)[-c(5,6,7)]
big_table_wk3 =c()
for (v in variants_round1) {
  subset = hsc_validation_old_to_use[grep(v,rownames(hsc_validation_old_to_use_wk3)),c(v,"H236H")]
  init_ratio = subset[grep("HDRT",rownames(subset)),v]/subset[grep("HDRT",rownames(subset)),"H236H"]
  subset["ratio"] = subset[v]/subset$H236H/init_ratio
  subset_to_use = subset[grep('HDRT|DP|DN|B|_M|NK',rownames(subset)),]
  plotting_df = data.frame(category=sapply(strsplit(rownames(subset_to_use), "_"), "[", 3), value=subset_to_use$ratio)
  plotting_df$category[1] = "HDRT"
  plotting_df$variant = rep(v,nrow(plotting_df))
  big_table_wk3 = rbind(big_table_wk3,plotting_df)
}



big_table_wk3$category[which(big_table_wk3$category=="B")] = "CD19"
big_table_wk3$category[which(big_table_wk3$category=="M")] = "CD33"
big_table_wk3$category[which(big_table_wk3$category=="NK")] = "CD56"
big_table_wk3$category[which(big_table_wk3$category=="DP")] = "CD5+/CD7+ T"
big_table_wk3$category[which(big_table_wk3$category=="DN")] = "CD5-/CD7- T"

hsc_validation_old_to_use_wk6 = rbind(hsc_validation_round0[grep("HDRT",rownames(hsc_validation_round0)),],
                                      hsc_validation_round1[grep("Wk6",rownames(hsc_validation_round1)),])
variants_round1 = colnames(hsc_validation_old_to_use)[-c(5,6,7)]
big_table_wk6 =c()
for (v in variants_round1) {
  subset = hsc_validation_old_to_use_wk6[grep(v,rownames(hsc_validation_old_to_use_wk6)),c(v,"H236H")]
  init_ratio = subset[grep("HDRT",rownames(subset)),v]/subset[grep("HDRT",rownames(subset)),"H236H"]
  subset["ratio"] = subset[v]/subset$H236H/init_ratio
  subset_to_use = subset[grep('HDRT|DP|DN|B|_M|NK',rownames(subset)),]
  plotting_df = data.frame(category=sapply(strsplit(rownames(subset_to_use), "_"), "[", 3), value=subset_to_use$ratio)
  plotting_df$category[1] = "HDRT"
  plotting_df$variant = rep(v,nrow(plotting_df))
  big_table_wk6 = rbind(big_table_wk6,plotting_df)
}

big_table_wk6$category[which(big_table_wk6$category=="DP")] = "CD4+/CD8+ T"
big_table_wk6$category[which(big_table_wk6$category=="DN")] = "CD4-/CD8- T"

big_table_wk3_6 = rbind(big_table_wk3,big_table_wk6)

big_table_wk3_6$category = factor(big_table_wk3_6$category,levels = c("HDRT","CD5-/CD7- T",
                                                                      "CD5+/CD7+ T","CD4-/CD8- T","CD4+/CD8+ T","CD19","CD33","CD56"))
big_table_wk3_6$variant = factor(big_table_wk3_6$variant,levels = c("G232R","A234P",
                                                                    "A234A","S233T"))
ggplot(big_table_wk3_6, aes(x = category, y = variant, fill = value)) +
  geom_tile() + 
  #scale_fill_gradient(low = "white", high = "black") +  # Grayscale gradient
  scale_fill_gradientn(
    colors = c("white", "black"),  # White to black gradient
    values = c(0, 1.5),       # Define stops in data range
    limits = c(0, 1.5),            # Clamp values at 1.5
    oob = scales::squish           # Anything above 1.5 stays black
  ) +
  theme_minimal() +
  labs(title = "",
       x = "Cell Type",
       y = "Variant",
       fill = "Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=40,face="bold"),
        axis.text.y = element_text(angle = 45, hjust = 1,size=40,face="bold"),
        axis.title.x = element_text(size = 25,face="bold"),
        axis.title.y = element_text(size = 25,face="bold"),
        legend.text=element_text(size=20,face="bold"),legend.title=element_text(size=20,face="bold"))

ggsave("figures9_heatmap.pdf",dev="pdf",height=10,width=10)

conservation_data = read.csv("data_5g_13g_combined_with_cons.csv")
conservation_data$Mutcode1 = paste0(str_sub(conservation_data$Mutcode2, 1, 6),str_sub(conservation_data$Mutcode2, -3, -1))
match_to_master = c()
for (i in 1:nrow(non_wild_types)) {
  mutcode = non_wild_types$Mutcode[i]
  pos = which(toupper(conservation_data$Mutcode1)==mutcode)
  if (length(pos)==0) { match_to_master = c(match_to_master,NA); next}
  match_to_master = c(match_to_master,conservation_data$exon_conservation_scores[pos])
}



## make sure to use nwt from previous screen
non_wild_types$conservation = match_to_master
ggplot(non_wild_types, aes(x = new_score_gdna_avg, y = conservation)) +
  geom_point(aes(color = Mutation.type)) +
  labs(x = "Screen Score", y = "Conservation Score") +
  scale_color_manual(values = c("Mis" = "#1f77b4", "Non" = "#e63f29", "Syn" = "#42bd0d"))+
  theme_minimal(base_size = 15) + xlim(-12,4)+
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid = element_blank() # Remove gridlines
  )+ theme(axis.text = element_text(size = 10),  # Adjust axis text size
           axis.ticks = element_line(size = 0.5), 
           panel.grid = element_blank(), axis.title.y.right = element_blank(), 
           axis.title.x = element_text(size = 15, , margin = margin(t = 15)),
           axis.title.y = element_text(size = 15, margin = margin(r = 15)),
           aspect.ratio = .7, plot.title = element_text(hjust = 0.5), 
           plot.caption = element_text(hjust = 0.5, size = 10, face = "italic", margin = margin(t = 10))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")
ggsave("screen_vs_conservation_4625.pdf",device = "pdf")


mutation_colors <- c(
  "Benign" = "#00ff00", 
  "Likely Benign" = "lightgreen", 
  "VUS" = "#C11C84", 
  "Likely Pathogenic" = "#de8bf9",  # A lighter red for Likely Pathogenic
  "X-SCID" = "#f94040", 
  "Unclassified" = "lightgrey") 
#new_annos_no_wt = new_annos_no_wt[order(new_annos_no_wt$`pro.#`),]
#new_annos_no_wt$conservation = conservation_data$exon_conservation_scores
non_wild_types$ClinVar[which(non_wild_types$ClinVar == "X-SCID HMGD")] = "X-SCID"
classified = non_wild_types[which(non_wild_types$ClinVar != "Unclassified"),]
ggplot(classified, aes(x = new_score_gdna_avg, y = conservation)) +
  geom_point(aes(color = ClinVar),size=4) +
  labs(x = "Screen Score", y = "Conservation Score") +
  scale_color_manual(values = mutation_colors)+
  theme_minimal(base_size = 15) + xlim(-12,4)+
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid = element_blank() # Remove gridlines
  )+ theme(axis.text = element_text(size = 15,face="bold"),  # Adjust axis text size
           axis.ticks = element_line(size = 0.5), 
           panel.grid = element_blank(), axis.title.y.right = element_blank(), 
           axis.title.x = element_text(size = 20,face="bold" , margin = margin(t = 15)),
           axis.title.y = element_text(size = 20,face="bold", margin = margin(r = 15)),
           aspect.ratio = .7, plot.title = element_text(hjust = 0.5), legend.text = element_text(size=15,face="bold"))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")
ggsave("figuresx_conservation_vs_screen_42226.pdf",device = "pdf",height=10,width=10)
###
ex5_master_table = read.xlsx("master_table_with_new_scores_and_classifications_and_clinvar_82625.xlsx")
non_wild_types_ex5 = ex5_master_table[which(ex5_master_table$Mutation.type %in% c("Mis","Syn","Non")),]
non_wild_types_ex5$Mutation.type = factor(non_wild_types_ex5$Mutation.type)
mutation_colors <- c("Syn" = "#BAE0AD", "Mis" = "#6D91BA", "Non" = "#B3776F")
#all 
non_wild_types_ex5$Mutation.type = factor(non_wild_types_ex5$Mutation.type,levels=c("Non","Syn","Mis"))
ggplot(non_wild_types_ex5, aes(x =new_score_gdna_avg)) +
  geom_histogram(aes(fill=Mutation.type,y=..density..),bins = 95, alpha = 1) +
  scale_fill_manual(values = mutation_colors) +
  theme_minimal() + labs(x="Score",y="",fill="Mutation type") + geom_density(aes(fill=Mutation.type))
ggsave("figuresx_bimodality_31026.pdf",device="pdf")

#split up
cutoff1 <-c(-1.469)
cutoff2 <- c(-0.834)
non_wild_types_ex5 = non_wild_types #check this!
syn = non_wild_types_ex5[which(non_wild_types_ex5$Mutation.type== "Syn"),]
syn$Mutation.type = factor(syn$Mutation.type)
ggplot(syn, aes(x =new_score_gdna_avg)) +
  geom_histogram(aes(fill=Mutation.type,y=..density..),bins = 95, alpha = 1) +
  scale_fill_manual(values = mutation_colors) +
  theme_minimal() + xlim(-11,3)+ylim(0,1.5) + labs(x="Score",y="",fill="Mutation type") + geom_density(aes(fill=Mutation.type))

ggplot(syn, aes(x = new_score_gdna_avg)) + 
  geom_histogram(aes(fill=Mutation.type),bins = 95, alpha = 1) + 
  theme_minimal() + xlim(-11,3) + ylim(0,35) + labs(x="Score",y="",fill="Mutation type") + 
  scale_fill_manual(values = mutation_colors) +
  geom_vline(xintercept = cutoff1, linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = cutoff2, linetype = "dashed", color = "blue")   # Add horizontal line at y = 0
  

mis = non_wild_types_ex5[which(non_wild_types_ex5$Mutation.type== "Mis"),]
#syn$Mutation.type = factor(syn$Mutation.type)
ggplot(mis, aes(x =new_score_gdna_avg)) +
  geom_histogram(aes(fill=Mutation.type,y=..density..),bins = 95, alpha = 1) +
  scale_fill_manual(values = mutation_colors) +
  theme_minimal() + xlim(-11,3)+ylim(0,1.5) + labs(x="Score",y="",fill="Mutation type") + geom_density(aes(fill=Mutation.type))

ggplot(mis, aes(x = new_score_gdna_avg)) + 
  geom_histogram(aes(fill=Mutation.type),bins = 95, alpha = 1) + 
  theme_minimal() + xlim(-11,3) + ylim(0,35) + labs(x="Score",y="",fill="Mutation type") + 
  scale_fill_manual(values = mutation_colors) +
  geom_vline(xintercept = cutoff1, linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = cutoff2, linetype = "dashed", color = "blue")   # Add horizontal line at y = 0


non = non_wild_types_ex5[which(non_wild_types_ex5$Mutation.type== "Non"),]
#syn$Mutation.type = factor(syn$Mutation.type)
ggplot(non, aes(x =new_score_gdna_avg)) +
  geom_histogram(aes(fill=Mutation.type,y=..density..),bins = 95, alpha = 1) +
  scale_fill_manual(values = mutation_colors) +
  theme_minimal() + xlim(-11,3)+ylim(0,1.5) + labs(x="Score",y="",fill="Mutation type") + geom_density(aes(fill=Mutation.type))

ggplot(non, aes(x = new_score_gdna_avg)) + 
  geom_histogram(aes(fill=Mutation.type),bins = 95, alpha = 1) + 
  theme_minimal() + xlim(-11,3) + ylim(0,35) + labs(x="Score",y="",fill="Mutation type") + 
  scale_fill_manual(values = mutation_colors) +
  geom_vline(xintercept = cutoff1, linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = cutoff2, linetype = "dashed", color = "blue")   # Add horizontal line at y = 0

write.xlsx(cleaned_table,"screen_data_formatted_42326.xlsx")
