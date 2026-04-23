library(openxlsx)
#### HEATMAP ###
#no longer excluding
og_annotation = read.xlsx("master_table_with_new_scores_and_classifications_and_clinvar_42226.xlsx")
og_annotation$clinical_numbering = og_annotation$`cDNA.#`-38
wild_types = og_annotation[which(og_annotation$Mutation.type=="WT_"),]
wild_types$orig_base = toupper(sapply(strsplit(wild_types$base, " "), "[", 2))
wild_types_compact = wild_types[,c(35,41,6,40,27)]
wild_types_compact$new_score_gdna_avg=0
wild_types_compact$face = "bold"
wild_types_compact$color = "black"
wild_types_compact$border = "black"
non_wild_types = og_annotation[-which(og_annotation$Mutation.type=="WT_"),]
non_wild_types =  non_wild_types[which(!is.na(non_wild_types$new_score_gdna_avg)),]
non_wild_types$orig_base = toupper(sapply(strsplit(non_wild_types$base, " "), "[", 2))
non_wild_types_compact = non_wild_types[,c(35,41,6,40,27)]
non_wild_types_compact$face = "plain"
non_wild_types_compact$color = "black"
non_wild_types_compact$border=NA
colnames(wild_types_compact) = c("Score","Base","AminoAcid","DNA_Position","AA_Position","Text","Color","Border")
colnames(non_wild_types_compact) = c("Score","Base","AminoAcid","DNA_Position","AA_Position","Text","Color","Border")
heatmap_df = as.data.frame(rbind(rbind(non_wild_types_compact,wild_types_compact)))
heatmap_df_part1 = heatmap_df[which(heatmap_df$DNA_Position<=679),] #717 originally
heatmap_df_part1$Base <- as.factor(heatmap_df_part1$Base)
aa_positions_df = heatmap_df_part1[order(heatmap_df_part1$DNA_Position),c(1,4,5)]
aa_positions_subset <- aa_positions_df[seq(1, nrow(aa_positions_df), by = 60), ]
# Extract wild-type base positions for annotation
wild_type_labels <- wild_types_compact[, c("Base", "DNA_Position","Score")]

wild_type_labels <- unique(wild_type_labels)  # Remove duplicates
#pam changes

wild_type_labels$y_position <- max(as.numeric(factor(wild_type_labels$Base))) + 1
wild_type_labels$Base[which(wild_type_labels$DNA_Position==645)] = "G^"
wild_type_labels$Base[which(wild_type_labels$DNA_Position==741)] = "G^"
# Adjust y-position to place text above the heatmap
wild_type_labels_part1 = wild_type_labels[which(wild_type_labels$DNA_Position<=679),]

#library(showtext)
#font_add_google("Arial") # Arimo is a metrically compatible alternative
#showtext_auto()

ggplot(heatmap_df_part1, aes(x = DNA_Position, y = Base, fill = Score)) +
  geom_tile(aes(colour = as.factor(Border)), linewidth = 0.5) +
  geom_text(aes(label = AminoAcid,fontface=Text), colour=heatmap_df_part1$Color,size = 25/.pt) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 1.5, 
  #                     name = "Score") +
  scale_fill_gradientn(
    colors = c("red","white", "lightblue", "blue"),
    values = c(0,0.4, 0.75, 0.9, 1),  # -2 → 0, 0 → 0.5, 0.75 → 0.6875, 2 → 1
    limits = c(-2.5, 2),
    oob = scales::squish,
    name = "Score"
  ) +
  scale_colour_identity() +
  labs(x = "Position", y = "Base Change") +
  geom_text(
    data = aa_positions_subset,
    aes(x = DNA_Position, y = max(as.numeric(heatmap_df_part1$Base)) + 0.95, label = AA_Position),
    vjust = 2, size =8, color = "black", fontface="italic"
  ) +coord_cartesian(clip="off") +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.4)))  +
  
  geom_text(
    data = wild_type_labels_part1,
    aes(x = DNA_Position, label = Base,y=0),  # Removed y aesthetic from aes()
    vjust = 0, size = 21/.pt, fontface = "bold", color = "black"
  ) +
  theme_void() +guides(fill = guide_colorbar(barwidth = 1, barheight = 10)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold", margin = margin(t = 10),),
    axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 10),angle = 90),
    axis.text = element_text(size = 25,face="bold"), legend.title = element_text(size=25,face="bold",margin = margin(b = 20)),legend.text = element_text(size=25,face="bold")
  )

#ggsave("heatmap_bigger_font.pdf",device = "pdf",width=32,height=7)
#ggsave("heatmap_all_variants_top_12226.eps",device = "eps",width=30,height=5)
ggsave("heatmap_all_variants_top_42226.pdf",device = "pdf",width=32,height=7)
ggsave("heatmap_all_variants_top_42226.eps",device = "eps",width=32,height=7)


heatmap_df_part2 = heatmap_df[which(heatmap_df$DNA_Position>679),]
heatmap_df_part2$Base <- as.factor(heatmap_df_part2$Base)
aa_positions_df = heatmap_df_part2[order(heatmap_df_part2$DNA_Position),c(1,4,5)]
aa_positions_subset <- aa_positions_df[seq(1, nrow(aa_positions_df), by = 60), ]
wild_type_labels_part2 = wild_type_labels[which(wild_type_labels$DNA_Position>679),]

ggplot(heatmap_df_part2, aes(x = DNA_Position, y = Base, fill = Score)) +
  geom_tile(aes(colour = as.factor(Border)), linewidth = 0.5) +
  geom_text(aes(label = AminoAcid,fontface=Text), colour=heatmap_df_part2$Color,size = 25/.pt) +
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 1.5, 
  #                     name = "Score") +
  scale_fill_gradientn(
    colors = c("red","white", "lightblue", "blue"),
    values = c(0,0.4, 0.75, 0.9, 1),  # -2 → 0, 0 → 0.5, 0.75 → 0.6875, 2 → 1
    limits = c(-2.5, 2),
    oob = scales::squish,
    name = "Score"
  ) +
  scale_colour_identity() +
  labs(x = "Position", y = "Base Change") +
  geom_text(
    data = aa_positions_subset,
    aes(x = DNA_Position, y = max(as.numeric(heatmap_df_part2$Base)) + 0.95, label = AA_Position),
    vjust = 2, size =8, color = "black", fontface="italic"
  ) +coord_cartesian(clip="off") +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.4)))  +
  
  geom_text(
    data = wild_type_labels_part2,
    aes(x = DNA_Position, label = Base,y=0),  # Removed y aesthetic from aes()
    vjust = 0, size = 21/.pt, fontface = "bold", color = "black"
  ) +
  theme_void() +guides(fill = guide_colorbar(barwidth = 1, barheight = 10)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold", margin = margin(t = 10),),
    axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 10),angle = 90),
    axis.text = element_text(size = 25,face="bold"), legend.title = element_text(size=25,face="bold",margin = margin(b = 20)),legend.text = element_text(size=25,face="bold")
  )

#ggsave("heatmap_all_variants_bottom_12226.pdf",device = "pdf",width=30,height=5)
#ggsave("heatmap_all_variants_bottom_12226.eps",device = "eps",width=30,height=5)
ggsave("heatmap_all_variants_bottom_42226.pdf",device = "pdf",width=32,height=7)
ggsave("heatmap_all_variants_bottom_42226.eps",device = "eps",width=32,height=7)





##
mutation_colors <- c("Syn" = "#BAE0AD", "Mis" = "#6D91BA", "Non" = "#B3776F")

summary_stats <- function(y) {
  r <- quantile(y, probs = c(0.25, 0.5, 0.75))
  names(r) <- c("ymin", "y", "ymax")
  r
}
library(ggpubr)
for_means = compare_means(new_score_gdna_avg ~ Mutation.type,  data = non_wild_types,method = "t.test")
my_comparisons <- list( c("Syn", "Mis"), c("Syn", "Non"), c("Mis", "Non") )
# ggplot(non_wild_types, aes(x = Mutation.type, y = new_score_gdna_avg)) +
#   geom_jitter(aes(color = Mutation.type), width = 0.2, size = 3, alpha = 0.6) +  # Add transparency to jittered points
#   stat_summary(fun.data = summary_stats, geom = "errorbar",width=0.25, color = "black") +  # Increase size of error bars
#   #stat_summary(fun = median, geom = "point", size = 1, color = "black") +  # Increase size of median points
#   #stat_compare_means(comparisons = my_comparisons,method = "t.test") +
#   scale_color_manual(values = mutation_colors) + stat_compare_means(comparisons=my_comparisons,method = "t.test") + 
#   labs(x = "", y = "Score", title = "") +  # Add title
#   theme_minimal(base_size = 10) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
#     axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
#     axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 10)),
#     axis.text = element_text(size = 15),
#     legend.position = "none",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.line = element_line(colour = "black"),
#     axis.ticks = element_line(colour = "black")
#   ) +
#   geom_hline(yintercept = 0, linetype = "dotted") +  # Add horizontal line at y = 0
#   scale_y_continuous(breaks = seq(-12, 6, by = 1), limits = c(-12, 6))  # Adjust y-axis limits and breaks
# ggsave("donor_average_scores_82625.pdf",device = "pdf")
# ggsave("donor_average_scores_82625.eps",device = "eps")
# ggsave("donor_average_scores_82625.png",device = "png")

### donor compmarison
# ggplot(non_wild_types, aes(x = new_score_gdna_d1, y = new_score_gdna_d2)) +
#   geom_point(aes(color = Mutation.type),size=3) +
#   labs(x = "Donor 1 Score", y = "Donor 2 Score",caption="r=0.92") +
#   scale_color_manual(values = c("Mis" = "#1f77b4", "Non" = "#e63f29", "Syn" = "#42bd0d"))+
#   theme_minimal(base_size = 15) + xlim(-12,4)+ geom_smooth(method='lm', formula= y~x,se=F,fullrange=T)+
#   theme(legend.title = element_blank(),
#         legend.position = "right",
#         panel.grid = element_blank(), legend.text = element_text(size=20)# Remove gridlines
#   )+ theme(axis.text = element_text(size = 15),  # Adjust axis text size
#            axis.ticks = element_line(size = 0.5), 
#            panel.grid = element_blank(), axis.title.y.right = element_blank(), 
#            axis.title.x = element_text(size = 20, , margin = margin(t = 15)),
#            axis.title.y = element_text(size = 20, margin = margin(r = 15)),
#            aspect.ratio = .7, plot.title = element_text(hjust = 0.5), 
#            plot.caption = element_text(hjust = 0.5, size = 15, face = "italic", margin = margin(t = 10))) +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   geom_vline(xintercept = 0, linetype = "dotted")
# 
# ggsave("donor_comparisons_with_trendline_82625.pdf",device = "pdf")
# ggsave("donor_comparisons_with_trendline_82625.eps",device = "eps")
# ggsave("donor_comparisons_with_trendline_82625.png",device = "png")


###donor comparison v2
###colors are classifications, shapes are mutation types
mutation_colors <- c("Syn" = "#BAE0AD", "Mis" = "#6D91BA", "Non" = "#B3776F")
non_wild_types$classification = tools::toTitleCase(non_wild_types$classification)

donorcomp_plot = ggplot(non_wild_types, aes(x = new_score_gdna_d1, y = new_score_gdna_d2)) +
  geom_point(aes(color = classification,shape=Mutation.type),size=1.5) +
  labs(x = "Donor 1 Score", y = "Donor 2 Score",caption="r=0.92") +
  scale_color_manual(values = c("Functional" = "lightblue1", "Non-Functional" = "red", "Intermediate" = "grey"))+
  scale_shape_manual(values = c("Mis" = 15, "Non" = 16, "Syn" = 17)) + 
  theme_minimal(base_size = 15) + xlim(-12,4)+ geom_smooth(method='lm', formula= y~x,se=F,fullrange=T)+
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid = element_blank(), legend.text = element_text(size=20)# Remove gridlines
  )+ guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)),shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + theme(axis.text = element_text(size = 20,face="bold"),  # Adjust axis text size
                                                                                                                                                   axis.ticks = element_line(size = 0.5), 
                                                                                                                                                   panel.grid = element_blank(), axis.title.y.right = element_blank(), 
                                                                                                                                                   axis.title.x = element_text(size = 25, face="bold", margin = margin(t = 15)),
                                                                                                                                                   axis.title.y = element_text(size = 25,face="bold", margin = margin(r = 15)),
                                                                                                                                                   aspect.ratio = .7, plot.title = element_text(hjust = 0.5), legend.text = element_text(size=20,face="bold"),
                                                                                                                                                   plot.caption = element_text(hjust = 0.5, size = 20, face = "italic", margin = margin(t = 10))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")


donorcomp_legend = get_legend(donorcomp_plot)
as_ggplot(donorcomp_legend)
ggsave("fig1d_legend.pdf",device = "pdf")
ggsave("fig1d_legend.eps",device = "eps")


ggplot(non_wild_types, aes(x = new_score_gdna_d1, y = new_score_gdna_d2)) +
  geom_point(aes(color = classification,shape=Mutation.type),size=1.5) +
  labs(x = "Donor 1 Score", y = "Donor 2 Score",caption="r=0.92") +
  scale_color_manual(values = c("Functional" = "lightblue1", "Non-Functional" = "red", "Intermediate" = "grey"))+
  scale_shape_manual(values = c("Mis" = 15, "Non" = 16, "Syn" = 17)) + 
  theme_minimal(base_size = 15) + xlim(-12,4)+ geom_smooth(method='lm', formula= y~x,se=F,fullrange=T)+
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(), legend.text = element_text(size=20)# Remove gridlines
  )+ guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)),shape = guide_legend(override.aes = list(size = 4, alpha = 1))) + theme(axis.text = element_text(size = 20,face="bold"),  # Adjust axis text size
           axis.ticks = element_line(size = 0.5), 
           panel.grid = element_blank(), axis.title.y.right = element_blank(), 
           axis.title.x = element_text(size = 25, face="bold", margin = margin(t = 15)),
           axis.title.y = element_text(size = 25,face="bold", margin = margin(r = 15)),
           aspect.ratio = .7, plot.title = element_text(hjust = 0.5), legend.text = element_text(size=20,face="bold"),
           plot.caption = element_text(hjust = 0.5, size = 20, face = "italic", margin = margin(t = 10))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")

ggsave("fig1d_donor_comparisons_with_trendline_12226.pdf",device = "pdf",width=5,height=5)
ggsave("fig1d_donor_comparisons_with_trendline_12226.eps",device = "eps",width=5,height=5)
ggsave("fig1d_donor_comparisons_with_trendline_11226.eps",device = "png")



# gdna vs cdna
ggplot(non_wild_types, aes(x = new_score_gdna_avg, y = new_score_cdna_avg)) +
  geom_point(aes(color = Mutation.type),size=1.5) +
  labs(x = "gDNA Score", y = "cDNA Score",caption="r=0.91") +
  scale_color_manual(values = c("Mis" = "#1f77b4", "Non" = "#e63f29", "Syn" = "#42bd0d"))+
  theme_minimal(base_size = 15) + xlim(-12,4)+ geom_smooth(method='lm', formula= y~x,se=F,fullrange=T)+
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid = element_blank(), legend.text = element_text(size=20,face="bold")# Remove gridlines
  )+ guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) + theme(axis.text = element_text(size = 20, face="bold"),  # Adjust axis text size
           axis.ticks = element_line(size = 0.5), 
           panel.grid = element_blank(), axis.title.y.right = element_blank(), 
           axis.title.x = element_text(size = 25,face="bold" , margin = margin(t = 15)),
           axis.title.y = element_text(size = 25,face="bold", margin = margin(r = 15)),
           aspect.ratio = .7, plot.title = element_text(hjust = 0.5), 
           plot.caption = element_text(hjust = 0.5, size = 20, face = "italic", margin = margin(t = 10))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") + scale_y_continuous(breaks = seq(-12, 7, by = 4), limits = c(-12, 7))

ggsave("fig1e_cdna_gdna_comparisons_with_trendline_12226.pdf",device = "pdf",width=5,height=5)
ggsave("fig1e_cdna_gdna_comparisons_with_trendline_12226.eps",device = "eps",width=5,height=5)
ggsave("fig1e_cdna_gdna_comparisons_with_trendline_11226.png",device = "png")


### cutoffs
cutoff1 <-c(-1.545)
cutoff2 <- c(-0.796)
# ggplot(non_wild_types, aes(x = new_score_gdna_d1, y = new_score_gdna_d2)) +
#   geom_point(aes(color = Mutation.type),size=3) +
#   labs( x = "Donor 1 Score", y = "Donor 2 Score", title = (N[f] == 336 * " , " ~ N[i] == 23 * " , " ~ N[nf] == 130 * " , " ~ r == 0.917)) +
#   # Add the segments for the red cutoff
#   geom_segment(aes(x = cutoff1, y = min(new_score_gdna_d2), xend = cutoff1, yend = cutoff1), color = "darkred", linetype = "dashed", size = .9) +
#   geom_segment(aes(x = min(new_score_gdna_d1), y = cutoff1, xend = cutoff1, yend = cutoff1), color = "darkred", linetype = "dashed", size = .9) +
#   # Add the segments for the blue cutoff, pointing in the opposite direction
#   geom_segment(aes(x = cutoff2, y = max(new_score_gdna_d2), xend = cutoff2, yend = cutoff2), color = "blue", linetype = "dashed", size = .9) +
#   geom_segment(aes(x = 3.5, y = cutoff2, xend = cutoff2, yend = cutoff2), color = "blue", linetype = "dashed", size = .9) +
#   scale_color_manual(values = c("Mis" = "#1f77b4", "Non" = "#e63f29", "Syn" = "#42bd0d"))+
#   theme_minimal(base_size = 15) + xlim(-12,7)+
#   theme(legend.title = element_blank(),
#         legend.position = "right",
#         panel.grid = element_blank() # Remove gridlines
#   )+ theme(axis.text = element_text(size = 15),  # Adjust axis text size
#            axis.ticks = element_line(size = 0.5), 
#            panel.grid = element_blank(), axis.title.y.right = element_blank(), 
#            axis.title.x = element_text(size = 20, , margin = margin(t = 15)),
#            axis.title.y = element_text(size = 20, margin = margin(r = 15)),
#            aspect.ratio = .7, plot.title = element_text(hjust = 0.5), 
#            plot.caption = element_text(hjust = 0.5, size = 15, face = "italic", margin = margin(t = 10)))+
#   annotate("text", x = -8, y = -12, label = "Non-functional", 
#            hjust = 0, vjust = 0, size = 4, color = "darkred",fontface = "bold")+
#   annotate("text", x = -10, y = -1, label = "Intermediate", 
#            hjust = 0, vjust = 0, size = 4, color = "darkgrey",fontface = "bold")+
#   annotate("text", x = 2, y = 1.5, label = "Functional", 
#            hjust = 0, vjust = 0, size = 4, color = "blue",fontface = "bold") +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   geom_vline(xintercept = 0, linetype = "dotted")
# ggsave("donor_comparisons_with_cutoffs_82625.pdf",device="pdf")
# ggsave("donor_comparisons_with_cutoffs_82625.eps",device="eps")
# ggsave("donor_comparisons_with_cutoffs_82625.png",device="png")
# 

## I added new annos in this spreadsheet:
## dont need to do this jul 2025
#new_annos = read.xlsx("master_table_with_new_scores_and_classifications_and_clinvar_31025.xlsx")

#new_annos$ClinVar[which(new_annos$ClinVar=="X-SCID HMGD")] = "X-SCID"
#new_annos_no_wt = new_annos[-which(new_annos$Mutation.type=="WT_"),]
mutation_colors <- c(
  "Benign" = "#00ff00", 
  "Likely Benign" = "lightgreen", 
  "VUS" = "#C11C84", 
  "Likely Pathogenic" = "#de8bf9",  # A lighter red for Likely Pathogenic
  "X-SCID" = "#f94040", 
  "Unclassified" = "lightgrey") 
summary_stats <- function(y) {
  r <- quantile(y, probs = c(0.25, 0.5, 0.75))
  names(r) <- c("ymin", "y", "ymax")
  r
}
library(ggpubr)
#for_means = compare_means(new_score_gdna_avg ~ Mutation.type,  data = non_wild_types,method = "t.test")
my_comparisons <- list( c("Syn", "Mis"), c("Syn", "Non"), c("Mis", "Non") )

stat_df <- compare_means(
  new_score_gdna_avg ~ Mutation.type,
  data = non_wild_types,
  comparisons = my_comparisons
)

# Create formatted label
stat_df$label <- paste0("p=", signif(stat_df$p, 2))
stat_df$y.position <- c(6, 5, 4)
#

non_wild_types$ClinVar[which(non_wild_types$ClinVar=="X-SCID HMGD")] = "X-SCID"
non_wild_types$ClinVar[which(non_wild_types$ClinVar=="Likely Pathogenic HMGD")] = "Likely Pathogenic"
non_wild_types$ClinVar = factor(non_wild_types$ClinVar,levels = c("Benign","Likely Benign","Likely Pathogenic","X-SCID","VUS","Unclassified"))

mutplot = ggplot(non_wild_types, aes(x = Mutation.type, y = new_score_gdna_avg)) +
  geom_jitter(aes(color = ClinVar), width = 0.2, size = 4,alpha=0.8) +  # Add transparency to jittered points
  stat_summary(fun.data = summary_stats, geom = "errorbar", width = 0.1, size = 0.5, color = "black") +  # Increase size of error bars
  #stat_summary(fun = median, geom = "point", size = 1, color = "black") +  # Increase size of median points
  scale_color_manual(values = mutation_colors) +
  labs(x = "", y = "Score", title = "") +  # Add title
  theme_minimal(base_size = 7.5) + stat_pvalue_manual(
    stat_df,
    label = "label",
    tip.length = 0.01, size=6
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 20,face = "bold"), legend.text = element_text(size=15,face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),legend.title=element_blank()#,legend.justification = c(1, 1)
  ) + 
  geom_hline(yintercept = cutoff1, linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_hline(yintercept = cutoff2, linetype = "dashed", color = "blue") +  # Add horizontal line at y = 0
  scale_y_continuous(breaks = seq(-12, 7, by = 2), limits = c(-12, 7))  # Adjust y-axis limits and breaks

mutplot_legend = get_legend(mutplot)
as_ggplot(mutplot_legend)
ggsave("fig1c_legend.pdf",dev="pdf")
ggsave("fig1c_legend.eps",dev="eps")


ggplot(non_wild_types, aes(x = Mutation.type, y = new_score_gdna_avg)) +
  geom_jitter(aes(color = ClinVar), width = 0.2, size = 4,alpha=0.8) +  # Add transparency to jittered points
  stat_summary(fun.data = summary_stats, geom = "errorbar", width = 0.1, size = 0.5, color = "black") +  # Increase size of error bars
  #stat_summary(fun = median, geom = "point", size = 1, color = "black") +  # Increase size of median points
  scale_color_manual(values = mutation_colors) +
  labs(x = "", y = "Score", title = "") +  # Add title
  theme_minimal(base_size = 7.5) + stat_pvalue_manual(
    stat_df,
    label = "label",
    tip.length = 0.01, size=6
  ) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 20,face = "bold"), legend.text = element_text(size=15,face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),legend.title=element_blank(),legend.position="none"#legend.justification = c(1, 1)
  ) + 
  geom_hline(yintercept = cutoff1, linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = 0,linetype="dotted") +
  geom_hline(yintercept = cutoff2, linetype = "dashed", color = "blue") +  # Add horizontal line at y = 0
  scale_y_continuous(breaks = seq(-12, 7, by = 2), limits = c(-12, 7))  # Adjust y-axis limits and breaks

ggsave("fig_1c_score_with_clinvar_and_cutoffs_12226.pdf",device="pdf",height = 10,width=10)
ggsave("fig_1c_score_with_clinvar_and_cutoffs_12226.eps",device="eps",height=10,width=10)
ggsave("fig_1c_score_with_clinvar_and_cutoffs_12226.png",device="png")


#cadd
cadd = read.delim("CADD IL2RG EXON5 - Sheet1.tsv")
cadd_interesting = cadd[,c(1,2,3,4,5,6,7,8,9,17,18,24,25,29,36,37,38,39,115,134,135)]
cadd_interesting_exon5 = cadd_interesting[which(cadd_interesting$AnnoType=="CodingTranscript"),]
#cadd numbering is off?
cadd_interesting_exon5$cDNApos = cadd_interesting_exon5$cDNApos-54
#cadd seems to be reading a different strand so we must flip the bases...
#too lazy to test online ways so will doi t manually
cadd_interesting_exon5 = cadd_interesting_exon5[order(cadd_interesting_exon5$cDNApos),]
new_ref = c()
new_alt = c()

for (i in 1:nrow(cadd_interesting_exon5)) { 
  ref = cadd_interesting_exon5$Ref[i]
  if (ref=="A") {new_ref = c(new_ref,"T")}
  if (ref=="T") {new_ref = c(new_ref,"A")}
  if (ref=="C") {new_ref = c(new_ref,"G")}
  if (ref=="G") {new_ref = c(new_ref,"C")}
  alt = cadd_interesting_exon5$Alt[i]
  if (alt=="A") {new_alt = c(new_alt,"T")}
  if (alt=="T") {new_alt = c(new_alt,"A")}
  if (alt=="C") {new_alt = c(new_alt,"G")}
  if (alt=="G") {new_alt = c(new_alt,"C")}
  
}
cadd_interesting_exon5$Ref = new_ref
cadd_interesting_exon5$Alt = new_alt
cadd_interesting_exon5$MutcodeDNA = paste(cadd_interesting_exon5$Ref,cadd_interesting_exon5$cDNApos,cadd_interesting_exon5$Alt,sep="")

##
library(Biostrings)
wt = "GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGACAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGAAGCAATACTTCAAAAG"
#ekansh_loess_2 = read.csv("d1_d2_day3_loess_then_scale.csv")
#new_annos$sequences = og_annotation$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG
#seq = ekansh_loess_2$sequences[1]
#muts = pwalign::mismatchTable(pwalign::pairwiseAlignment(seq,wt))
#orig_base = muts$SubjectSubstring
#new_base = muts$PatternSubstring
mutcodes_dna = c()
library(pwalign)
for (seq in non_wild_types$GAACAATCAGTGGATTATAGACATAAGTTCTCCTTGCCTAGTGTGGATGGaCAGAAACGCTACACGTTTCGTGTTCGGAGCCGCTTTAACCCACTCTGTGGAAGTGCTCAGCATTGGAGTGAATGGAGCCACCCAATCCACTGGGGaAGCAATACTTCAAAAG){
  muts = pwalign::mismatchTable(pwalign::pairwiseAlignment(toupper(seq),wt))
  mutcode_curr = paste(muts$SubjectSubstring,muts$PatternStart+632,muts$PatternSubstring,sep="")
  mutcodes_dna = c(mutcodes_dna,mutcode_curr)
}
non_wild_types$MutcodeDNA = mutcodes_dna
screen_w_cadd = merge(non_wild_types,cadd_interesting_exon5,by="MutcodeDNA")
screen_w_cadd_classified = screen_w_cadd[-which(screen_w_cadd$ClinVar=="Unclassified"),]

am_with_nucleotides = read.delim("AlphaMissense_IL2RG_with_header.tsv")
am_with_nucleotides_exon5 = am_with_nucleotides[757:1121,]

new_ref = c()
new_alt = c()

for (i in 1:nrow(am_with_nucleotides_exon5)) { 
  ref = am_with_nucleotides_exon5$REF[i]
  if (ref=="A") {new_ref = c(new_ref,"T")}
  if (ref=="T") {new_ref = c(new_ref,"A")}
  if (ref=="C") {new_ref = c(new_ref,"G")}
  if (ref=="G") {new_ref = c(new_ref,"C")}
  alt = am_with_nucleotides_exon5$ALT[i]
  if (alt=="A") {new_alt = c(new_alt,"T")}
  if (alt=="T") {new_alt = c(new_alt,"A")}
  if (alt=="C") {new_alt = c(new_alt,"G")}
  if (alt=="G") {new_alt = c(new_alt,"C")}
  
}
am_with_nucleotides_exon5$REF = new_ref
am_with_nucleotides_exon5$ALT = new_alt
#am_with_nucleotides_exon5$MutcodeDNA = paste(cadd_interesting_exon5$REF,cadd_interesting_exon5$,cadd_interesting_exon5$ALT,sep="")

screen_w_cadd$pos_ref_alt = paste(screen_w_cadd$Pos,screen_w_cadd$Ref,screen_w_cadd$Alt,sep="_")

am_with_nucleotides_exon5$pos_ref_alt = paste(am_with_nucleotides_exon5$POS,am_with_nucleotides_exon5$REF,am_with_nucleotides_exon5$ALT,sep="_")

screen_w_cadd_and_am = merge(screen_w_cadd,am_with_nucleotides_exon5,by="pos_ref_alt")
#updated colors
mutation_colors <- c(
  "Benign" = "#00ff00", 
  "Likely Benign" = "lightgreen", 
  "VUS" = "#C11C84", 
  "Likely Pathogenic" = "#de8bf9",  # A lighter red for Likely Pathogenic
  "X-SCID" = "#f94040", 
  "Unclassified" = "lightgrey") 

screen_w_cadd$ClinVar = factor(screen_w_cadd$ClinVar,levels = c("Benign","Likely Benign","Likely Pathogenic","X-SCID","VUS","Unclassified"))

ggplot(screen_w_cadd,aes(x=new_score_gdna_avg,y=PHRED)) + geom_point(aes(color = ClinVar), size = 4, alpha = 0.8) +
  scale_color_manual(values = mutation_colors) +
  labs(x = "Screen Score", y = "CADD Score", title = "") +  # Add title
  theme_minimal(base_size = 7.5) +
  #geom_smooth(method = "loess", formula = y ~ x, color = "blue", se = FALSE,size=0.5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 20,face="bold"), legend.text = element_text(size=15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),legend.title=element_blank(), legend.position = "none"
  ) 
ggsave("fig2a_compare_to_cadd_12226.pdf",dev="pdf",width=5,height=5)
ggsave("fig2a_compare_to_cadd_12226.eps",dev="eps",width=5,height=5)
ggsave("fig2a_compare_to_cadd_12226.png",dev="png")


ggplot(screen_w_cadd,aes(x=new_score_gdna_avg,y=EsmScoreMissense)) + geom_point(aes(color = ClinVar), size = 4, alpha = 0.8) +
  scale_color_manual(values = mutation_colors) +
  labs(x = "Screen Score", y = "ESM1b Score", title = "") +  # Add title
  theme_minimal(base_size = 7.5) +
  #geom_smooth(method = "loess", formula = y ~ x, color = "blue", se = FALSE,size=0.5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 20,face="bold"), legend.text = element_text(size=15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),legend.title=element_blank(), legend.position = "none"
  ) 
ggsave("fig2a_compare_to_esm_12226.pdf",dev="pdf",width=5,height=5)
ggsave("fig2a_compare_to_esm_12226.eps",dev="eps",width=5,height=5)
ggsave("fig2a_compare_to_esm_11426.png",dev="png")



am_plot = ggplot(screen_w_cadd_and_am,aes(x=new_score_gdna_avg,y=am_pathogenicity)) + geom_point(aes(color = ClinVar), size = 4, alpha = 0.8) +
  scale_color_manual(values = mutation_colors) +
  labs(x = "Screen Score", y = "AlphaMissense Score", title = "") +  # Add title
  #geom_smooth(method = "loess", formula = y ~ x, color = "blue", se = FALSE,size=0.5) +
  theme_minimal(base_size = 7.5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 20,face="bold"),legend.text = element_text(size=20,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),legend.title=element_blank()
  ) 


am_legend = get_legend(am_plot)
as_ggplot(am_legend)
ggsave("fig2a_legend.pdf",dev="pdf")
ggsave("fig2a_legend.eps",dev="eps")

ggplot(screen_w_cadd_and_am,aes(x=new_score_gdna_avg,y=am_pathogenicity)) + geom_point(aes(color = ClinVar), size = 4, alpha = 0.8) +
  scale_color_manual(values = mutation_colors) +
  labs(x = "Screen Score", y = "AlphaMissense Score", title = "") +  # Add title
  #geom_smooth(method = "loess", formula = y ~ x, color = "blue", se = FALSE,size=0.5) +
  theme_minimal(base_size = 7.5) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 25, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 20,face="bold"),legend.text = element_text(size=15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),legend.position="none",legend.title=element_blank()
  ) 
ggsave("fig2a_compare_to_am_12226.pdf",dev="pdf",width=5,height=5)
ggsave("fig2a_compare_to_am_12226.eps",dev="eps",width=5,height=5)
ggsave("fig2a_compare_to_am_12226.png",dev="png",width=10,height=10)

### validation!
library(dplyr)
library(ggplot2)
library(tidyr)


## heatmap version of hsc data
master_with_variants_to_test = read.xlsx("master_table_mark_validation_82025.xlsx")
master_with_variants_to_test$validation[which(is.na(master_with_variants_to_test$validation))] = " "
master_with_variants_to_test_nwt = master_with_variants_to_test[which(master_with_variants_to_test$Mutation.type %in% c("Syn","Mis","Non")),]

match_to_master_validation_code  = c()
match_to_master_validation_num  = c()

for (m in non_wild_types$Mutcode) {
  validation_code = master_with_variants_to_test_nwt$validation[which(master_with_variants_to_test_nwt$Mutcode==m)]
  validation_num = master_with_variants_to_test_nwt$validation_num[which(master_with_variants_to_test_nwt$Mutcode==m)]
  match_to_master_validation_code = c(match_to_master_validation_code,validation_code)
  match_to_master_validation_num  = c(match_to_master_validation_num,validation_num)
}

non_wild_types$validation = match_to_master_validation_code
non_wild_types$validation_num = match_to_master_validation_num

mutation_colors <- c(
  "Benign" = "#00ff00", 
  "Likely Benign" = "lightgreen", 
  "Predicted functional" = "#f0ff00",
  "VUS" = "orange", 
  "Likely Pathogenic" = "#de8bf9",  # A lighter red for Likely Pathogenic
  "X-SCID" = "#f94040", 
  "Predicted non-functional" = "brown",
  " " = "darkgrey") 

non_wild_types$validation = factor(non_wild_types$validation,levels = c(" ","Benign","Likely Benign","Likely Pathogenic","Predicted functional","Predicted non-functional","VUS","X-SCID"))
# hack to make the grey ones get plotted before the colored ones
non_wild_types = non_wild_types[order(non_wild_types$validation),]
#tested_only = master_with_variants_to_test_nwt[which(!is.na(master_with_variants_to_test_nwt$validation)),]
# 
# ggplot(non_wild_types, aes(x = new_score_gdna_d1, y = new_score_gdna_d2)) +
#   geom_point(aes(color = validation),size=3) +
#   scale_color_manual(values = mutation_colors)+
#   labs(x = "Donor 1 Score", y = "Donor 2 Score",) +
#   theme_minimal(base_size = 15) + xlim(-12,4) +
#   theme(legend.title = element_blank(),
#         legend.position = "right",
#         panel.grid = element_blank(), legend.text = element_text(size=15)# Remove gridlines
#   )+ theme(axis.text = element_text(size = 7),  # Adjust axis text size
#            axis.ticks = element_line(size = 0.5), 
#            panel.grid = element_blank(), axis.title.y.right = element_blank(), 
#            axis.title.x = element_text(size = 15, , margin = margin(t = 15)),
#            axis.title.y = element_text(size = 15, margin = margin(r = 15)),
#            aspect.ratio = .7, plot.title = element_text(hjust = 0.5), 
#            plot.caption = element_text(hjust = 0.5, size = 15, face = "italic", margin = margin(t = 10))) +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   geom_vline(xintercept = 0, linetype = "dotted")
# 
# ggsave("validated_hsc_variants_71725.pdf",dev="pdf")

#master_with_variants_to_test = read.xlsx("master_table_mark_validation.xlsx")
#master_with_variants_to_test$validation_num[which(is.na(master_with_variants_to_test$validation_num))] = ""
#master_with_variants_to_test_nwt = master_with_variants_to_test[which(master_with_variants_to_test$Mutation.type %in% c("Syn","Mis","Non")),]
library(ggrepel)
ggplot(non_wild_types, aes(x = new_score_gdna_d1, y = new_score_gdna_d2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_point(aes(shape=Mutation.type),size=3) +
  #scale_color_manual(values = mutation_colors)+
  labs(x = "Donor 1 Score", y = "Donor 2 Score",) +
  geom_label_repel(size=8,data = subset(master_with_variants_to_test_nwt, !is.na(validation_num)),aes(label=validation_num), max.overlaps=Inf)+
  scale_shape_manual(values = c("Mis" = 15, "Non" = 16, "Syn" = 17)) + 
  theme_minimal(base_size = 15) + xlim(-12,4) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid = element_blank(), legend.text = element_text(size=20,face="bold")# Remove gridlines
  )+ theme(axis.text = element_text(size = 20,face="bold"),  # Adjust axis text size
           axis.ticks = element_line(size = 0.5), 
           panel.grid = element_blank(), axis.title.y.right = element_blank(), 
           axis.title.x = element_text(size = 25, , face="bold",margin = margin(t = 15)),
           axis.title.y = element_text(size = 25, face="bold",margin = margin(r = 15)),
           aspect.ratio = .7, plot.title = element_text(hjust = 0.5), 
           plot.caption = element_text(hjust = 0.5, size = 20, face = "italic", margin = margin(t = 10))) 

ggsave("fig2b_validated_hsc_variants_numbered_12226.pdf",dev="pdf",width=10,height=10)
ggsave("fig2b_validated_hsc_variants_numbered_12226.eps",dev="eps",width=10,height=10)



ggplot(non_wild_types, aes(x = new_score_gdna_d1, y = new_score_gdna_d2)) +
  geom_point(aes(color = Mutation.type),size=3) +
  geom_label_repel(data = subset(master_with_variants_to_test_nwt, !is.na(validation_num)),aes(label=validation_num), max.overlaps=Inf)+
  labs( x = "Donor 1 Score", y = "Donor 2 Score", title = (N[f] == 336 * " , " ~ N[i] == 19 * " , " ~ N[nf] == 131 * " , " ~ r == 0.917)) +
  # Add the segments for the red cutoff
  geom_segment(aes(x = cutoff1, y = min(new_score_gdna_d2), xend = cutoff1, yend = cutoff1), color = "darkred", linetype = "dashed", size = .9) +
  geom_segment(aes(x = min(new_score_gdna_d1), y = cutoff1, xend = cutoff1, yend = cutoff1), color = "darkred", linetype = "dashed", size = .9) +
  # Add the segments for the blue cutoff, pointing in the opposite direction
  geom_segment(aes(x = cutoff2, y = max(new_score_gdna_d2), xend = cutoff2, yend = cutoff2), color = "blue", linetype = "dashed", size = .9) +
  geom_segment(aes(x = 3.5, y = cutoff2, xend = cutoff2, yend = cutoff2), color = "blue", linetype = "dashed", size = .9) +
  scale_color_manual(values = c("Mis" = "#1f77b4", "Non" = "#e63f29", "Syn" = "#42bd0d"))+
  theme_minimal(base_size = 15) + xlim(-12,7)+
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid = element_blank() # Remove gridlines
  )+ theme(axis.text = element_text(size = 15),  # Adjust axis text size
           axis.ticks = element_line(size = 0.5), 
           panel.grid = element_blank(), axis.title.y.right = element_blank(), 
           axis.title.x = element_text(size = 20, , margin = margin(t = 15)),
           axis.title.y = element_text(size = 20, margin = margin(r = 15)),
           aspect.ratio = .7, plot.title = element_text(hjust = 0.5), 
           plot.caption = element_text(hjust = 0.5, size = 15, face = "italic", margin = margin(t = 10)))+
  annotate("text", x = -8, y = -12, label = "Non-functional", 
           hjust = 0, vjust = 0, size = 4, color = "darkred",fontface = "bold")+
  annotate("text", x = -10, y = -1, label = "Intermediate", 
           hjust = 0, vjust = 0, size = 4, color = "darkgrey",fontface = "bold")+
  annotate("text", x = 2, y = 1.5, label = "Functional", 
           hjust = 0, vjust = 0, size = 4, color = "blue",fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")



hsc_validation_round2 = read.table("all_sequence_counts_21525.txt",header=T)
hsc_validation_round1 = read.table("all_sequence_counts_10824.txt",header=T)
hsc_validation_round0 = read.table("all_sequence_counts_9524.txt",header=T)

hsc_validation_old_to_use = rbind(hsc_validation_round0[grep("HDRT",rownames(hsc_validation_round0)),],
                                  hsc_validation_round1[grep("Wk3",rownames(hsc_validation_round1)),])

variants_round2 = colnames(hsc_validation_round2)[-c(11,12,13)] # should remove H236H, WT, 2PAM
variants_round1 = colnames(hsc_validation_old_to_use)[-c(5,6,7)]
big_table = c()
for (v in variants_round2) {
  subset = hsc_validation_round2[grep(v,rownames(hsc_validation_round2)),c(v,"H236H")]
  init_ratio = subset[grep("HDRT",rownames(subset)),v]/subset[grep("HDRT",rownames(subset)),"H236H"]
  subset["ratio"] = subset[v]/subset$H236H/init_ratio
  subset_to_use = subset[grep('HDRT|DP|DN|B|_M|NK',rownames(subset)),]
  plotting_df = data.frame(category=sapply(strsplit(rownames(subset_to_use), "_"), "[", 2), value=subset_to_use$ratio)
  plotting_df$category[1] = "HDRT"
  plotting_df$variant = rep(v,nrow(plotting_df))
  big_table = rbind(big_table,plotting_df)
}

for (v in variants_round1) {
  subset = hsc_validation_old_to_use[grep(v,rownames(hsc_validation_old_to_use)),c(v,"H236H")]
  init_ratio = subset[grep("HDRT",rownames(subset)),v]/subset[grep("HDRT",rownames(subset)),"H236H"]
  subset["ratio"] = subset[v]/subset$H236H/init_ratio
  subset_to_use = subset[grep('HDRT|DP|DN|B|_M|NK',rownames(subset)),]
  plotting_df = data.frame(category=sapply(strsplit(rownames(subset_to_use), "_"), "[", 3), value=subset_to_use$ratio)
  plotting_df$category[1] = "HDRT"
  plotting_df$variant = rep(v,nrow(plotting_df))
  big_table = rbind(big_table,plotting_df)
}



big_table$category[which(big_table$category=="B")] = "CD19+"
big_table$category[which(big_table$category=="M")] = "CD33+"
big_table$category[which(big_table$category=="NK")] = "CD56+"
big_table$category[which(big_table$category=="DP")] = "CD5+/CD7+ T"
big_table$category[which(big_table$category=="DN")] = "CD5-/CD7- T"
big_table$category = factor(big_table$category,levels = c("HDRT","CD5+/CD7+ T",
                                                          "CD5-/CD7- T","CD19+","CD33+","CD56+"))
big_table$variant = factor(big_table$variant,levels = c("R226H","R222C","G232R",
                                                        "S238R","T220P","R222G","A234P","I244N",
                                                        "A234A","T220M","I244V","S233T",
                                                        "V223I","P211R"))
ggplot(big_table, aes(x = variant, y = category, fill = value)) +
  geom_tile() + 
  #scale_fill_gradient(low = "white", high = "black") +  # Grayscale gradient
  scale_fill_gradientn(
    colors = c("white", "black"),  # White to black gradient
    values = c(0, 1.25),       # Define stops in data range
    limits = c(0, 1.25),            # Clamp values at 1.5
    oob = scales::squish           # Anything above 1.5 stays black
  ) +
  theme_minimal() +
  labs(title = "",
       x = "Variant",
       y = "Cell Type",
       fill = "Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=15),
        axis.text.y = element_text(angle = 45, hjust = 1,size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text=element_text(size=10))
ggsave("hsc_differentiation_all_conditions_32825.pdf",dev="pdf",width=10)
ggsave("hsc_differentiation_all_conditions_32825.eps",dev="pdf",width=10)

big_table_for_paper = big_table[which(big_table$category %in% c("HDRT","CD5+/CD7+ T","CD19+")),]
ggplot(big_table_for_paper, aes(x = variant, y = category, fill = value)) +
  geom_tile() + 
  #scale_fill_gradient(low = "white", high = "black") +  # Grayscale gradient
  scale_fill_gradientn(
    colors = c("white", "black"),  # White to black gradient
    values = c(0, 1.25),       # Define stops in data range
    limits = c(0, 1.25),            # Clamp values at 1.5
    oob = scales::squish           # Anything above 1.5 stays black
  ) +
  
  theme_minimal() +
  labs(title = "",
       x = "Variant",
       y = "Cell Type",
       fill = "Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=15),
        axis.text.y = element_text(angle = 45, hjust = 1,size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text=element_text(size=10))
ggsave("hsc_differentiation_three_conditions_4325.pdf",dev="pdf",width=10)
ggsave("hsc_differentiation_three_conditions_4325.eps",dev="pdf",width=10)
ggsave("hsc_differentiation_three_conditions_4325.png",dev="png",width=10)
write.xlsx(big_table_for_paper,"hsc_validation_data_supp.xlsx")

big_table_for_paper$numbering = as.character(big_table_for_paper$variant)
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="T220P")] = 1
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="R222C")] = 2
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="R222G")] = 3
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="R226H")] = 4
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="A234P")] = 5
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="S238R")] = 6
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="I244N")] = 7
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="T220M")] = 8
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="S233T")] = 9
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="A234A")] = 10
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="I244V")] = 11
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="P211R")] = 12
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="V223I")] = 13
big_table_for_paper$numbering[which(big_table_for_paper$numbering=="G232R")] = 14
#big_table_for_paper$numbering = factor(big_table_for_paper$numbering,levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14))
big_table_for_paper$numbering = factor(big_table_for_paper$numbering,levels = c(14,13,12,11,10,9,8,7,6,5,4,3,2,1))
ggplot(big_table_for_paper, aes(x = category, y = numbering, fill = value)) +
  geom_tile() + 
  #scale_fill_gradient(low = "white", high = "black") +  # Grayscale gradient
  scale_fill_gradientn(
    colors = c("white", "black"),  # White to black gradient
    values = c(0, 1.25),       # Define stops in data range
    limits = c(0, 1.25),            # Clamp values at 1.5
    oob = scales::squish           # Anything above 1.5 stays black
  ) +
  
  theme_minimal() +
  labs(title = "",
       x = "Cell Type",
       y = "Variant",
       fill = "Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=25,face="bold"),
        axis.text.y = element_text(angle = 0, hjust = 1,size=25,face="bold"),
        axis.title.x = element_text(size = 25, face="bold"),
        axis.title.y = element_text(size = 25,face="bold"),
        legend.text=element_text(size=20,face="bold"),legend.title=element_text(size=20,face="bold"))
ggsave("fig2d_hsc_differentiation_three_conditions_12226.pdf",dev="pdf",width=5,height=10)
ggsave("fig2d_hsc_differentiation_three_conditions_12226.eps",dev="pdf",width=5,height=10)
ggsave("fig2d_hsc_differentiation_three_conditions_12226.png",dev="png",width=3)



## just the first screen iteration all conditions
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
ggplot(big_table_wk3_6, aes(x = variant, y = category, fill = value)) +
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
       x = "Variant",
       y = "Cell Type",
       fill = "Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=15),
        axis.text.y = element_text(angle = 45, hjust = 1,size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text=element_text(size=10))
ggsave("hsc_differentiation_first_donor_to_48_dp_4325.pdf",dev="pdf",width=10)
#ggsave("hsc_differentiation_all_conditions_32825.eps",dev="pdf",width=10)



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
big_table_t$variant = factor(big_table_t$variant,levels = c("F227C","G232R","R222P","R226C","R222G","R226H","S225R","W237*",
                                                            "W237R","W237S","A234P","E239*","R226P","R224G",
                                                            "R224P","W246S","R226G","A234A","G215E","H245D",
                                                            "H245L","H245H","H245N","I244I","I244V","K252R","S233T","G247R"))
ggplot(big_table_t, aes(x = variant, y = category, fill = values)) +
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
       x = "Variant",
       y = "Cell Type",
       fill = "Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=15),
        axis.text.y = element_text(angle = 45, hjust = 1,size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text=element_text(size=10))
ggsave("t_cell_differentiation_4325.pdf",dev="pdf",width=15)
#ggsave("t_cell_differentiation_31025.eps",dev="eps",width=15)


##shortened version

ratios = read.csv("~/Downloads/count files/ratios_test_benign_32_mutations.csv")
variants = ratios$aa_change
#G247E and G247A did not have a working hdrt value so I will exclude them
big_table_t= c()
for (v in variants){
  if (v %in% c("H236H","G247E","G247A")) {next}
  subset = ratios[which(ratios$aa_change==v),-c(1,2,3)]
  values = c(subset$ratio_hdrt,
             mean(c(subset$ratio_day3_d1,subset$ratio_day3_d2),na.rm=T),
             mean(c(subset$ratio_day5_il2_d1,subset$ratio_day5_il2_d2),na.rm=T))
  values = values/values[1]
  category = c("HDRT","Day 3","Day 5 IL2+")
  variant = rep(v,length(values)) 
  variant_out = data.frame(category,values,variant)
  big_table_t = rbind(big_table_t,variant_out)
}

big_table_t$category = factor(big_table_t$category,levels = c("HDRT","Day 3","Day 5 IL2+"))
big_table_t$variant = factor(big_table_t$variant,levels = c("F227C","G232R","R222P","R226C","R222G","R226H","S225R","W237*",
                                                            "W237R","W237S","G247G","A234P","E239*","R226P","R224G",
                                                            "R224P","W246S","R226G","A234A","G215E","H245D",
                                                            "H245L","H245H","H245N","I244I","I244V","K252R","S233T","G247R"))
ggplot(big_table_t, aes(x = variant, y = category, fill = values)) +
  geom_tile() + 
  scale_fill_gradientn(
    colors = c("white", "black"),  # White to black gradient
    values = c(0, 2.5, 5),       # Define stops in data range
    limits = c(0, 2.5),            # Clamp values at 1.5
    oob = scales::squish           # Anything above 1.5 stays black
  ) +  theme_minimal() +
  labs(title = "",
       x = "Variant",
       y = "Cell Type",
       fill = "Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=15),
        axis.text.y = element_text(angle = 45, hjust = 1,size=15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text=element_text(size=10))
ggsave("t_cell_differentiation_condensed_32825.pdf",dev="pdf",width=15)
#ggsave("t_cell_differentiation_31025.eps",dev="eps",width=15)
