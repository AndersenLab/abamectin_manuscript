# ABAMECTIN MANUSCRIPT FIGURES #

###############################
# IMPORTANT!!!!!
# Set working directory
setwd("~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/")
###############################

library(tidyverse)
library(linkagemapping)
source("scripts/20201210_functions.R")
tsize <- 12
options(stringsAsFactors=FALSE)


#####################################
#       Fig 1 - mappings            #
#####################################
t <- "mean.EXT"

# load gwer linkage map
riailmap <- read.csv("data/FileS5_riailmap.csv")

# drug cross
linkagemapping::load_cross_obj("N2xCB4856cross_full")

# load pheno
riailpheno <- read.csv("data/FileS4_riailpheno.csv")

# plot linkage
lod <- plot_lods(riailmap %>%
                     dplyr::filter(trait == glue::glue("abamectin_{t}")) %>%
                     dplyr::mutate(condtrt = trait)) 

# plot gwas
wimap <- read.csv("data/FileS3_wimap.csv")
trt <- "mean.TOF"

# GWAS mapping
gwas_map <- wimap %>%
    dplyr::filter(trait == glue::glue("abamectin_{trt}") & CHROM != "MtDNA") %>%
    dplyr::distinct(marker, .keep_all = T) %>%
    ggplot2::ggplot(.) +
    ggplot2::aes(x = POS/1e6, y = log10p) +
    ggplot2::scale_color_manual(values = c("0" = "black", "1" = "deeppink2", "2" = "red")) +
    ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6,    # this is the plot boundary for LD and gene plots
                                    xmax = endPOS/1e6,    # this is the plot boundary for LD and gene plots
                                    ymin = 0, 
                                    ymax = Inf, 
                                    fill = "palevioletred1"), 
                       color = "palevioletred1",fill = "palevioletred1",linetype = 2, alpha=.3)+
    ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN), color = "gray", alpha = .75, size = 1) +
    ggplot2::geom_point(ggplot2::aes(color= factor(aboveEIGEN)), size = 0.7) +
    ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
    ggplot2::theme_bw(12) +
    ggplot2::theme(
        axis.text = ggplot2::element_text(color = "black", face = "bold"),
        axis.title = ggplot2::element_text(face = "bold", color = "black"),
        strip.text = ggplot2::element_text(face = "bold", color = "black"),
        plot.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.position = "none",
        panel.background = ggplot2::element_rect(color = NA, size = 0.6)) +
    ggplot2::labs(x = "Genomic position (Mb)",
                  y = expression(-log[10](italic(p))))


gwaspeaks <- wimap %>%
    dplyr::filter(trait == "abamectin_mean.TOF") %>%
    dplyr::distinct(startPOS, endPOS)


# fine map
variants_chrV <- read.csv("data/FileS7_variants_chrV.csv")

# linkage mapping CI
lm_peaks <- riailmap %>%
    na.omit() %>%
    dplyr::filter(trait == "abamectin_mean.EXT", chr == "V") %>%
    dplyr::mutate(peak_marker = dplyr::case_when(pos < 3e6 ~ "VL",
                                                 pos > 10e6 ~ "VR",
                                                 TRUE ~ "VC"))

# plot where glc1 variants are
glc_test <- variants_chrV %>%
    dplyr::filter(gene_name == "WBGene00001591",
                  strain == "CB4856") %>%
    dplyr::filter(log10p == max(log10p))%>%
    dplyr::mutate(peak_marker = "VR")

# gwas data
gwaspeaks <- wimap %>%
    na.omit() %>%
    dplyr::filter(trait == "abamectin_mean.TOF") %>%
    dplyr::distinct(startPOS, peakPOS, endPOS, log10p, EIGEN) %>%
    dplyr::mutate(yval = dplyr::case_when(startPOS < 3e6 ~ 8,
                                          startPOS > 10e6 ~ 6,
                                          TRUE ~ 5.6),
                  peak_marker = dplyr::case_when(startPOS < 3e6 ~ "VL",
                                                 startPOS > 10e6 ~ "VR",
                                                 TRUE ~ "VC"))

# plot all variants, color CB variants
finemap <- variants_chrV %>%
    dplyr::mutate(gt = ifelse(allele == REF, "REF", "ALT"),
                  peak = ifelse(POS %in% gwaspeaks$peakPOS, T, F), 
                  chr = "chrV") %>%
    tidyr::drop_na(gt) %>%
    ggplot(.) +
    aes(x = POS/1e6, y = log10p, fill = gt, color = gt) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("REF" = "grey", "ALT" = "blue"), name = "CB4856 genotype") +
    scale_fill_manual(values = c("REF" = "grey", "ALT" = "blue"), name = "CB4856 genotype") +
    theme_bw(tsize) +
    geom_vline(data = gwaspeaks, aes(xintercept = startPOS/1e6), linetype = "dashed", color = "palevioletred1") +
    geom_vline(data = gwaspeaks, aes(xintercept = endPOS/1e6), linetype = "dashed", color = "palevioletred1") +
    geom_hline(data = gwaspeaks, aes(yintercept = EIGEN), color = "grey", size = 0.5) +
    geom_rect(data = gwaspeaks, aes(xmin = startPOS/1e6, xmax = endPOS/1e6, ymin = -Inf, ymax = Inf), fill = "palevioletred1", alpha = 0.3, inherit.aes = F) +
    geom_rect(data = lm_peaks, aes(xmin = ci_l_pos/1e6, xmax = ci_r_pos/1e6, ymin = -Inf, ymax = Inf), fill = "skyblue", alpha = 0.3, inherit.aes = F) +
    # show QTL VL, VC, and VR
    geom_segment(aes(x = 1.757246, xend = 4.333001, y = 8, yend = 8), color = "black", size = 1) +
    geom_segment(aes(x = 6.118360, xend = 7.342129, y = 8, yend = 8), color = "black", size = 1) +
    geom_segment(aes(x = 15.983112, xend = 16.599066, y = 8, yend = 8), color = "black", size = 1) +
    geom_text(aes(x = 3.045124, y = 8.5), label = "VL", color = "black", fontface = "bold") +
    geom_text(aes(x = 6.730244, y = 8.5), label = "VC", color = "black", fontface = "bold") +
    geom_text(aes(x = 16.29109, y = 8.5), label = "VR", color = "black", fontface = "bold") +
    labs(x = "Genomic position (Mb)", y = expression(-log[10](italic(p)))) +
    theme(axis.text = element_text(color = "black", face = "bold"),
          axis.title = element_text(color = "black", face = "bold"),
          strip.text = element_text(color = "black", face = "bold"),
          panel.grid = element_blank(),
          legend.position = "none") +
    # show peak marker from gwas
    geom_point(data = glc_test, aes(x = POS/1e6, y = log10p), fill = "red", size = 2, shape = 23, color = "black") +
    facet_grid(~CHROM, scales = "free", space = "free")

cowplot::plot_grid(gwas_map, lod, finemap, nrow = 3, labels = c("A", "B", "C"), rel_heights = c(1, 1, 1.5))
ggsave("figures/Fig1_QTLmaps.png", width = 7.5, height = 8)

#####################################
#       Fig 2 - NILpheno            #
#####################################

# load data
NILpruned <- read.csv("data/FileS12_NILpheno.csv")
riailmap <- read.csv("data/FileS5_riailmap.csv")
nil_genotypes <- read.csv("data/FileS9_nilgenotypes.csv")
HTAstats <- read.csv("data/FileS11_HTAstats.csv")

trt <- "mean.EXT"

# get peaks
peaks <- riailmap %>%
    na.omit() %>%
    dplyr::filter(trait == glue::glue("abamectin_{trt}"),
                  chr == "V") %>%
    dplyr::mutate(pos = pos / 1e6)


# all
strainset <- c("N2", "CB4856", "ECA1059", "ECA1065", "ECA377", "ECA232")

# regress
regressed <- easysorter::regress(NILpruned %>%
                                     dplyr::filter(strain %in% strainset))%>%
    dplyr::mutate(strain_fill = dplyr::case_when(strain %in% c("N2") ~ "N2",
                                                 strain %in% c("CB4856") ~ "CB",
                                                 TRUE ~ "NIL"),
                  strain_geno = dplyr::case_when(strain %in% c("N2", "ECA232") ~ "NNN",
                                                 strain == "CB4856" ~ "CCC",
                                                 strain %in% c("ECA1059") ~ "CCN",
                                                 strain == "ECA1066" ~ "NCC",
                                                 strain == "ECA1065" ~ "CNN",
                                                 strain == "ECA377" ~ "NCN")) %>%
    dplyr::filter(trait == trt) %>%
    dplyr::mutate(rel_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::mutate(phenotype = rel_pheno)

# add stats
stats <- HTAstats %>%
    dplyr::filter(exp == "NILpruned",
                  grepl("N2", comparison),
                  comparison != "N2-CB4856",
                  trait == trt)%>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(strain = stringr::str_split_fixed(comparison, "N2-", 2)[,2],
                  groups = dplyr::case_when(grepl("CB4856", comparison) ~ "CB",
                                            grepl("N2", comparison) ~ "N2",
                                            TRUE ~ "NIL"),
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"))

plots <- plot_nil(regressed, nil_genotypes, stats, strainset, "V", ylab = "Mean optical density")

cowplot::plot_grid(plots[[1]] +
                       geom_vline(data = peaks, aes(xintercept = pos), color = "red"),
                   plots[[2]],
                   plots[[3]] +
                       geom_text(aes(x = strain, y = phen + 0.15, label = strain_geno), size = tsize/4), 
                   nrow = 1, ncol = 3, rel_widths = c(3, 1, 4.5), align = "h", axis = "b", labels = c("A", "", "B"))

# save
ggsave("figures/Fig2_nilpheno.png", width = 7.5, height = 4)

#####################################
#       Fig 3 - grace NIL           #
#####################################

# load data
nil_genotypes <- read.csv("data/FileS9_nilgenotypes.csv")
VC_NIL <- read.csv("data/FileS14_VCNIL.csv")
HTAstats <- read.csv("data/FileS11_HTAstats.csv")

trt <- "mean.TOF"

# strains
strainset <- c("N2", "CB4856", "ECA629", "ECA634", "ECA636", "ECA632")

# regress
regressed <- easysorter::regress(VC_NIL %>%
                                     dplyr::filter(strain %in% strainset))%>%
    dplyr::mutate(strain_fill = dplyr::case_when(strain %in% c("N2") ~ "N2",
                                                 strain %in% c("CB4856") ~ "CB",
                                                 TRUE ~ "NIL")) %>%
    dplyr::filter(trait == trt) %>%
    dplyr::mutate(rel_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::mutate(phenotype = rel_pheno)

# add stats
stats <- HTAstats %>%
    dplyr::filter(exp == "VC_NIL",
                  grepl("N2", comparison),
                  comparison != "N2-CB4856",
                  trait == trt)%>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(strain = stringr::str_split_fixed(comparison, "N2-", 2)[,2],
                  groups = dplyr::case_when(grepl("CB4856", comparison) ~ "CB",
                                            grepl("N2", comparison) ~ "N2",
                                            TRUE ~ "NIL"),
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"))

plots <- plot_nil(regressed, nil_genotypes, stats, strainset, "V", ylab = "Mean length")

# 4446729-7374928 - NIL defined interval
cowplot::plot_grid(plots[[1]] +
                       geom_vline(xintercept = c(4.446729,7.374928), linetype = "dashed"),
                   plots[[2]],
                   plots[[3]], 
                   nrow = 1, ncol = 3, rel_widths = c(3, 1, 4.5), align = "h", axis = "b", labels = c("A", "", "B"))

# save
ggsave("figures/Fig3_VCNIL.png", height = 4, width = 7.5)


#####################################
#       Fig 4 - lgc-54 del          #
##################################### 

# load data
lgc54del <- read.csv("data/FileS17_lgc54del.csv")
nil_genotypes <- read.csv("data/FileS9_nilgenotypes.csv")
HTAstats <- read.csv("data/FileS11_HTAstats.csv")

trt <- "mean.EXT"

# genome
# fake FX genos
fx1 <- nil_genotypes %>%
    dplyr::filter(sample == "N2") %>%
    dplyr::mutate(sample = "FX3518")

# fake FX geno
fx2 <- nil_genotypes %>%
    dplyr::filter(sample == "N2") %>%
    dplyr::mutate(sample = "FX3448")

# add to nil genotypes
nilgeno <- nil_genotypes %>%
    dplyr::filter(sample %in% c("N2", "CB4856")) %>%
    dplyr::bind_rows(fx1, fx2)

# plot chrV genotsype
chrVgeno <- nilgeno %>%
    dplyr::filter(chrom == "V") %>%
    dplyr::mutate(chrom = "chrV",
                  sample = factor(sample, levels = rev(c("N2", "CB4856", "FX3448", "FX3518")))) %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2), alpha = 0.7)+
    facet_grid(~chrom, scales = "free",  space = "free")+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_bw(tsize) +
    theme(
        axis.text.x = element_text(face="bold", color="black"),
        axis.text.y = element_text(face="bold", color="black"),
        axis.title.x = element_text(face="bold", color="black"),
        axis.title.y = element_text(face="bold", color="black"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face="bold"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
    labs(x = "Genomic position (Mb)", y = "") +
    geom_point(data = nilgeno %>%
                   dplyr::filter(grepl("FX", sample), chrom == "V") %>%
                   dplyr::mutate(chrom = "chrV"),
               aes(x = 6.8, y = sample), shape = 24, color = "black", fill = "grey", size = 5)

# genome plot
genome_plot <- nilgeno %>%
    dplyr::filter(chrom == "V") %>%
    dplyr::mutate(chrom = "Genome",
                  sample = factor(sample, levels = rev(c("N2", "CB4856", "FX3448", "FX3518")))) %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2), alpha = 0.7)+
    facet_grid(~chrom, scales = "free",  space = "free")+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_bw(tsize) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())+
    labs(x = "Genomic position (Mb)", y = "")


# first just aba
# regress
aba <- easysorter::regress(lgc54del %>%
                               dplyr::filter(strain %in% c("N2", "CB4856", "FX3518", "FX3448"))) %>%
    dplyr::select(date:phenotype) %>%
    # dplyr::mutate(condition = "abamectin") %>%
    dplyr::filter(trait == trt) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(relative_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::mutate(phenotype = relative_pheno) %>%
    dplyr::mutate(strain_fill = dplyr::case_when(strain %in% c("N2") ~ "N2",
                                                 strain %in% c("CB4856") ~ "CB",
                                                 TRUE ~ "NIL"),
                  groups = dplyr::case_when(strain == "CB4856" ~ "CB",
                                            TRUE ~ "N2"))

# add stats
stats <- HTAstats %>%
    dplyr::filter(exp == "lgc54del",
                  condition == "abamectin",
                  grepl("N2", comparison),
                  comparison != "N2-CB4856",
                  trait == trt)%>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(strain = dplyr::case_when(comparison == "N2-FX3448" ~ "FX3448",
                                            comparison == "N2-FX3518" ~ "FX3518",
                                            TRUE ~ "NA"),
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"))

# add stats to dataframe
test <- aba %>%
    dplyr::left_join(stats) %>%
    dplyr::group_by(strain, condition) %>%
    dplyr::mutate(phen = max(phenotype, na.rm = T) + 0.1)

test$strain <- factor(test$strain, levels = rev(c("N2", "CB4856", "FX3448", "FX3518")), ordered = T)

pheno <- test %>%
    ggplot2::ggplot(.) +
    ggplot2::aes(x = strain, y = phenotype, fill = strain_fill) +
    ggplot2::geom_jitter(width = 0.1, size = 0.5, alpha = 0.8) +
    ggplot2::geom_boxplot(outlier.color = NA, alpha = 0.5, size = 0.3) +
    ggplot2::geom_text(aes(label = sig, y = phen, color = groups), size = tsize/3, angle = -90) +
    ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB" = "blue", "NIL" = "grey")) +
    ggplot2::scale_color_manual(values = c("N2" = "orange", "CB" = "blue")) +
    ggplot2::theme_bw(tsize) +
    ggplot2::theme(axis.text.x = element_text(face = "bold", color = "black"),
                   axis.title.x = element_text(face = "bold", color = "black"),
                   strip.text = element_text(face = "bold"),
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   legend.position = "none",
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank()) +
    ggplot2::labs(x = "", y = "Mean optical density") +
    ggplot2::facet_grid(~condition, scales = "free") +
    ggplot2::coord_flip() +
    scale_y_continuous(limits = c(0, 1.1),
                       breaks = seq(0,1, 0.25))

abaplot <- cowplot::plot_grid(chrVgeno, genome_plot, pheno, nrow = 1, align = "h", axis = "b", rel_widths = c(1, 0.3, 1))

# show just DMSO
dmso <- lgc54del %>%
    dplyr::filter(strain %in% c("N2", "CB4856", "FX3518", "FX3448"),
                  condition == "DMSO") %>%
    dplyr::mutate(condition = ifelse(condition == "DMSO", "DMSO control", condition)) %>%
    dplyr::select(date:phenotype) %>%
    dplyr::filter(trait == trt) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(relative_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::mutate(phenotype = relative_pheno) %>%
    dplyr::mutate(strain_fill = dplyr::case_when(strain %in% c("N2") ~ "N2",
                                                 strain %in% c("CB4856") ~ "CB",
                                                 TRUE ~ "NIL"),
                  groups = dplyr::case_when(strain == "CB4856" ~ "CB",
                                            TRUE ~ "N2"))

# add stats
stats <- HTAstats %>%
    dplyr::filter(exp == "lgc54del",
                  condition == "DMSO",
                  grepl("N2", comparison),
                  comparison != "N2-CB4856",
                  trait == trt)%>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(strain = dplyr::case_when(comparison == "N2-FX3448" ~ "FX3448",
                                            comparison == "N2-FX3518" ~ "FX3518",
                                            TRUE ~ "NA"),
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"),
                  condition = "DMSO control")

# add stats to dataframe
test <- dmso %>%
    dplyr::left_join(stats) %>%
    dplyr::group_by(strain, condition) %>%
    dplyr::mutate(phen = max(phenotype, na.rm = T) + 0.1)

test$strain <- factor(test$strain, levels = rev(c("N2", "CB4856", "FX3448", "FX3518")), ordered = T)

pheno <- test %>%
    ggplot2::ggplot(.) +
    ggplot2::aes(x = strain, y = phenotype, fill = strain_fill) +
    ggplot2::geom_jitter(width = 0.1, size = 0.5, alpha = 0.8) +
    ggplot2::geom_boxplot(outlier.color = NA, alpha = 0.5, size = 0.3) +
    ggplot2::geom_text(aes(label = sig, y = phen, color = groups), size = tsize/3, angle = 90) +
    ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB" = "blue", "NIL" = "grey")) +
    ggplot2::scale_color_manual(values = c("N2" = "orange", "CB" = "blue")) +
    ggplot2::theme_bw(tsize) +
    ggplot2::theme(axis.text.x = element_text(face = "bold", color = "black"),
                   axis.title.x = element_text(face = "bold", color = "black"),
                   strip.text = element_text(face = "bold"),
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   legend.position = "none",
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank()) +
    ggplot2::labs(x = "", y = "Mean optical density") +
    ggplot2::facet_grid(~condition, scales = "free") +
    ggplot2::coord_flip() +
    scale_y_continuous(limits = c(0, 1.1),
                       breaks = seq(0,1, 0.25))

dmsoplot <- cowplot::plot_grid(chrVgeno, genome_plot, pheno, nrow = 1, align = "h", axis = "b", rel_widths = c(1, 0.3, 1))

cowplot::plot_grid(abaplot, dmsoplot, nrow = 2, align = "v", axis = "l", labels = c("A", "B"))

ggsave("figures/Fig4_lgc54.png", height = 6, width = 7.5)


#####################################
#   Fig 5 - HC overlap plot         #
#####################################

# chr lengths to define plot paramters
HCON_chrV_length = 48868368
CELE_chrV_length = 20924180
difference = (HCON_chrV_length - CELE_chrV_length)/2

# read in table and force HCON start/stop coordinates to be numeric
new_table <- read.csv("data/FileS18_hc_overlap.csv")

filtered_table <- new_table %>%
    tidyr::drop_na(HCON_start) %>%
    dplyr::filter(is.na(VR_nil)) %>%
    dplyr::mutate(CELE_position = (start + stop)/2) %>%
    dplyr::mutate(HCON_position = (HCON_start + HCON_stop)/2) %>%
    dplyr::mutate(HCON_interval = if_else(interval_1_status == "Y" | interval_2_status == "Y", "Y", "N")) %>%
    dplyr::select(transcriptID, CELE_position, HCON_ortho, HCON_position, HCON_interval)

# colors for figure
ce_color <- "#e07403"
hc_color <- "#397a4c"
line_color <- "black"

ggplot() + 
    geom_segment(data=filtered_table %>%
                     dplyr::filter(is.na(HCON_interval)), 
                 aes(y = 1, yend=10, x=HCON_position, xend=CELE_position+difference),
                 size = 0.25, color = "grey") + 
    geom_segment(data=filtered_table %>%
                     dplyr::filter(HCON_interval == "Y"), 
                 aes(y = 1, yend=10, x=HCON_position, xend=CELE_position+difference),
                 size = 0.25, color = line_color) + 
    theme_bw(12) +
    # header and footer boxes
    geom_rect(aes(ymin=10, ymax=10.7, xmin=difference, xmax=CELE_chrV_length+difference), fill="#fcbc7a", colour=NA) + # grey90
    geom_rect(aes(ymin=0.3, ymax=1, xmin=0, xmax=HCON_chrV_length), fill="#bedb92", colour=NA) + # grey60
    # VL region
    geom_rect(aes(ymin=10, ymax=10.7, xmin=1+difference, xmax=4333001+difference), fill=ce_color, alpha = 0.8,colour=NA, size=0.3) +
    geom_text(aes(x=difference+2166500, y=10.35), label= "VL", color = "black", size = 3) +
    # VC region
    geom_rect(aes(ymin=10, ymax=10.7, xmin=5260997+difference, xmax=7342129+difference), fill=ce_color, alpha = 0.8,colour=NA, size=0.3) +
    geom_text(aes(x=difference+6301563, y=10.35), label= "VC", color = "black" , size = 3) +
    # header and footer text
    # geom_text(aes(x=difference+(CELE_chrV_length/2), y=10.3), label=expression(~bolditalic('C. elegans')~bold(''))) +
    # geom_text(aes(x=difference+(CELE_chrV_length/2), y=0.7), label=expression(~bolditalic('H. contortus')~bold(''))) +
    geom_text(aes(x=difference+(CELE_chrV_length*0.7), y=10.35), label=expression(~bolditalic('C. elegans')~bold('chrV'))) +
    geom_text(aes(x=(CELE_chrV_length)-difference, y=0.65), label=expression(~bolditalic('H. contortus')~bold('chrV'))) + ## CAN  YOU SEE TEXT IF NOT WHITE?????
    # VR region
    # geom_rect(aes(ymin=10, ymax=10.6, xmin=15933659+difference, xmax=16336743+difference), fill="skyblue", alpha = 0.3,colour="black", size=0.3) +
    # position of H. contortus intervals
    geom_rect(aes(ymin=0.3, ymax=1, xmin=37000000, xmax=42000000), fill=hc_color, color=NA, alpha = 0.8) +
    geom_rect(aes(ymin=0.3, ymax=1, xmin=45000000, xmax=48000000), fill=hc_color, color=NA, alpha = 0.8) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          axis.title.x = element_text(face = "bold")) +
    # header and footer outline
    geom_rect(aes(ymin=10, ymax=10.7, xmin=difference, xmax=CELE_chrV_length+difference), fill=NA, colour="black", size = 0.3) +
    geom_rect(aes(ymin=0.3, ymax=1, xmin=0, xmax=HCON_chrV_length), fill=NA, colour="black", size = 0.3) +
    # add tick lines?
    geom_rect(aes(ymin = 0.2, ymax = 0.3, xmin = 0, xmax = 0), color = "black") +
    geom_text(aes(y = 0.01, x = 0), label = "0", size = 3) +
    geom_rect(aes(ymin = 0.2, ymax = 0.3, xmin = 5e6, xmax = 5e6), color = "black") +
    geom_text(aes(y = 0.01, x = 5e6), label = "5", size = 3) +
    geom_rect(aes(ymin = 0.2, ymax = 0.3, xmin = 10e6, xmax = 10e6), color = "black") +
    geom_rect(aes(ymin = 0.2, ymax = 0.3, xmin = 15e6, xmax = 15e6), color = "black") +
    geom_rect(aes(ymin = 0.2, ymax = 0.3, xmin = 20e6, xmax = 20e6), color = "black") +
    geom_rect(aes(ymin = 0.2, ymax = 0.3, xmin = 25e6, xmax = 25e6), color = "black") +
    geom_rect(aes(ymin = 0.2, ymax = 0.3, xmin = 30e6, xmax = 30e6), color = "black") +
    geom_rect(aes(ymin = 0.2, ymax = 0.3, xmin = 35e6, xmax = 35e6), color = "black") +
    geom_rect(aes(ymin = 0.2, ymax = 0.3, xmin = 40e6, xmax = 40e6), color = "black") +
    geom_rect(aes(ymin = 0.2, ymax = 0.3, xmin = 45e6, xmax = 45e6), color = "black") +
    geom_text(aes(y = 0.01, x = 10e6), label = "10", size = 3) +
    geom_text(aes(y = 0.01, x = 15e6), label = "15", size = 3) +
    geom_text(aes(y = 0.01, x = 20e6), label = "20", size = 3) +
    geom_text(aes(y = 0.01, x = 25e6), label = "25", size = 3) +
    geom_text(aes(y = 0.01, x = 35e6), label = "35", size = 3) +
    geom_text(aes(y = 0.01, x = 40e6), label = "40", size = 3) +
    geom_text(aes(y = 0.01, x = 45e6), label = "45", size = 3) +
    geom_text(aes(y = 0.01, x = 30e6), label = "30", size = 3) +
    # add c.e. tick lines
    geom_rect(aes(ymin = 10.7, ymax = 10.8, xmin = difference, xmax = difference), color = "black") +
    geom_rect(aes(ymin = 10.7, ymax = 10.8, xmin = 5e6+difference, xmax = 5e6+difference), color = "black") +
    geom_rect(aes(ymin = 10.7, ymax = 10.8, xmin = 10e6+difference, xmax = 10e6+difference), color = "black") +
    geom_rect(aes(ymin = 10.7, ymax = 10.8, xmin = 15e6+difference, xmax = 15e6+difference), color = "black") +
    geom_rect(aes(ymin = 10.7, ymax = 10.8, xmin = 20e6+difference, xmax = 20e6+difference), color = "black") +
    geom_text(aes(y = 11, x = difference), label = "0", size = 3) +
    geom_text(aes(y = 11, x = 5e6+difference), label = "5", size = 3) +
    geom_text(aes(y = 11, x = 10e6+difference), label = "10", size = 3) +
    geom_text(aes(y = 11, x = 15e6+difference), label = "15", size = 3) +
    geom_text(aes(y = 11, x = 20e6+difference), label = "20", size = 3) +
    labs(x = "Genomic position (Mb)")

ggsave("figures/Fig5_hc_overlap.png", height = 3.5, width = 7.5)

