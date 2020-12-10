# ABAMECTIN MANUSCRIPT SUPP #

###############################
# IMPORTANT!!!!!
# Set working directory
setwd("~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/")
###############################

library(tidyverse)
library(linkagemapping)
source("scripts/20201210_functions.R")

tsize = 12

#####################################
#       Fig S1 - dose 2017          #
#####################################


newdose <- read.csv("data/FileS1_newdose.csv")

# adjust for differences in control with subtraction
dosepruned_sub <- NULL
for(t in unique(newdose$trait)) {
    dosepruned2 <- newdose %>%
        dplyr::filter(trait == t) %>%
        dplyr::rename(dose = dose_uM)
    
    # get average control phenotype for each strain
    dosepruned3 <- dosepruned2 %>%
        dplyr::group_by(strain, dose, trait) %>%
        dplyr::filter(., dose == 0) %>%
        dplyr::mutate(mean_ctrl = mean(phenotype)) %>%
        dplyr::ungroup() %>%
        dplyr::select(condition, strain, mean_ctrl) %>%
        dplyr::full_join(dosepruned2, by = c("condition", "strain")) %>%
        dplyr::distinct() %>%
        dplyr::mutate(newpheno = phenotype - mean_ctrl)
    
    # add to final df
    dosepruned_sub <- rbind(dosepruned_sub, dosepruned3)
    
}

dosepruned_sub %>%
    dplyr::mutate(choose = ifelse(dose == 7.5, "*", NA)) %>%
    dplyr::group_by(condition, trait, dose) %>%
    dplyr::mutate(phen = max(newpheno, na.rm = T)) %>%
    dplyr::group_by(condition, trait, dose, strain) %>%
    dplyr::mutate(avg = mean(newpheno, na.rm = T),
                  sd = sd(newpheno, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dose = as.numeric(dose)) %>%
    tidyr::drop_na(condition) %>%
    dplyr::arrange(dose) %>%
    dplyr::distinct(condition, trait, dose, strain, .keep_all = T) %>%
    ggplot2::ggplot(.) +
    ggplot2::aes(x = dose, color = strain, group = interaction(strain, dose)) +
    ggplot2::geom_point(aes(y = avg)) +
    # ggplot2::geom_errorbar(aes(x = dose, ymin = avg - sd, ymax = avg + sd), width = 15) + # why is width so large? weird...
    ggplot2::geom_line(aes(y = avg, group = strain)) +
    ggplot2::geom_text(aes(y = phen, label = choose), color = "red", size = 6) +
    ggplot2::theme_bw(12) +
    ggplot2::scale_x_continuous(limits = c(0, 60),
                                breaks = unique(dosepruned_sub$dose),
                                expand=expand_scale(mult=c(0,0))) +
    ggplot2::theme(panel.grid = element_blank(),
                   strip.text = element_text(face = "bold"),
                   axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1),
                   axis.text.y = element_text(face = "bold", color = "black"),
                   axis.title = element_text(face = "bold", color = "black"),
                   legend.position = "bottom",
                   legend.text = element_text(face = "bold"),
                   legend.title = element_text(face = "bold")) +
    ggplot2::labs(x = "Dose (nM)", y = "Relative phenotype") +
    ggplot2::scale_color_manual(values = c("N2" = "orange", "CB4856" = "blue",
                                           "DL238" = "green", "JU775" = "purple"),
                                name = "Strain") +
    facet_wrap(~trait, scales = "free")
ggsave("figures/FigS1_dose2017.png", height = 5, width = 6)


#####################################
#       Fig S2 - WI/GWAS            #
#####################################

# make plots like linkage
wipheno <- read.csv("data/FileS2_wipheno.csv")
wimap <- read.csv("data/FileS3_wimap.csv")

# plot gwas for all four traits
i <- 1
plot_list <- list()
for(t in c("mean.TOF", "mean.EXT", "mean.norm.EXT", "norm.n")){
    
    ### plot WI phenos
    pheno <- wipheno %>%
        dplyr::filter(trait == t) %>%
        dplyr::group_by(strain) %>%
        dplyr::summarise(phenotype = mean(phenotype)) %>%
        dplyr::mutate(norm_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
        dplyr::arrange(norm_pheno) %>%
        na.omit() 
    
    pheno$strain <- factor(pheno$strain, levels = unique(pheno$strain))
    
    wi_pheno <- pheno %>%
        dplyr::mutate(strain_fill = dplyr::case_when(strain == "N2" ~ "N2",
                                                     strain == "CB4856" ~ "CB4856",
                                                     TRUE ~ "WT")) %>%
        dplyr::mutate(strain_fill = factor(strain_fill, levels = c("WT", "N2", "CB4856"))) %>%
        ggplot2::ggplot(.) +
        ggplot2::aes(x = strain, y = norm_pheno) +
        ggplot2::geom_bar(stat = "identity", fill = "grey", color = "grey") +
        theme_bw(tsize) +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              panel.grid = element_blank(),
              axis.title = element_text(face = "bold", color = "black"),
              axis.text.y = element_text(face = "bold", color = "black"),
              legend.position = "none",
              legend.title = element_blank()) +
        labs(x = "Strain", y = t)
    
    # GWAS mapping
    gwas_map <- wimap %>%
        dplyr::filter(trait == glue::glue("abamectin_{t}") & CHROM != "MtDNA") %>%
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
        ggplot2::geom_point( ggplot2::aes(color= factor(aboveEIGEN)) ) +
        ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
        ggplot2::theme_bw(tsize) +
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
    
    # pxg from gwas
    peaks <- wimap %>%
        dplyr::filter(trait == glue::glue("abamectin_{t}") & CHROM != "MtDNA") %>%
        na.omit()
    if(nrow(peaks) > 0) {
        if(t == "norm.n") {
            theme_size <- 9
        } else {
            theme_size <- tsize
        }
        pxg <- peaks %>% 
            dplyr::mutate(value=as.numeric(value),
                          allele = factor(allele, levels = c(-1,1), labels = c("REF","ALT"))) %>% 
            dplyr::mutate(norm_pheno = ((value - min(value, na.rm = T)) / (max(value, na.rm = T) - min(value, na.rm = T))),
                          value = norm_pheno) %>%
            dplyr::group_by(allele) %>%
            ggplot()+
            aes(x = allele, y = value, fill = allele)+
            geom_jitter(size = 0.5, width = 0.1) +
            geom_boxplot(alpha = 0.8, outlier.color = NA) +
            theme_bw(theme_size)+
            scale_fill_manual(values = c("REF" = "grey", "ALT" = "palevioletred1")) +
            ggplot2::facet_grid(~factor(marker, unique(peaks$marker)), scales = "free") + 
            ggplot2::theme_bw(tsize) + 
            ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", color = "black"), 
                           axis.text.y = ggplot2::element_text(face = "bold", color = "black"), 
                           axis.title.x = ggplot2::element_text(face = "bold", color = "black", vjust = -0.3),
                           axis.title.y = ggplot2::element_text(face = "bold", color = "black"), 
                           strip.text.x = ggplot2::element_text(face = "bold", color = "black"),
                           strip.text.y = ggplot2::element_text(face = "bold", color = "black"),
                           plot.title = ggplot2::element_blank(), 
                           legend.position = "none", 
                           panel.grid = element_blank(),
                           panel.background = ggplot2::element_rect(color = NA, size = 0.6)) +
            ggplot2::labs(x = "Genotype at QTL", y = t)
        
        # cowplot
        plot_list[[i]] <- (cowplot::plot_grid(wi_pheno, gwas_map, pxg, nrow = 3, align = "v", axis = "lr", labels = c("A", "B", "C")))
    } else {
        pxg <- ggplot2::ggplot() +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.border = element_rect(color = "grey99", size = 0.0001)) + ### THIS IS RIDICULOUS
            ggplot2::labs(x = " ", y = " ")
        # cowplot
        plot_list[[i]] <- (cowplot::plot_grid(wi_pheno, gwas_map, pxg, nrow = 3, align = "v", axis = "lr", labels = c("A", "B", "")))
    }
    
    # next plot
    i <- i + 1
}

ggsave(plot_list[[1]], file = "figures/FigS2-1_gwas.png", width = 7.5, height = 8)
ggsave(plot_list[[2]], file = "figures/FigS2-2_gwas.png", width = 7.5, height = 8)
ggsave(plot_list[[3]], file = "figures/FigS2-3_gwas.png", width = 7.5, height = 8)
ggsave(plot_list[[4]], file = "figures/FigS2-4_gwas.png", width = 7.5, height = 8)



#####################################
#       Fig S3 - RIAIL/LM           #
#####################################


riailpheno <- read.csv("data/FileS4_riailpheno.csv")
linkagemapping::load_cross_obj("N2xCB4856cross_full")
riailmap <- read.csv("data/FileS5_riailmap.csv")

i <- 1
plot_list <- list()
for(t in c("mean.TOF", "mean.EXT", "mean.norm.EXT", "norm.n")){
    
    ### plot riail phenos
    pheno <- riailpheno %>%
        dplyr::filter(trait == t) %>%
        dplyr::mutate(strain_fill = dplyr::case_when(strain == "N2" ~ "n2",
                                                     strain == "CB4856" ~ "cb",
                                                     TRUE ~ "RIL")) %>%
        dplyr::group_by(strain, trait) %>%
        dplyr::mutate(avg_phen = mean(phenotype, na.rm = T)) %>%
        dplyr::distinct(strain, trait, .keep_all = T) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(norm_pheno = ((avg_phen - min(avg_phen, na.rm = T)) / (max(avg_phen, na.rm = T) - min(avg_phen, na.rm = T)))) %>%
        dplyr::arrange(norm_pheno)
    
    pheno$strain <- factor(pheno$strain, levels = unique(pheno$strain))
    
    riail <- pheno %>%
        ggplot2::ggplot(.) +
        ggplot2::aes(x = strain, y = norm_pheno, fill = strain_fill, color = strain_fill) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = c("n2" = "orange", "cb" = "blue", "RIL" = "grey")) +
        ggplot2::scale_color_manual(values = c("n2" = "orange", "cb" = "blue", "RIL" = "grey")) +
        theme_bw(tsize) +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              panel.grid = element_blank(),
              axis.title = element_text(face = "bold", color = "black"),
              axis.text.y = element_text(face = "bold", color = "black"),
              legend.position = "none") +
        labs(x = "Strain", y = t)
    
    ### plot lod
    
    # plot LOD
    test <- riailmap %>%
        dplyr::mutate(condtrt = trait,
                      trait = stringr::str_split_fixed(trait, "_", 2)[,2]) %>%
        dplyr::filter(trait == t)
    lod <- plot_lods(test)

    peak <- test %>%
        na.omit()
    
    if(nrow(peak) > 0) {
        # relative riail pheno
        rph <- riailpheno %>% 
            dplyr::filter(trait == t) %>%
            dplyr::mutate(rel_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
            dplyr::mutate(phenotype = rel_pheno)
        
        # drugcross
        drugcross <- linkagemapping::mergepheno(N2xCB4856cross_full2, rph, set = 2)
        
        test2 <- test %>%
            dplyr::mutate(condition = stringr::str_split_fixed(condtrt, "_", 2)[,1])
        
        # plot pxg
        pxg <- plot_pxg(drugcross, test2, yaxis = t)
        
        # cowplot
        plot_list[[i]] <- (cowplot::plot_grid(riail, lod, pxg, nrow = 3, align = "v", axis = "lr", labels = c("A", "B", "C")))
    } else {
        pxg <- ggplot2::ggplot() +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.border = element_rect(color = "grey99", size = 0.0001)) + ### THIS IS RIDICULOUS
            ggplot2::labs(x = " ", y = " ")
        # cowplot
        plot_list[[i]] <- (cowplot::plot_grid(riail, lod, pxg, nrow = 3, align = "v", axis = "lr", labels = c("A", "B", "")))
    }
    
    # next plot
    i <- i + 1
    
}

ggsave(plot_list[[1]], file = "figures/FigS3-1_linkage.png", width = 7.5, height = 8)
ggsave(plot_list[[2]], file = "figures/FigS3-2_linkage.png", width = 7.5, height = 8)
ggsave(plot_list[[3]], file = "figures/FigS3-3_linkage.png", width = 7.5, height = 8)
ggsave(plot_list[[4]], file = "figures/FigS3-4_linkage.png", width = 7.5, height = 8)

#####################################
#   Fig S4 - mapping overview       #
#####################################

riailmap <- read.csv("data/FileS5_riailmap.csv")
wimap <- read.csv("data/FileS3_wimap.csv")

# get chromosome lengths
chr_lens <- data.frame(chr = c("I", "II", "III", "IV", "V", "X"), 
                       start = rep(1,6), 
                       end = c(14972282,15173999,13829314,17450860,20914693,17748731),
                       condition = "abamectin",
                       trait = "mean.TOF",
                       method = "Linkage")

# convert linkage and gwa to same columns
new_linkage <- riailmap %>%
    na.omit() %>%
    dplyr::select(chr, pos, ci_l_pos, ci_r_pos, trait, sig = lod) %>%
    dplyr::mutate(method = "Linkage")

new_gwas <- wimap %>%
    na.omit() %>%
    dplyr::select(chr = CHROM, pos = POS, ci_l_pos = startPOS, ci_r_pos = endPOS, trait, sig = log10p) %>%
    dplyr::distinct() %>%
    dplyr::mutate(method = "Association")

newmap <- new_linkage %>%
    dplyr::bind_rows(new_gwas) %>%
    tidyr::separate(trait, c("condition", "trait"), sep = "_")


ggplot(newmap)+
    aes(x=pos/1E6, y=trait)+
    theme_bw(tsize) +
    viridis::scale_fill_viridis(name = "Significance") + 
    viridis::scale_color_viridis(name = "Significance") +
    geom_segment(aes(x = ci_l_pos/1e6, y = trait, xend = ci_r_pos/1e6, yend = trait, color = sig), size = 2, alpha = 1) +
    geom_segment(data=chr_lens, aes(x =  start/1e6, xend = end/1e6, y = trait, yend=trait), color='transparent', size =0.1) +
    geom_point(aes(fill = sig), color = "black", shape = 21, size = 2, alpha = 1)+
    # scale_shape_manual(values = c("yes" = 24, "no" = 25)) +
    xlab("Genomic position (Mb)") + ylab("") +
    guides(shape = FALSE) +
    theme(axis.text.x = element_text(size=10, face="bold", color="black"),
          # axis.ticks.y = element_blank(),
          legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10),
          legend.key.size = unit(.75, "cm"),
          panel.grid.major.x = element_line(),
          panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_text(size = 10, face = "bold", color = "black"),
          axis.title.x = element_text(size=12, face="bold", color= "black"),
          axis.title.y = element_blank(),
          strip.text.x = element_text(size=12, face="bold", color="black"),
          strip.text.y = element_text(size=12, face="bold", color="black", angle = 0),
          strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
          plot.title = element_text(size=12, face="bold")) +
    facet_grid(factor(method, levels = c("Linkage", "Association")) ~ chr, scales = "free_x", space = "free")
ggsave("figures/FigS4_mappingsummary.png", width = 7.5, height = 3)



#####################################
#         Fig S5 - scan2            #
#####################################

riailpheno <- read.csv("data/FileS4_riailpheno.csv")
linkagemapping::load_cross_obj("N2xCB4856cross_full")
riailmap <- read.csv("data/FileS5_riailmap.csv")

# set trait
trt <- "mean.EXT"

# filter for trait
drugtrait <- riailpheno %>%
    dplyr::filter(!strain %in% c("N2", "CB4856")) %>%
    dplyr::filter(trait == trt, condition == "abamectin") %>%
    dplyr::select(strain, condition, trait, phenotype)

# merge cross object
drugcross <- linkagemapping::mergepheno(N2xCB4856cross_full2, drugtrait, set = 2)

# run scan2
scan2 <- qtl::scantwo(drugcross, pheno.col=3, method="mr")

# riail marker conversion
mappos <- qtl::pull.map(N2xCB4856cross_full2, as.table = TRUE) %>%
    dplyr::mutate(marker = rownames(.),
                  cM = as.character(pos)) %>%
    dplyr::select(-pos, -chr) %>%
    dplyr::distinct(cM, .keep_all = T)  %>%
    dplyr::mutate(pos = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2]))

# change markers in scan2 from cM to genomic position
scan2$map$pos <- mappos$pos

# plot
png("figures/FigS5_scantwo.png")
plot(scan2) 
dev.off()


#####################################
#          Fig S6 - CSSV            #
#####################################

trt <- "mean.TOF"

# load data
cssVpruned <- read.csv("data/FileS10_cssV.csv")
HTAstats <- read.csv("data/FileS11_HTAstats.csv")

# all
strainset <- c("N2", "CB4856", "ECA573", "ECA554")

# regress
regressed <- easysorter::regress(cssVpruned %>%
                                     dplyr::filter(strain %in% strainset))%>%
    dplyr::mutate(strain_fill = dplyr::case_when(strain %in% c("N2") ~ "N2",
                                                 strain %in% c("CB4856") ~ "CB",
                                                 TRUE ~ "NIL")) %>%
    dplyr::filter(trait == trt) %>%
    dplyr::mutate(rel_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::mutate(phenotype = rel_pheno)

# add stats
stats <- HTAstats %>%
    dplyr::filter(exp == "cssV",
                  comparison %in% c("N2-ECA573", "ECA554-CB4856"),
                  condition == "abamectin",
                  trait == trt)%>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(strain = dplyr::case_when(comparison == "N2-ECA573" ~ "ECA573",
                                            comparison == "ECA554-CB4856" ~ "ECA554",
                                            TRUE ~ "NA"),
                  groups = dplyr::case_when(grepl("CB4856", comparison) ~ "CB",
                                            grepl("N2", comparison) ~ "N2",
                                            TRUE ~ "NIL"),
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"))

plots <- plot_nil(regressed, nil_genotypes, stats, strainset, "V", ylab = "Mean length")

cowplot::plot_grid(plots[[1]],
                   plots[[2]],
                   plots[[3]],
                   nrow = 1, ncol = 3, rel_widths = c(3, 1, 4.5), align = "h", axis = "b", labels = c("A", "", "B"))

# save
ggsave("figures/FigS6_cssV.png", height = 4, width = 7.5)


#####################################
#  Fig S7 -fine map with NIL        #
#####################################

variants_chrV <- read.csv("data/FileS7_variants_chrV.csv")
wimap <- read.csv("data/FileS3_wimap.csv")
riailmap <- read.csv("data/FileS5_riailmap.csv")
nil_genotypes <- read.csv("data/FileS9_nilgenotypes.csv")


# peaks
gwaspeaks <- wimap %>%
    dplyr::filter(trait == "abamectin_mean.TOF") %>%
    dplyr::distinct(startPOS, endPOS)

lm_peaks <- riailmap %>%
    na.omit() %>%
    dplyr::filter(trait == "abamectin_mean.EXT", chr == "V") %>%
    dplyr::mutate(peak_marker = dplyr::case_when(pos < 3e6 ~ "VL",
                                                 pos > 10e6 ~ "VR",
                                                 TRUE ~ "VC"))

glc_test <- variants_chrV %>%
    dplyr::filter(gene_name == "WBGene00001591",
                  strain == "CB4856") %>%
    dplyr::filter(log10p == max(log10p))%>%
    dplyr::mutate(peak_marker = "VR")

lgc_test <- variants_chrV %>%
    dplyr::filter(gene_name == "WBGene00020528",
                  strain == "CB4856") %>%
    dplyr::filter(log10p == max(log10p))

glc3_test <- variants_chrV %>%
    dplyr::filter(gene_name == "WBGene00001593",
                  strain == "CB4856") %>%
    dplyr::filter(log10p == max(log10p))

# plot all variants, color CB variants
chrV_variant_plot <- variants_chrV %>%
    dplyr::mutate(gt = ifelse(allele == REF, "REF", "ALT"),
                  # peak = ifelse(POS %in% gwaspeaks$peakPOS, T, F), 
                  chr = "chrV") %>%
    tidyr::drop_na(gt) %>%
    ggplot(.) +
    aes(x = POS/1e6, y = log10p, fill = gt, color = gt) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("REF" = "grey", "ALT" = "blue"), name = "CB4856 genotype") +
    scale_fill_manual(values = c("REF" = "grey", "ALT" = "blue"), name = "CB4856 genotype") +
    theme_bw(tsize) +
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
    theme(axis.text.y = element_text(color = "black", face = "bold"),
          axis.title.y = element_text(color = "black", face = "bold"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(color = "black", face = "bold"),
          panel.grid = element_blank(),
          legend.position = "none") +
    # show peak marker from gwas
    geom_point(data = glc_test, aes(x = POS/1e6, y = log10p), fill = "red", size = 3, shape = 23, color = "black") +
    # add lgc-54
    geom_point(data = lgc_test, aes(x = POS/1e6, y = log10p), fill = "red", size = 3, shape = 22, color = "black") +
    # add glc-3
    geom_point(data = glc3_test, aes(x = POS/1e6, y = log10p), fill = "red", size = 3, shape = 21, color = "black") +
    facet_grid(~CHROM, scales = "free", space = "free") +
    geom_vline(xintercept = c(0, 3.120167)) +
    geom_vline(xintercept = c(5.260997, 5.906132), linetype = "dashed") +
    geom_vline(xintercept = c(13.678801, 19.303558), linetype = "dotted")

# NIL plot
nils <- nil_plot(rev(c("ECA1059", "ECA1065", "ECA377", "ECA232", "ECA632")), "V", order = F)[[1]] +
    theme(strip.background = element_blank(),
          strip.text = element_blank()) +    
    geom_vline(xintercept = c(0, 3.120167)) +
    geom_vline(xintercept = c(5.260997, 5.906132), linetype = "dashed") +
    geom_vline(xintercept = c(13.678801, 19.303558), linetype = "dotted")

plot <- cowplot::plot_grid(chrV_variant_plot, nils, nrow = 2, align = "v", axis = "l", labels = c("A", "B"), rel_heights = c(2,1))
ggsave(plot, filename = "figures/FigS7_finemap_NIL.png", height = 5, width = 7.5)

#####################################
#       Fig S8 - mediation          #
#####################################

aba_mediation <- read.csv("data/FileS15_mediation.csv")

# PLOT
aba_mediation %>%
    tidyr::separate(eqtl_marker, c("chrom", "pos"), "_") %>%
    dplyr::mutate(pos = as.numeric(pos),
                  qtl = ifelse(pheno_marker == "V_2895960", "VL", "VC"),
                  gene = dplyr::case_when(probe == "A_12_P115868" ~ "ugt-53",
                                          probe == "A_12_P106882" ~ "nuo-5",
                                          probe == "A_12_P102287" ~ "H23N18.4")) %>%
    dplyr::filter(var == "med") %>%
    dplyr::arrange(pos) %>%
    dplyr::group_by(pheno_marker) %>%
    dplyr::mutate(q90 = quantile(abs_est, probs = 0.9)[[1]]) %>%
    ggplot(.) +
    aes(x = pos/1e6, y = abs_est, text = probe) +
    geom_point(aes(alpha = pval), color = "black") +
    scale_alpha_continuous(range = c(1, 0.1)) +
    theme_bw(tsize) +
    labs(x = "Genomic position (Mb)", y = "Mediation estimate") +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          axis.text = element_text(face="bold", color="black"),
          axis.title = element_text(face="bold", color="black"),
          strip.text = element_text(face = "bold", color = "black")) +
    geom_hline(aes(yintercept = q90), color = "grey") +
    facet_grid(~ factor(qtl, levels = c("VL", "VC")), scales = "free", space = "free") +
    geom_text(aes(label = gene, x = pos/1e6 + 0.4, y = abs_est + 0.01), fontface = "italic")
ggsave("figures/FigS8_mediation.png", height = 4, width = 7.5)


#####################################
#       Fig S9 - grace NIL          #
#####################################

nil_genotypes <- read.csv("data/FileS9_nilgenotypes.csv")
VC_NIL <- read.csv("data/FileS16_VCNIL.csv")
HTAstats <- read.csv("data/FileS11_HTAstats.csv")
trt <- "mean.EXT"

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

plots <- plot_nil(regressed, nil_genotypes, stats, strainset, "V", ylab = "Mean optical density")

# 4446729-7374928 - NIL defined interval
cowplot::plot_grid(plots[[1]] +
                       geom_vline(xintercept = c(4.446729,7.374928), linetype = "dashed"),
                   plots[[2]],
                   plots[[3]], 
                   nrow = 1, ncol = 3, rel_widths = c(3, 1, 4.5), align = "h", axis = "b", labels = c("A", "", "B"))

# save
ggsave("figures/FigS9_VCNIL.png", height = 4, width = 7.5)

