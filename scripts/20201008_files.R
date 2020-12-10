# ABAMECTIN MANUSCRIPT FILES #

library(tidyverse)
library(linkagemapping)

#####################################
#       File S1 - dose 2017         #
#####################################

# dose response from briana 2017
newdose <- easysorter::read_data("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/dose/brianadose/") %>%
    easysorter::remove_contamination() %>%
    easysorter::sumplate(directories = FALSE, quantiles = TRUE) %>%
    easysorter::bioprune() %>%
    tidyr::gather(trait, phenotype, -c(date:col))%>%
    dplyr::ungroup() %>%
    dplyr::filter(trait %in% c("mean.TOF", "mean.EXT", "norm.n", "mean.norm.EXT")) %>%
    dplyr::mutate(control = ifelse(condition == "DMSO", "None", control)) %>%
    dplyr::mutate(dose = paste0("0.", stringr::str_split_fixed(condition, "abamectin",2)[,2]),
                  dose = ifelse(condition == "DMSO", "0", dose),
                  dose_uM = as.numeric(dose) * 1000, 
                  condition = "abamectin",
                  date = "20170501",
                  experiment = "dose9a",
                  round = 9,
                  assay = "a") %>%
    dplyr::select(date:condition, dose_uM, control:phenotype)
save(newdose, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS1_newdose.Rda")


#####################################
#       File S2 - WI pheno          #
#####################################

# # DON'T LIKE THAT THERE IS NO METADATA ASSOCIATED HERE... DATE, ROW, COL, BLEACH, ETC.
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/GWAS/abamectin_pheno.Rda")

wipheno <- abamectin_pheno %>%
    dplyr::filter(trait %in% c("norm.n", "mean.TOF", "mean.EXT", "mean.norm.EXT"))
save(wipheno, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS2_wipheno.Rda")

vcf_samples <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/GWAS/vcf_samples.txt") %>%
    dplyr::select(-(`#CHROM`:FORMAT))
strains <- names(vcf_samples)

length(unique(wipheno$strain))

# to cegwas
wipheno2 <- wipheno %>%
    dplyr::filter(strain %in% strains) %>%
    dplyr::mutate(trait = paste0(condition, "_", trait)) %>%
    dplyr::select(-condition) %>%
    tidyr::spread(trait, phenotype)

length(unique(wipheno2$strain))

test <- wipheno2 %>%
    na.omit()
length(unique(test$strain))

readr::write_tsv(wipheno2, "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS2_wipheno.tsv")


#####################################
#       File S3 - WI map            #
#####################################

# prep data for cegwas2
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS2_wipheno.Rda")

# average value per strain
pheno <- wipheno %>%
    dplyr::group_by(trait, strain) %>%
    dplyr::summarize(value = mean(phenotype)) %>%
    tidyr::spread(trait, value)

readr::write_tsv(pheno, "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/GWAS/20200728_gwas_pheno.tsv")

# RUN CEGWAS

# new data 20201008
meanEXT <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/GWAS/Analysis_Results-20201006/Mappings/Data/abamectin_mean.EXT_processed_mapping.tsv")
meanTOF <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/GWAS/Analysis_Results-20201006/Mappings/Data/abamectin_mean.TOF_processed_mapping.tsv")
normn <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/GWAS/Analysis_Results-20201006/Mappings/Data/abamectin_norm.n_processed_mapping.tsv")
meannormEXT <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/GWAS/Analysis_Results-20201006/Mappings/Data/abamectin_mean.norm.EXT_processed_mapping.tsv")

wimap <- meanEXT %>%
    dplyr::bind_rows(meanTOF, normn, meannormEXT) %>%
    dplyr::rename(EIGEN = BF, aboveEIGEN = aboveBF)

save(wimap, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS3_wimap.Rda")


#####################################
#       File S4 - riailpheno        #
#####################################
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/linkage/allRIAILsregressed.Rda")
data("N2xCB4856cross")

sets <- N2xCB4856cross$pheno %>%
    dplyr::filter(set == 2)

riailpheno <- allRIAILsregressed %>%
    dplyr::filter(condition == "abamectin",
                  trait %in% c("norm.n", "mean.TOF", "mean.EXT", "mean.norm.EXT"),
                  strain %in% c("N2", "CB4856", as.character(sets$strain))) %>%
    dplyr::select(date:phenotype)
save(riailpheno, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS4_riailpheno.Rda")


#####################################
#       File S5 - LM                #
#####################################
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/linkage/abamectin-GWER.proximal.annotated.Rda")

riailmap <- annotatedmap %>%
    dplyr::mutate(condition = stringr::str_split_fixed(trait, "\\.", 2)[,1],
                  trait = stringr::str_split_fixed(trait, "\\.", 2)[,2])%>%
    dplyr::filter(trait %in% c("norm.n", "mean.TOF", "mean.EXT", "mean.norm.EXT")) %>%
    tidyr::unite(trait, c(condition, trait), sep = "_")
save(riailmap, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS5_riailmap.Rda")

#####################################
#     File S6 - scan2summary        #
#####################################

# read in data from quest
scan2summary <- readr::read_tsv("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/linkage/20200728_scan2summary.tsv")
save(scan2summary, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS6_scan2summary.Rda")

#####################################
#       File S7 - chrV var          #
#####################################

# OG script: 20200727_finemapQTL_cegwas2.R

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS5_riailmap.Rda")

# peaks <- riailmap %>%
#     na.omit() %>%
#     dplyr::filter(trait == "abamectin_mean.EXT", chr == "V")
peaks <- data.frame(chr = c("V", "V", "V"),
                    ci_l_pos = c(2182093,4983265,15933659),
                    ci_r_pos = c(3698359,7342129,17863725),
                    marker = c("VL", "VC", "VR"),
                    peak_marker = c("VL", "VC", "VR"))

# load processed mapping df
pr_map <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/GWAS/Analysis_Results-20201006/Mappings/Data/abamectin_mean.TOF_processed_mapping.tsv")

map_pheno <- na.omit(pr_map) %>%
    dplyr::distinct(strain, value) %>%
    as.data.frame()

map_pheno$value <- as.numeric(map_pheno$value)

# load genotype matrix
complete_matrix <- readr::read_tsv("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/GWAS/Analysis_Results-20201006/Genotype_Matrix/Genotype_Matrix.tsv") %>%
    na.omit()


chrV_matrix <- complete_matrix %>%
    dplyr::filter(CHROM == "V")

kinship_matrix <- rrBLUP::A.mat(t(complete_matrix[,5:ncol(complete_matrix)]), 
                                n.core = 4)

kinship_matrix <- kinship_matrix[row.names(kinship_matrix)%in%map_pheno$strain,
                                 colnames(kinship_matrix)%in%map_pheno$strain]

gwa_mapping <- function (data, 
                         cores = 4, 
                         kin_matrix = kinship_matrix, 
                         snpset = NULL, 
                         min.MAF = 0.05,
                         p3d = FALSE) {
    x <- data
    
    y <- snpset %>% dplyr::mutate(marker = paste0(CHROM, "_", POS)) %>% 
        dplyr::select(marker, everything(), -REF, -ALT) %>% 
        as.data.frame()
    
    kin <- as.matrix(kin_matrix)
    
    pmap <- rrBLUP::GWAS(pheno = x, 
                         geno = y, 
                         K = kin, 
                         min.MAF = min.MAF, 
                         n.core = cores, 
                         P3D = as.logical(p3d), 
                         plot = FALSE)
    
    return(pmap)
}

# perform fine mapping for entire chromosome, 1 mb at a time
allmaps <- NULL
for(i in 1:21) {
    
    # look 1 Mb at a time
    ROI_matrix <- chrV_matrix %>%
        dplyr::filter(POS < i*1e6, POS > (i-1)*1e6)
    
    ROI_matrix <- ROI_matrix[colnames(ROI_matrix)%in%c("CHROM","POS","REF","ALT",map_pheno$strain)]
    
    roi_mapping <- gwa_mapping(data = map_pheno,
                               snpset = ROI_matrix,
                               kin_matrix = kinship_matrix)
    
    allmaps <- rbind(allmaps, roi_mapping)
    
}

allmaps <- allmaps %>%
    dplyr::rename(log10p = value)
# save(allmaps, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/allmaps.Rda")


# without running ALL snpeff annotations, just need to know which vars for CB?
# use new VCF - hard filtered
q_vcf <- glue::glue("http://storage.googleapis.com/elegansvariation.org/releases/20200815/variation/WI.20200815.hard-filter.vcf.gz")

# snpeff annotations
snpeff_out <- list()

# Now add snpeff annotations for all variants on chromosome V
for(i in 1:21) {
    
    snpeff_out[[i]] <- cegwas2::query_vcf(glue::glue("V:{i-1}000000-{i}000000"),
                                          impact = "ALL",
                                          vcf = q_vcf,
                                          samples = "CB4856") %>%
        tidyr::unite(marker, CHROM, POS, sep = "_") %>%
        dplyr::left_join(allmaps, ., by = "marker") %>% 
        tidyr::drop_na(REF)
    
}

# bind all rows together
snpeff_df <- dplyr::bind_rows(snpeff_out)

variants_chrV <- snpeff_df %>%
    dplyr::select(marker, CHROM:ALT, strain = SAMPLE, allele = a1, effect:transcript_biotype, nt_change:aa_change)%>%
    dplyr::distinct(marker, CHROM, POS, strain, REF, ALT, .keep_all = T)
save(variants_chrV, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS7_variants_chrV.Rda")


#####################################
#       File S8 - NILgeno(VCF)      #
#####################################

# need to make still


#####################################
#       File S9 - NILgeno           #
#####################################

load("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/nilgeno.Rda")

nil_genotypes <- nilgeno %>%
    dplyr::filter(sample %in% c("N2", "CB4856", "ECA1059", "ECA1066", "ECA1065", "ECA377", "ECA232", "ECA629", "ECA634", "ECA636", "ECA632", 'ECA554','ECA573'))

save(nil_genotypes, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS9_nilgenotypes.Rda")

#####################################
#       File S10- CSSV              #
#####################################
# load data

# dose response from briana 2017
cssVpruned <- easysorter::read_data("~/Dropbox/HTA/Results/20170123_cssVa/") %>%
    easysorter::remove_contamination() %>%
    easysorter::sumplate(directories = FALSE, quantiles = TRUE) %>%
    easysorter::bioprune() %>%
    tidyr::gather(trait, phenotype, -c(date:col))%>%
    dplyr::ungroup() %>%
    dplyr::filter(trait %in% c("mean.TOF", "mean.EXT", "norm.n"),
                  condition %in% c("abamectin", "DMSO")) %>%
    easysorter::prune_outliers() %>%
    dplyr::mutate(experiment = "CSSV",
                  round = 1,
                  assay = "a")

save(cssVpruned, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS10_cssV.Rda")

#####################################
#       File S11 - stats            #
#####################################

source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_stats.R")

# load data
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS12_NILpheno.Rda")
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS16_VCNIL.Rda")
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS17_lgc54del.Rda")
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS10_cssV.Rda")

# regress
regressed_nil <- easysorter::regress(NILpruned)
regressed_grace <- easysorter::regress(VC_NIL)
regressed_lgc <- easysorter::regress(lgc54del)
regressed_cssV <- easysorter::regress(cssVpruned)

# calculate stats
stats_nil <- get_stats(regressed_nil) %>%
    dplyr::mutate(exp = "NILpruned")
stats_grace <- get_stats(regressed_grace) %>%
    dplyr::mutate(exp = "VC_NIL")
stats_lgc <- get_stats(regressed_lgc) %>%
    dplyr::mutate(exp = "lgc54del")
stats_dmso <- get_stats(lgc54del %>% dplyr::filter(condition == "DMSO")) %>%
    dplyr::mutate(exp = "lgc54del")
stats_cssV <- get_stats(regressed_cssV) %>%
    dplyr::mutate(exp = "cssV")

# combine
HTAstats <- stats_nil %>%
    dplyr::bind_rows(stats_grace, stats_lgc, stats_dmso, stats_cssV)
    
# save
save(HTAstats, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS11_HTAstats.Rda")

#####################################
#       File S12 - NILpheno         #
#####################################

# use assay b to show
load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/NILs/abaVnil_pruned1b.Rda")

NILpruned <- abaVnil_pruned1b %>%
    dplyr::mutate(experiment = "NILpruned") %>%
    dplyr::filter(trait %in% c("mean.EXT", "mean.TOF", "norm.n"),
                  strain != "ECA1093")

save(NILpruned, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS12_NILpheno.Rda")

#####################################
#      File S13 - candidate genes   #
#####################################

source("~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/NIL_narrowing.R")

VC <- query_genes("V:5260997-5906132", v = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/genome/WI.20200815.hard-filter.vcf.gz")
VC <- VC %>%
    dplyr::mutate(QTL_region = "VC")

VC_genes <- VC %>%
    dplyr::select(-go_term, -go_name, -go_description, -go_annotation) %>%
    dplyr::distinct() %>%
    dplyr::mutate(QTL_region = "VC")

VL <- query_genes("V:2629324-3076312", v = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/genome/WI.20200815.hard-filter.vcf.gz")
VL <- VL %>%
    dplyr::mutate(QTL_region = "VL")
VL_genes <- VL %>%
    dplyr::select(-go_term, -go_name, -go_description, -go_annotation) %>%
    dplyr::distinct() %>%
    dplyr::mutate(QTL_region = "VL")

candidategenes <- rbind(VL, VC)
save(candidategenes, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS13_candidategenes.Rda")

#####################################
#      File S14 - set1 RIAIL        #
#####################################

# load abamectin set 1 phenotypes
linkagemapping::load_cross_obj("N2xCB4856cross_full")
load("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/RIAIL_data/allRIAILsregressed.Rda")

set1_strains <- N2xCB4856cross_full2$pheno %>%
    dplyr::filter(set == 1)

set1pheno <- allRIAILsregressed %>%
    dplyr::filter(condition == "abamectin", 
                  trait %in% c("mean.EXT", "mean.TOF", "mean.norm.EXT", "norm.n"),
                  strain %in% set1_strains$strain)
save(set1pheno, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS14_set1pheno.Rda")

#####################################
#      File S15 - mediation         #
#####################################

# load abamectin set 1 phenotypes
linkagemapping::load_cross_obj("N2xCB4856cross_full")
data("eQTLpeaks")
load("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/RIAIL_data/allRIAILsregressed.Rda")

set1_strains <- N2xCB4856cross_full2$pheno %>%
    dplyr::filter(set == 1)

# VL - mean.EXT
ext <- allRIAILsregressed %>%
    dplyr::filter(condition == "abamectin", 
                  trait == "mean.EXT",
                  strain %in% set1_strains$strain) %>%
    dplyr::mutate(newpheno = (phenotype - mean(phenotype, na.rm = T)) / sd(phenotype, na.rm = T)) %>%
    dplyr::select(strain, trait, phenotype = newpheno)

eqtl_VL <- eQTLpeaks %>%
    dplyr::filter(chr == "V",
                  ci_l_pos < 3076312, 
                  ci_r_pos > 2629324)

med_VL <- NULL
for(i in unique(eqtl_VL$trait)) {
    med_VL <- med_VL %>%
        dplyr::bind_rows(calc_mediation(peak = "V_2895960", expression_probe = i, phenodf = ext, cross = N2xCB4856cross_full2))
}

med_VL <-  med_VL %>%
    dplyr::left_join(eqtl_VL %>% dplyr::distinct(marker, trait), by = c("probe" = "trait")) %>%
    dplyr::rename(eqtl_marker = marker) %>%
    dplyr::mutate(abs_est = abs(estimate))

# VC - mean.EXT (no QTL in set1 though...)
eqtl_VC <- eQTLpeaks %>%
    dplyr::filter(chr == "V",
                  ci_l_pos < 5906132, 
                  ci_r_pos > 5260997)

# # need to find a marker...
# markers <- data.frame(N2xCB4856cross_full2$geno$V$data) %>%
#     tidyr::gather(marker, genotype) %>%
#     dplyr::distinct(marker) %>%
#     tidyr::separate(marker, into = c("chr", "pos")) %>%
#     dplyr::mutate(pos = as.numeric(pos)) %>%
#     dplyr::filter(pos < 5906132, pos > 5260997)

med_VC <- NULL
for(i in unique(eqtl_VC$trait)) {
    med_VC <- med_VC %>%
        dplyr::bind_rows(calc_mediation(peak = "V_5585423", expression_probe = i, phenodf = ext, cross = N2xCB4856cross_full2))
}

med_VC <-  med_VC %>%
    dplyr::left_join(eqtl_VC %>% dplyr::distinct(marker, trait), by = c("probe" = "trait")) %>%
    dplyr::rename(eqtl_marker = marker) %>%
    dplyr::mutate(abs_est = abs(estimate))

aba_mediation <- med_VL %>%
    dplyr::bind_rows(med_VC)
save(aba_mediation, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS15_mediation.Rda")



#####################################
#      File S16 - grace NIL         #
#####################################

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/NILs/20190313_abaV/abaVnil_pruned1a.Rda")

VC_NIL <- abaVnil_pruned1a %>%
    dplyr::filter(strain %in% c("N2", "CB4856", "ECA629", "ECA634", "ECA636", "ECA632"),
                  trait %in% c("mean.TOF", "mean.EXT", "norm.n")) %>%
    dplyr::mutate(round = 1, assay = "a", experiment = "VC_NIL")

save(VC_NIL, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS16_VCNIL.Rda")


#####################################
#      File S17 - lgc-54 del        #
#####################################

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/data/NILs/abamectinNILs_pruned1.Rda")

lgc54del <- abamectinNILs_pruned1 %>%
    dplyr::filter(strain %in% c("N2", "CB4856", "FX3518", "FX3448"),
                  trait %in% c("norm.n", "mean.TOF", "mean.EXT")) %>%
    dplyr::mutate(round = 1, assay = "a",
                  condition = ifelse(condition == "Abamectin", "abamectin", "DMSO")) %>%
    dplyr::select(date:phenotype)
save(lgc54del, file = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS17_lgc54del.Rda")

#####################################
#    File S18 - parasite overlap    #
#####################################

hc_overlap <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1uPsF5NrXSSpjhLaOFSARog-hKPi8exhqtG3oZpzSO5M/edit#gid=146627604")
hc_overlap[hc_overlap == "-"] <- NA
hc_overlap$HCON_start <- as.numeric(as.character(hc_overlap$HCON_start))
hc_overlap$HCON_stop <- as.numeric(as.character(hc_overlap$HCON_stop))
readr::write_tsv(hc_overlap, "~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS18_hc_overlap.tsv")





#####################################
#       Supplement Files - csv      #
#####################################

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS1_newdose.Rda")
write.csv(newdose, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS1_newdose.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS2_wipheno.Rda")
write.csv(wipheno, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS2_wipheno.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS3_wimap.Rda")
write.csv(wimap, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS3_wimap.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS4_riailpheno.Rda")
write.csv(riailpheno,"~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS4_riailpheno.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS5_riailmap.Rda")
write.csv(riailmap, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS5_riailmap.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS6_scan2summary.Rda")
write.csv(scan2summary, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS6_scan2summary.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS7_variants_chrV.Rda")
write.csv(variants_chrV, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS7_variants_chrV.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS9_nilgenotypes.Rda")
write.csv(nil_genotypes, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS9_nilgenotypes.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS10_cssV.Rda")
write.csv(cssVpruned, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS10_cssV.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS11_HTAstats.Rda")
write.csv(HTAstats, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS11_HTAstats.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS12_NILpheno.Rda")
write.csv(NILpruned, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS12_NILpheno.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS13_candidategenes.Rda")
write.csv(candidategenes, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS13_candidategenes.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS14_set1pheno.Rda")
write.csv(set1pheno, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS14_set1pheno.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS15_mediation.Rda")
write.csv(aba_mediation, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS15_mediation.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS16_VCNIL.Rda")
write.csv(VC_NIL, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS16_VCNIL.csv")

load("~/Dropbox/AndersenLab/LabFolders/Katie/projects/abamectin/manuscript/data/FileS17_lgc54del.Rda")
write.csv(lgc54del, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS17_lgc54del.csv")

hc_overlap <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1uPsF5NrXSSpjhLaOFSARog-hKPi8exhqtG3oZpzSO5M/edit#gid=146627604")
hc_overlap[hc_overlap == "-"] <- NA
write.csv(hc_overlap, "~/Dropbox/AndersenLab/LabFolders/Katie/git/abamectin_manuscript/data/FileS18_hc_overlap.csv")
