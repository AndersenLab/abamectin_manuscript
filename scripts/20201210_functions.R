# functions for plotting
library(tidyverse)
library(linkagemapping)
library(cegwas2)


data("eQTLpeaks")
data("eqtlpheno")
load("data/extra/gene_annotations.Rda")

# linkage mapping function 
plot_lods <- function(map, tsize = 12) {
    map1 <- map %>% 
        dplyr::group_by(marker, condtrt) %>% 
        dplyr::filter(lod == max(lod))
    cis <- map %>% 
        dplyr::group_by(marker, condtrt) %>% 
        dplyr::mutate(maxlod = max(lod)) %>%
        dplyr::group_by(iteration, condtrt) %>% 
        dplyr::filter(!is.na(var_exp)) %>%
        dplyr::do(head(., n = 1))
    
    totalmap <- NULL
    if(nrow(cis) > 0) {
        for(i in unique(cis$condtrt)) {
            drugci <- cis %>%
                dplyr::filter(condtrt == i)
            drugmap <- map1 %>%
                dplyr::filter(condtrt == i)
            map2 <- linkagemapping:::cidefiner(drugci, drugmap)
            totalmap <- rbind(totalmap, map2)
        }
        
        totalmap$condtrt <- gsub("_", "\n", totalmap$condtrt)
        cis$condtrt <- gsub("_", "\n", cis$condtrt)
        
        ggplot2::ggplot(totalmap) + 
            ggplot2::aes(x = pos/1e+06,y = lod) + 
            ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e+06,ymin = 0, ymax = ci_lod), fill = "skyblue", alpha = 0.5) +
            ggplot2::geom_point(data = cis, ggplot2::aes(x = pos/1e+06,y = (1.05 * maxlod)), fill = "red", shape = 25,
                                size = tsize/7, show.legend = FALSE, color = "red") + 
            # ggplot2::geom_text(data = cis, ggplot2::aes(x = pos/1e+06, y = (1.2 * maxlod), label = paste0(100 *round(var_exp, digits = 3), "%")),
                               # size = tsize/5, colour = "black") +
            ggrepel::geom_text_repel(data = cis, ggplot2::aes(x = pos/1e+06, y = (1.2* maxlod), label = paste0(100 *round(var_exp, digits = 3), "%")),
                               size = tsize/5, colour = "black") +
            ggplot2::geom_line(size = tsize/25, alpha = 0.85) +
            ggplot2::facet_grid(~ chr, scales = "free", space = "free_x") +
            ggplot2::labs(x = "Genomic position (Mb)", y = "LOD") +
            ggplot2::scale_colour_discrete(name = "Mapping\nIteration") +
            # ggplot2::scale_x_continuous(expand = c(0,0)) + 
            # ggplot2::scale_y_continuous(expand = expand_scale(mult=c(0,0.1))) +
            ggplot2::theme_bw(tsize) +
            ggplot2::theme(
                axis.text = ggplot2::element_text(color = "black", face = "bold"),
                axis.title = ggplot2::element_text(face = "bold", color = "black"),
                strip.text = ggplot2::element_text(face = "bold", color = "black"),
                plot.title = ggplot2::element_blank(),
                panel.grid = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(color = NA, size = 0.6))
    } else {
        totalmap <- map1
        totalmap$condtrt <- gsub("_", "\n", totalmap$condtrt)
        
        ggplot2::ggplot(totalmap) + 
            ggplot2::aes(x = pos/1e+06,y = lod) + 
            ggplot2::geom_line(size = tsize/25, alpha = 0.85) +
            ggplot2::facet_grid( ~ chr, scales = "free", space = "free_x") +
            ggplot2::labs(x = "Genomic position (Mb)", y = "LOD") +
            ggplot2::theme_bw(tsize) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(color = "black", face = "bold"),
                axis.text.y = ggplot2::element_text(face = "bold", color = "black"),
                axis.title.x = ggplot2::element_text(face = "bold", color = "black"),
                axis.title.y = ggplot2::element_text(face = "bold", color = "black"),
                strip.text.x = ggplot2::element_text(face = "bold", color = "black"),
                strip.text.y = ggplot2::element_text(face = "bold", color = "black"),
                plot.title = ggplot2::element_blank(),
                panel.grid = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(color = NA, size = 0.6))
    }
}

# function for PxG 
plot_pxg <- function(cross, map, tsize = 12, yaxis = "Relative animal length") {
    
    # get unique QTL peaks
    peaks <- map %>% 
        dplyr::group_by(iteration, condition) %>% 
        dplyr::filter(!is.na(var_exp)) %>% 
        dplyr::do(head(., n = 1))
    
    # clean the markers and column names
    uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
    colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
    colnames(cross$pheno) <- stringr::str_replace(colnames(cross$pheno), "\\.", "_")
    
    # get only the traits of interest
    pheno <- cross$pheno %>% 
        dplyr::select(dplyr::one_of(map$condtrt))
    
    # get the genotype for the RIAILs and add to pheno
    geno <- data.frame(linkagemapping:::extract_genotype(cross)) %>% 
        dplyr::select(which(colnames(.) %in% uniquemarkers)) %>% 
        data.frame(., pheno)
    
    # reorder data and change -1 and 1 to N2 and CB and plot!
    df <- geno %>%
        tidyr::gather(marker, genotype, -dplyr::one_of(map$condtrt)) %>%
        dplyr::mutate(genotype = dplyr::case_when(genotype == -1 ~ "N2",
                                                  genotype == 1 ~ "CB4856",
                                                  TRUE ~ "NA")) %>%
        tidyr::gather(trait, phenotype, dplyr::one_of(map$condtrt)) %>%
        dplyr::mutate(genotype = factor(genotype, levels = c("N2", "CB4856"), labels= c("N2", "CB4856"))) %>%
        dplyr::left_join(peaks, by = c("marker", "trait" = "condtrt")) %>%
        tidyr::drop_na(lod) %>%
        dplyr::mutate(marker = stringr::str_replace(marker, "_", ":"),
                      trait = stringr::str_split_fixed(trait, "_", 2)[,2],
                      pos = as.numeric(stringr::str_split_fixed(marker, ":", 2)[,2])) %>%
        dplyr::arrange(chr, pos) %>%
        tidyr::drop_na(genotype)
    df %>%
        ggplot2::ggplot(.) + 
        ggplot2::aes(x = genotype, y = phenotype) +
        ggplot2::geom_jitter(width = 0.1, size = 0.07, alpha = 0.5) + 
        ggplot2::geom_boxplot(ggplot2::aes(fill = genotype, alpha = 0.5), size = 0.2, outlier.shape = NA) + 
        ggplot2::scale_fill_manual(values = c(`N2` = "orange", `CB4856` = "blue")) + 
        ggplot2::facet_grid(~factor(marker, unique(df$marker)), scales = "free") + 
        ggplot2::theme_bw(tsize) + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", color = "black"), 
                       axis.text.y = ggplot2::element_text(face = "bold", color = "black"), 
                       axis.title.x = ggplot2::element_text(face = "bold", color = "black", vjust = -0.3),
                       axis.title.y = ggplot2::element_text(face = "bold", color = "black"), 
                       strip.text.x = ggplot2::element_text(face = "bold", color = "black"),
                       strip.text.y = ggplot2::element_text(face = "bold", color = "black"),
                       # axis.title.x = element_blank(),
                       plot.title = ggplot2::element_blank(), 
                       legend.position = "none", 
                       panel.grid = element_blank(),
                       panel.background = ggplot2::element_rect(color = NA, size = 0.6)) +
        ggplot2::labs(x = "Genotype at QTL", y = yaxis)
}

# plot NIL pheno and geno with stats
plot_nil <- function(phenodf, genodf, statdf, strains, chr, tsize = 12, ylab = "Animal length") {

    # plot phenotype
    pheno <- phenodf %>%
        dplyr::group_by(strain, condition) %>%
        dplyr::mutate(phen = max(phenotype) + 0.1) %>%
        dplyr::ungroup() %>%
        dplyr::full_join(statdf, by = c("strain", "condition", "trait")) %>%
        dplyr::mutate(strain = factor(strain, 
                                      levels = rev(strains))) %>%
        dplyr::filter(!is.na(strain),
                      strain %in% strains) %>%
        ggplot(.) +
        aes(x = strain, y = phenotype, fill = strain_fill) +
        geom_jitter(width = 0.1, size = 0.05) +
        geom_boxplot(outlier.color = NA, alpha = 0.5, size = 0.2) +
        ggplot2::geom_text(aes(label = sig, y = phen, color = groups), size = tsize/4, angle = -90) +
        scale_fill_manual(values = c("N2" = "orange", "CB" = "blue", "NIL" = "grey")) +
        scale_color_manual(values = c("N2" = "orange", "CB" = "blue")) +
        theme_bw(tsize) +
        theme(axis.text.x = element_text(face="bold", color="black"),
              axis.title.x = element_text(face="bold", color="black"),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              strip.text = element_text(face = "bold", color = "black"),
              legend.position = "none",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank()) +
        coord_flip() +
        labs(x = " ", y = ylab)  +
        facet_grid(~condition)
    
    # plot chrV genotype
    chrgeno <- genodf %>%
        dplyr::filter(chrom == chr,
                      sample %in% strainset) %>%
        dplyr::mutate(chrom = paste0("chr", chr),
                      sample = factor(sample, levels = rev(strainset))) %>%
        dplyr::distinct() %>%
        ggplot(.)+
        geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2), alpha = 0.7)+
        facet_grid(~chrom, scales = "free",  space = "free")+
        scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
        theme_bw(tsize) +
        theme(axis.text.x = element_text(face="bold", color="black"),
              axis.text.y = element_text(face="bold", color="black"),
              axis.ticks.y = element_blank(),
              # axis.text.y = element_blank(),
              axis.title.x = element_text(face="bold", color="black"),
              axis.title.y = element_text(face="bold", color="black"),
              strip.text = element_text(face = "bold", color = "black"),
              plot.title = element_text(face="bold"),
              legend.position = "none",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())+
        labs(x = "Genomic position (Mb)", y = "") 
    
    # plot genome background
    back <- genodf %>%
        dplyr::filter(sample %in% strainset,
                      chrom == ifelse(chr == "III", "II", "III")) %>%
        dplyr::mutate(bp = end - start) %>%
        dplyr::group_by(sample, gt_name) %>%
        dplyr::summarize(max = sum(bp)) %>%
        dplyr::arrange(desc(max)) %>%
        dplyr::filter(max == first(max)) %>%
        dplyr::mutate(start = 0,
                      end = 15e6,
                      chr = "Genome")  %>%
        dplyr::ungroup() %>%
        dplyr::mutate(sample = factor(sample, levels = rev(strainset))) %>%
        ggplot(.)+
        geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2), alpha = 0.7)+
        facet_grid(~chr, scales = "free",  space = "free")+
        scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
        theme_bw(tsize) +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              legend.position = "none",
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              strip.text.x = element_text(face = "bold", color = "black"),
              strip.text.y = element_blank())+
        labs(x = "Genomic position (Mb)", y = "")
    
    # return objects
    return(list(chrgeno, back, pheno))
}

# Look for genes in interval
query_genes <- function(region, GO = NULL, strain = "CB4856", v = "~/Dropbox/AndersenLab/LabFolders/Katie/projects/genome/WI.20200815.hard-filter.vcf.gz") {
    
    # filter eqtl to > 5% VE
    eqtlmap2 <- eQTLpeaks %>%
        dplyr::filter(var_exp >= 0.05)
    
    # how many genes are in the interval?
    all_genes <- cegwas2::query_vcf(region, impact = "ALL", samples = strain, vcf = v)
    print(glue::glue("There are {length(unique(all_genes$gene_id))} genes in the interval {region}"))
    
    # how many eQTL map to this region?
    chrom <- stringr::str_split_fixed(region, ":", 2)[,1]
    left_pos <- as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(region, ":", 2)[,2], "-", 2)[,1])
    right_pos <- as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(region, ":", 2)[,2], "-", 2)[,2])
    
    all_eQTL <- eqtlmap2 %>%
        dplyr::filter(chr == chrom,
                      ci_l_pos < right_pos,
                      ci_r_pos > left_pos)
    print(glue::glue("There are {nrow(all_eQTL)} eQTL ({length(unique(all_eQTL$trait))} traits) with VE > 5% that map to {region}"))
    
    # all eQTL probes
    all_eQTL_probes <- probe_info %>%
        dplyr::filter(probe %in% all_eQTL$trait)
    
    ##############################
    # if wbgene is NA - try to fix
    ##############################
    
    # filter na
    na1 <- all_eQTL_probes %>%
        dplyr::group_by(probe) %>%
        dplyr::mutate(num_na = sum(is.na(wbgene))/length(wbgene)) %>%
        dplyr::filter(num_na == 1)
    
    unique_probes <- paste(unique(na1$probe), collapse = ",")
    print(glue::glue("There are {nrow(na1)} genes with an eQTL that need to be hand curated: {unique_probes}"))
    
    ##################################
    
    # which of the eQTL are overlapping with genes in interval?
    eQTL_outside_CI <- all_eQTL_probes %>%
        dplyr::filter(!wbgene %in% all_genes$gene_id)
    print(glue::glue("There are {nrow(all_eQTL)-length(unique(eQTL_outside_CI$wbgene))-nrow(na1)} genes in the region with an eQTL, {length(unique(eQTL_outside_CI$wbgene))} genes outside the region with an eQTL, and {nrow(na1)} unknown"))
    
    # Total genes of interest:
    print(glue::glue("There are at least {length(unique(all_genes$gene_id)) + length(unique(eQTL_outside_CI$wbgene))} total genes of interest."))
    
    # how many of the genes in interval have variation?
    vars <- all_genes %>%
        dplyr::mutate(GT = ifelse(a1 == REF, "ref", "alt")) %>%
        dplyr::filter(GT == "alt")
    
    # genes with protein coding vars
    proteincode <- vars %>%
        dplyr::filter(impact %in% c("MODERATE", "HIGH"))
    print(glue::glue("There are {length(unique(vars$gene_id))}/{length(unique(all_genes$gene_id))} genes in interval with genetic variation, {length(unique(proteincode$gene_id))}/{length(unique(vars$gene_id))} have protein-coding variation"))
    
    
    # should I look at GO annotations?
    if(!is.null(GO)) {
        # total genes with GO annotations
        go_genes <- gene_annotations %>%
            dplyr::filter(wbgene %in% c(all_genes$gene_id, eQTL_outside_CI$wbgene)) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_genes$wbgene))}/{length(unique(all_genes$gene_id)) + length(unique(eQTL_outside_CI$wbgene))} genes with {GO} annotation"))
        
        # genes with GO annotations and variation
        go_var <- gene_annotations %>%
            dplyr::filter(wbgene %in% vars$gene_id) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_var$wbgene))}/{length(unique(go_genes$wbgene))} genes with {GO} annotation AND genetic variation"))
        
        # genes with GO annotation and protein-coding variation
        go_pcvar <- gene_annotations %>%
            dplyr::filter(wbgene %in% proteincode$gene_id) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_pcvar$wbgene))}/{length(unique(go_genes$wbgene))} genes with {GO} annotation AND protein-coding genetic variation"))
        
        # genes with GO annotation and eQTL
        go_eqtl <- gene_annotations %>%
            dplyr::filter(wbgene %in% all_eQTL_probes$wbgene) %>%
            dplyr::filter_all(any_vars(stringr::str_detect(., pattern = GO)))
        print(glue::glue("There are {length(unique(go_eqtl$wbgene))}/{length(unique(go_genes$wbgene))} genes with {GO} annotation AND eQTL"))
        
        # return final dataframe with all info (might be off, only has 133 instead of 134?)
        total_genes <- gene_annotations %>%
            dplyr::filter(wbgene %in% c(all_genes$gene_id, eQTL_outside_CI$wbgene)) %>%
            dplyr::mutate(inside_CI = ifelse(wbgene %in% all_genes$gene_id, T, F),
                          eqtl = ifelse(wbgene %in% all_eQTL_probes$wbgene, T, F),
                          vars = ifelse(wbgene %in% vars$gene_id, T, F),
                          pc_vars = ifelse(wbgene %in% proteincode$gene_id, T, F),
                          go_annotation = ifelse(wbgene %in% go_genes$wbgene, T, F))
        
        distinct <- total_genes %>%
            dplyr::distinct(wbgene, .keep_all = T)
        
        # pc alone
        pc <- nrow(distinct %>% dplyr::filter(pc_vars == T, eqtl == F))
        
        # pc + eQTL
        pc_eqtl <- nrow(distinct %>% dplyr::filter(pc_vars == T, eqtl == T))
        
        # eQTL alone
        e <- nrow(distinct %>% dplyr::filter(pc_vars == F, eqtl == T))
        
        print(glue::glue("There are {pc + pc_eqtl + e} genes with protein-coding variation and/or an eQTL (top priority)"))
    } else {
        
        # return final dataframe with all info (might be off, only has 133 instead of 134?)
        total_genes <- gene_annotations %>%
            dplyr::filter(wbgene %in% c(all_genes$gene_id, eQTL_outside_CI$wbgene)) %>%
            dplyr::mutate(inside_CI = ifelse(wbgene %in% all_genes$gene_id, T, F),
                          eqtl = ifelse(wbgene %in% all_eQTL_probes$wbgene, T, F),
                          vars = ifelse(wbgene %in% vars$gene_id, T, F),
                          pc_vars = ifelse(wbgene %in% proteincode$gene_id, T, F),
                          go_annotation = NA)
        
        distinct <- total_genes %>%
            dplyr::distinct(wbgene, .keep_all = T)
        
        # pc alone
        pc <- nrow(distinct %>% dplyr::filter(pc_vars == T, eqtl == F))
        
        # pc + eQTL
        pc_eqtl <- nrow(distinct %>% dplyr::filter(pc_vars == T, eqtl == T))
        
        # eQTL alone
        e <- nrow(distinct %>% dplyr::filter(pc_vars == F, eqtl == T))
        
        print(glue::glue("There are {pc + pc_eqtl + e} genes with protein-coding variation and/or an eQTL (top priority)"))
    }
    
    return(total_genes)
}    

# NIL genotype plots
nilgeno <- nil_genotypes
nil_plot <- function(strains, chr, left.cb = 0, left.n2 = 0, left.bound = 0, right.bound = 19e6, scan.range = 2e4, 
                     all.chr=F, section = "all", background = F, ci = NA, elements = F, order = T, n2_color = "orange", cb_color = "blue",
                     tsize = 12){
    
    # copy nilgeno
    nilgeno2 <- nilgeno
    
    # check to see if all strains are represented in the nilgeno2 file
    df <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::distinct(sample)
    
    if(nrow(df) != length(strains)) {
        unknown_strains <- setdiff(strains, df$sample)
        
        # make fake genotype for unknown strains
        for(s in unknown_strains) {
            fake_geno <- nilgeno2 %>%
                dplyr::filter(sample == "N2")%>%
                dplyr::distinct(chrom, .keep_all = T) %>%
                dplyr::mutate(sample = s,
                              gt = 3,
                              gt_name = "unknown")
            
            # add fake geno to nil geno
            nilgeno2 <- rbind(nilgeno2, fake_geno)
        }
    }
    
    # # # determine if NILs are CB or N2
    nilsII_sort_type <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::group_by(sample)%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::group_by(sample, start)%>%
        dplyr::mutate(gt_ct = sum(size))%>%
        dplyr::ungroup()%>%
        dplyr::group_by(sample, gt)%>%
        dplyr::mutate(major_gt = sum(gt_ct))%>%
        dplyr::ungroup()%>%
        dplyr::arrange(desc(major_gt))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::mutate(nil_type = ifelse(gt == 1, "CB", "N2"))%>%
        dplyr::select(sample, nil_type)
    
    # # # keep NILs that lost NIL genotype on right side
    nilsII_left <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::filter((start > left.cb - .3e6 & start < left.cb + .3e6 & nil_type == "CB") |
                          (start > left.n2 - .3e6 & start < left.n2 + .3e6 & nil_type == "N2") )%>%
        dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
        dplyr::mutate(side = "LEFT")%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::ungroup()%>%
        dplyr::arrange(desc(size))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::filter(size > scan.range)%>% # # # remove small (likely wrong calls) around interval site
        dplyr::select(sample, side)
    
    # # # keep NILs that lost NIL genotype on left side
    nilsII_right <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::filter(!sample %in% nilsII_left$sample)%>%
        dplyr::mutate(side = "RIGHT")%>%
        dplyr::distinct(sample,.keep_all = T)%>%
        dplyr::select(sample, side)
    
    nil_sides <- bind_rows(nilsII_left,nilsII_right)
    
    
    nilsII_sort_left <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::left_join(.,nil_sides, by = "sample")%>%
        dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
        dplyr::group_by(sample)%>%
        dplyr::filter(start > left.bound & start < right.bound)%>%
        dplyr::filter(end > left.bound | end < right.bound)%>%
        dplyr::filter(end > left.bound | end < right.bound)%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::group_by(sample, start)%>%
        dplyr::mutate(gt_ct = sum(size))%>%
        dplyr::ungroup()%>%
        dplyr::filter(side == "LEFT")%>%
        dplyr::arrange( desc(gt_ct))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::arrange(nil_type, desc(gt_ct))
    
    #Here
    nilsII_sort_right <- nilgeno2 %>%
        dplyr::filter(sample %in% strains) %>%
        dplyr::filter(chrom == chr)%>%
        dplyr::left_join(.,nilsII_sort_type, by = "sample")%>%
        dplyr::left_join(.,nil_sides, by = "sample")%>%
        dplyr::filter(side == "RIGHT")%>%
        dplyr::filter(gt == 2 & nil_type == "CB" | gt == 1 & nil_type == "N2")%>%
        dplyr::group_by(sample)%>%
        dplyr::filter(start > left.bound & start < right.bound)%>%
        dplyr::filter(end > left.bound | end < right.bound)%>%
        dplyr::mutate(size = end - start )%>%
        dplyr::group_by(sample, start)%>%
        dplyr::mutate(gt_ct = sum(size))%>%
        dplyr::ungroup()%>%
        dplyr::arrange( desc(gt_ct))%>%
        dplyr::distinct(sample, .keep_all = T)%>%
        dplyr::ungroup()%>%
        dplyr::arrange( nil_type, gt_ct)
    
    N2CB <- nilgeno2 %>%
        dplyr::filter(sample %in% c("N2", "CB4856"), chrom == chr) %>%
        dplyr::mutate(nil_type = "parent", side = NA, size = NA, gt_ct = NA)
    
    nilsII_sort <- bind_rows(nilsII_sort_right, nilsII_sort_left) %>%
        arrange(desc(nil_type), desc(side), size)
    
    nilsII_sort <- rbind(nilsII_sort, N2CB)
    
    if(nrow(df) != length(strains)) {
        nilsII_sort <- rbind(nilsII_sort, nilgeno2 %>%
                                 dplyr::filter(sample %in% unknown_strains,
                                               chrom == chr) %>%
                                 dplyr::mutate(nil_type = "unkown", side = NA, size = NA, gt_ct = NA))
        
    }
    
    if (all.chr == T){
        
        nilsII <- nilgeno2 %>%
            dplyr::filter(sample %in% strains) %>%
            dplyr::filter(chrom != "MtDNA")%>%
            dplyr::left_join(.,nilsII_sort_type, by = "sample")
        
        if(nrow(df) != length(strains)) {
            nilsII <- nilsII %>%
                dplyr::mutate(nil_type = ifelse(gt_name == "unknown", "unknown", nil_type))
        }
        
        if(order == T) {
            nilsII$sample <- factor(nilsII$sample, levels = unique(nilsII_sort$sample), labels = unique(nilsII_sort$sample), ordered = T)
        } else {
            nilsII$sample <- factor(nilsII$sample, levels = unique(strains), ordered = T)
        }
        nilsII$gt <- as.character(nilsII$gt)
    } else {
        nilsII <- nilgeno2 %>%
            dplyr::filter(sample %in% strains) %>%
            dplyr::filter(chrom == chr, chrom != "MtDNA")%>%
            dplyr::left_join(.,nilsII_sort_type, by = "sample")
        
        if(nrow(df) != length(strains)) {
            nilsII <- nilsII %>%
                dplyr::mutate(nil_type = ifelse(gt_name == "unknown", "unknown", nil_type))
        }
        
        if(order == T) {
            nilsII$sample <- factor(nilsII$sample, levels = unique(nilsII_sort$sample), labels = unique(nilsII_sort$sample), ordered = T)
        } else {
            nilsII$sample <- factor(nilsII$sample, levels = unique(strains), ordered = T)
        }
        nilsII$gt <- as.character(nilsII$gt)
    }
    
    # make fake 'background' genome
    if(background == T) {
        # cannot show all chromosomes and the background
        if(all.chr == T) {
            bgplot <- NULL
        } else {
            bg <- data.frame(chrom = NA, start = NA, end = NA, sample = NA, gt = NA, gt_name = NA, supporting_sites = NA, sites = NA,
                             DP = NA, switches = NA, nil_type = NA)
            for(i in unique(nilsII$sample)) {
                # get the NIL type
                type <- nilsII %>%
                    dplyr::filter(sample == i) %>%
                    dplyr::distinct(nil_type) %>%
                    dplyr::pull(nil_type)
                
                # add a new row with the background genotype in a "new" chromosome
                bg <- rbind(bg, c("Genome", 1, 3e6, i, 
                                  ifelse(type == "N2", 2,
                                         ifelse(type == "CB", 1,
                                                ifelse(type == "unknown", 3, NA))),
                                  ifelse(type == "N2", "CB4856",
                                         ifelse(type == "CB", "N2",
                                                ifelse(type == "unknown", "unknown", NA))),
                                  NA, NA, NA, NA, type))
            }
            bg <- bg %>%
                tidyr::drop_na(sample) %>%
                dplyr::mutate(start = as.numeric(start), end = as.numeric(end))
            
            # add "chrom" to chromosome for plotting
            nils2 <- nilsII %>%
                dplyr::mutate(chrom = paste0("chr", chrom))
            
            # background plot
            bgplot <- ggplot(bg)+
                geom_segment(aes(x = start/1e6, y = factor(sample, levels = levels(nils2$sample)), xend = end/1e6, yend = sample, color = gt_name, size = 2))+
                scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
                # scale_color_manual(values=c("1"="#F0BA51","2"="#484DA0"))+
                facet_grid(~chrom, scales = "free",  space = "free")+
                theme_bw(tsize) +
                theme(axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      strip.text = element_text(face = "bold", color = "black"),
                      plot.title = element_text(face="bold"),
                      legend.position = "none",
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank(),
                      axis.ticks = element_blank())
        }
    } else { bgplot <- NULL }
    
    # only plot the N2-NILs
    if(section == "N2-NILs") {
        nl.pl <- ggplot(nilsII %>% dplyr::filter(nil_type == "CB"))+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
            theme_bw(tsize) +
            theme(axis.text.x = element_text(face="bold", color="black"),
                  axis.text.y = element_text(face="bold", color="black"),
                  axis.title.x = element_text(face="bold", color="black"),
                  axis.title.y = element_text(face="bold", color="black"),
                  strip.text = element_text(face = "bold", color = "black"),
                  plot.title = element_text(face="bold"),
                  legend.position = "none",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())+
            geom_vline(data = filter(nilsII, chrom == chr) %>%
                           dplyr::distinct(sample, chrom) %>%
                           dplyr::mutate(xint = paste(ci, collapse = ",")) %>%
                           tidyr::separate_rows(xint), 
                       aes(xintercept = as.numeric(xint)), color = "red") +
            labs(x = "Genomic position (Mb)", y = "")
    } else if(section == "CB-NILs") {
        # only plot the CB-NILs
        nl.pl <- ggplot(nilsII %>% dplyr::filter(nil_type == "CB"))+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
            theme_bw(tsize) +
            theme(axis.text.x = element_text(face="bold", color="black"),
                  axis.text.y = element_text(face="bold", color="black"),
                  axis.title.x = element_text(face="bold", color="black"),
                  axis.title.y = element_text(face="bold", color="black"),
                  strip.text = element_text(face = "bold", color = "black"),
                  plot.title = element_text(face="bold"),
                  legend.position = "none",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())+
            geom_vline(data = filter(nilsII, chrom == chr) %>%
                           dplyr::distinct(sample, chrom) %>%
                           dplyr::mutate(xint = paste(ci, collapse = ",")) %>%
                           tidyr::separate_rows(xint), 
                       aes(xintercept = as.numeric(xint)), color = "red") +
            labs(x = "Genomic position (Mb)", y = "")
    } else {
        
        # plot all NILs
        nl.pl <- ggplot(nilsII)+
            geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2))+
            facet_grid(~chrom, scales = "free",  space = "free")+
            scale_color_manual(values=c("N2"=n2_color,"CB4856"=cb_color, "unknown" = "grey"))+
            theme_bw(tsize) +
            theme(axis.text.x = element_text(face="bold", color="black"),
                  axis.text.y = element_text(face="bold", color="black"),
                  axis.title.x = element_text(face="bold", color="black"),
                  axis.title.y = element_text(face="bold", color="black"),
                  strip.text = element_text(face = "bold", color = "black"),
                  plot.title = element_text(face="bold"),
                  legend.position = "none",
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())+
            geom_vline(data = filter(nilsII, chrom == chr) %>%
                           dplyr::distinct(sample, chrom) %>%
                           dplyr::mutate(xint = paste(ci, collapse = ",")) %>%
                           tidyr::separate_rows(xint), 
                       aes(xintercept = as.numeric(xint)), color = "red") +
            labs(x = "Genomic position (Mb)", y = "")
    }
    
    # return plots and dataframe
    if(is.null(bgplot)) {
        return(list(nl.pl, nilsII_sort, nilsII))
    } else {
        if(elements == T) {
            return(list(nl.pl, bgplot, nilsII_sort))
        } else {
            return(list(cowplot::plot_grid(nl.pl, bgplot, nrow = 1, ncol = 2, align = "h", axis = "b", rel_widths = c(1, 0.3)),
                        nilsII_sort,
                        nilsII))
            # return(list(nl.pl,
            #             nilsII_sort,
            #             nilsII,
            #             bgplot))
        }
    }
}
