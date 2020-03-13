library(nichenetr)
library(Seurat)
library(tidyverse)
library(vroom)
library(clusterProfiler)

ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
saveRDS(ligand_target_matrix, "../data/ligand_target_matrix.rds")
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
saveRDS(weighted_networks, "../data/weighted_networks.rds")

# ligand-receptor interaction
## weighted_networks    <- readRDS("weighted_networks.rds")
## # ligand-target interaction
## ligand_target_matrix <- readRDS("ligand_target_matrix.rds")

weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network %>% 
                                                                distinct(from, to), by = c("from", "to"))

primary1_obj <- readRDS('../results/seurat_sara/20696_seurat-object.rds')
Idents(primary1_obj) <- primary1_obj$RNA_snn_res.0.8

library(dplyr)
library(tidyverse)

for (receiver in c('9', '2', '1')) {
    for (sender in c('3', '11', '15')) {
        expressed_genes_receiver = get_expressed_genes(receiver, primary1_obj, pct = 0.10)
        background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
        list_expressed_genes_sender = sender %>% unique() %>% lapply(get_expressed_genes, primary1_obj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
        expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
        seurat_obj_receiver = subset(primary1_obj, idents=receiver)
        allmarkers <- read.table(Sys.glob('../results/seurat_sara/20696*SCT*'))
        geneset_oi = allmarkers[allmarkers$cluster==receiver, 'genename']
        geneset_oi = as.vector(geneset_oi %>% .[. %in% rownames(ligand_target_matrix)])
        ligands = lr_network %>% pull(from) %>% unique()
        receptors = lr_network %>% pull(to) %>% unique()
        expressed_ligands = intersect(ligands, expressed_genes_sender)
        expressed_receptors = intersect(receptors, expressed_genes_receiver)
        potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
        ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
        ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
        best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
        DotPlot(primary1_obj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
        ggsave(paste0('sender', paste(sender, collapse='_'), 'receiver', receiver, '.pdf'))

        active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
        active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
        order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
        order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
        rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
        colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

        vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

        p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0, 0.001, 0.002, 0.005, 0.01))
        p_ligand_target_network
        ggsave('20696_ligand_target_prediction.pdf', width=19, height=7.5)

        lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
        best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
        lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

        lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
        lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

        dist_receptors = dist(lr_network_top_matrix, method = "binary")
        hclust_receptors = hclust(dist_receptors, method = "ward.D2")
        order_receptors = hclust_receptors$labels[hclust_receptors$order]
        
        dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
        hclust_ligands = hclust(dist_ligands, method = "ward.D2")
        order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

        order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
        order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

        vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
        rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
        colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

        p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
        p_ligand_receptor_network
        ggsave('20696_ligand_receptor_prediction.pdf', width=19, height=7.5)

        lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
        ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
        receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

        lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
        lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

        lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
        lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

        dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
        hclust_receptors = hclust(dist_receptors, method = "ward.D2")
        order_receptors = hclust_receptors$labels[hclust_receptors$order]

        dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
        hclust_ligands = hclust(dist_ligands, method = "ward.D2")
        order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

        order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
        order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

        vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
        rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
        colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

        p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
        p_ligand_receptor_network_strict
        ggsave('20696_ligand_receptor_prediction_restrict.pdf', width=19, height=7.5)


                                        # combined heatmap: overlay ligand activities with target genes
        ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

        rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
        colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

        vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
        p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange", legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))

        figures_without_legend = cowplot::plot_grid(p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                                                    p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                                    align = "hv",
                                                    nrow = 1,
                                                    rel_widths = c(ncol(vis_ligand_pearson) + 10, ncol(vis_ligand_target)))
        legends = cowplot::plot_grid(
                               ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
                               ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
                               nrow = 1,
                               align = "h")
        combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,2), nrow = 2, align = "hv")
        combined_plot
        ggsave(paste0('sender', paste(sender, collapse='_'), 'receiver', receiver, '_summary.pdf'), width=25, height=7.5)
        ## ggsave('20696_ligand_summary.pdf', width=25, height=7.5)
    }
}


library(circlize)
pemt_geneset = readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
