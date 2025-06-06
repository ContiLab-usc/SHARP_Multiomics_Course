
# Get sankey dataframe
get_sankey_df <- function(x,
                          G_color = "dimgray", 
                          X_color = "#eb8c30",
                          Z_color = "#2fa4da", 
                          Y_color = "#afa58e", 
                          pos_link_color = "#67928b", 
                          neg_link_color = "#d1e5eb", 
                          fontsize = 7) {
  K <- x$K
  var.names <- x$var.names
  pars <- x$pars
  dimG <- length(var.names$Gnames)
  dimZ <- length(var.names$Znames)
  valueGtoX <- as.vector(t(x$res_Beta[, -1]))
  valueXtoZ <- as.vector(t(x$res_Mu))
  valueXtoY <- as.vector(x$res_Gamma$beta)[1:K]
  valueXtoY[1] <- 0
  
  # GtoX
  GtoX <- data.frame(
    source = rep(x$var.names$Gnames, K), 
    target = paste0("Latent Cluster", 
                    as.vector(sapply(1:K, function(x) rep(x, dimG)))), 
    value = abs(valueGtoX), 
    group = as.factor(valueGtoX > 0))
  
  # XtoZ
  XtoZ <- data.frame(
    source = paste0("Latent Cluster", 
                    as.vector(sapply(1:K, 
                                     function(x) rep(x, dimZ)))), 
    target = rep(var.names$Znames, 
                 K), value = abs(valueXtoZ),
    group = as.factor(valueXtoZ > 
                        0))
  # XtoY
  XtoY <- data.frame(source = paste0("Latent Cluster", 1:K), 
                     target = rep(var.names$Ynames, K), value = abs(valueXtoY), 
                     group = as.factor(valueXtoY > 0))
  
  links <- rbind(GtoX, XtoZ, XtoY)
  
  nodes <- data.frame(
    name = unique(c(as.character(links$source), 
                    as.character(links$target))), 
    group = as.factor(c(rep("exposure", 
                            dimG), rep("lc", K), rep("biomarker", dimZ), "outcome")))
  
  ## the following two lines were used to exclude covars from the plot #HW added
  links <- links %>% filter(!grepl("cohort", source) & 
                              !grepl("age", source) & 
                              !grepl("fish", source) &
                              !grepl("gender", source))
  nodes <- nodes %>% filter(!grepl("cohort", name) &
                              !grepl("age", name) & 
                              !grepl("fish", name) &
                              !grepl("gender", name))  
  
  
  links$IDsource <- match(links$source, nodes$name) - 1
  links$IDtarget <- match(links$target, nodes$name) - 1
  
  color_scale <- data.frame(
    domain = c("exposure", "lc", "biomarker", 
               "outcome", "TRUE", "FALSE"), 
    range = c(G_color, X_color, 
              Z_color, Y_color, pos_link_color, neg_link_color))
  
  sankey_df = list(links = links, 
                   nodes = nodes)
  
  # p <- sankeyNetwork(
  #   Links = sankey_df$links, 
  #   Nodes = sankey_df$nodes, 
  #   Source = "IDsource", 
  #   Target = "IDtarget",
  #   Value = "value", 
  #   NodeID = "name", 
  #   colourScale = JS(sprintf("d3.scaleOrdinal()\n .domain(%s)\n .range(%s)\n ", 
  #                            jsonlite::toJSON(color_scale$domain), 
  #                            jsonlite::toJSON(color_scale$range))), 
  #   LinkGroup = "group", 
  #   NodeGroup = "group", 
  #   sinksRight = FALSE, 
  #   fontSize = fontsize)
  # p
  return(sankey_df)
}


# sankey_in_serial Function ----
sankey_in_serial <- function(lucid_fit, color_pal_sankey, text_size = 15) {
  
  lucid_fit1 = lucid_fit$submodel[[1]]
  lucid_fit2 = lucid_fit$submodel[[2]]
  lucid_fit3 = lucid_fit$submodel[[3]]
  
  # 1. Get sankey dataframes ----
  sankey_dat1 <- get_sankey_df(lucid_fit1)
  sankey_dat2 <- get_sankey_df(lucid_fit2)
  sankey_dat3 <- get_sankey_df(lucid_fit3)
  
  n_omics_1 <- length(lucid_fit1$var.names$Znames)
  n_omics_2 <- length(lucid_fit2$var.names$Znames)
  n_omics_3 <- length(lucid_fit3$var.names$Znames)
  
  # combine link data
  lnks1_methylation <- sankey_dat1[["links"]] %>% mutate(analysis = "1_methylation")
  lnks2_miRNA  <- sankey_dat2[["links"]] %>% mutate(analysis = "2_miRNA")
  lnks3_RNA <- sankey_dat3[["links"]] %>% mutate(analysis = "3_RNA")
  links <- bind_rows(lnks1_methylation, lnks2_miRNA, lnks3_RNA )
  
  # combine node data
  nodes1_methylation <- sankey_dat1[["nodes"]] %>% mutate(analysis = "1_methylation")
  nodes2_miRNA  <- sankey_dat2[["nodes"]] %>% mutate(analysis = "2_miRNA")
  nodes3_RNA <- sankey_dat3[["nodes"]] %>% mutate(analysis = "3_RNA")
  nodes <- bind_rows(nodes1_methylation, nodes2_miRNA, nodes3_RNA)
  
  
  # 2. Modify analysis 1 ----
  # For analysis 1, latent clusters need to be renamed to names from analysis 2:
  ## 2.1 Get new and original latent cluster names (from the next analysis) ----
  names_clusters_1 <- data.frame(
    name_og = c("Latent Cluster1", "Latent Cluster2"), 
    name_new = c("<b>Serum metabolites\nProfile 0</b>", "<b>Serum metabolites\nProfile 1</b>"))
  
  ## 2.2 Change link names ----
  # Change link names and 
  lnks1_methylation_new <- sankey_dat1[["links"]] %>%
    mutate(
      analysis = "1_methylation",
      source = case_when(
        source == names_clusters_1$name_og[1] ~ names_clusters_1$name_new[1],
        source == names_clusters_1$name_og[2] ~ names_clusters_1$name_new[2],
        TRUE ~ source),
      target = case_when(
        target == names_clusters_1$name_og[1] ~ names_clusters_1$name_new[1],
        target == names_clusters_1$name_og[2] ~ names_clusters_1$name_new[2],
        TRUE ~ target)) %>%
    filter(target != "outcome")
  
  ## 2.3 Change node names ----
  # first, change latent cluster names to analysis specific cluster names
  nodes1_methylation_new <- sankey_dat1[["nodes"]] %>%
    mutate(
      name = case_when(
        name == names_clusters_1$name_og[1] ~ names_clusters_1$name_new[1],
        name == names_clusters_1$name_og[2] ~ names_clusters_1$name_new[2],
        TRUE ~ name), 
      group = if_else(group == "biomarker", "Methylation", as.character(group))) %>%
    filter(group != "outcome")
  
  
  # Visualize
  # sankeyNetwork(
  #   Links = lnks1_methylation_new,
  #   Nodes = nodes1_methylation_new,
  #   Source = "IDsource", Target = "IDtarget",
  #   Value = "value", NodeID = "name", LinkGroup = "group", NodeGroup = "group",
  #   sinksRight = FALSE)
  
  
  # 3. Modify analysis 2 ----
  # For analysis 2, latent clusters need to be renamed to names from analysis 3:
  ## 3.1 Get new and og latent cluster names ----
  names_clusters_2 <- data.frame(
    name_og = c("Latent Cluster1", "Latent Cluster2"), 
    name_new = c("<b>Urine metabolites\nProfile 0</b>", "<b>Urine metabolites\nProfile 1</b>"))
  
  ## 3.2 Change cluster names ----
  lnks2_miRNA_new <- sankey_dat2[["links"]] %>% 
    mutate(
      analysis = "2_miRNA", 
      source = case_when(
        source == names_clusters_2$name_og[1] ~ names_clusters_2$name_new[1], 
        source == names_clusters_2$name_og[2] ~ names_clusters_2$name_new[2], 
        TRUE ~ source), 
      target = case_when(
        target == names_clusters_2$name_og[1] ~ names_clusters_2$name_new[1], 
        target == names_clusters_2$name_og[2] ~ names_clusters_2$name_new[2], 
        TRUE ~ target)) %>%
    filter(target != "outcome")
  
  lnks2_miRNA_new[1,1]  = "<b>Serum metabolites\nProfile 0</b>"
  lnks2_miRNA_new[2,1]  = "<b>Serum metabolites\nProfile 1</b>"
  ## 3.3 Change node names ----
  nodes2_miRNA_new <- sankey_dat2[["nodes"]] %>% 
    mutate(
      name = case_when(
        name == names_clusters_2$name_og[1] ~ names_clusters_2$name_new[1], 
        name == names_clusters_2$name_og[2] ~ names_clusters_2$name_new[2], 
        TRUE ~ name), 
      group = case_when(group == "exposure" ~ "lc", 
                        group == "biomarker" ~ "miRNA",
                        TRUE ~ as.character(group))) %>%
    filter(name != "outcome")
  
  # Visualize
  # sankeyNetwork(
  #   Links = lnks2_transcript_new, 
  #   Nodes = nodes2_transcript_new,
  #   Source = "IDsource", Target = "IDtarget",
  #   Value = "value", NodeID = "name", 
  #   LinkGroup = "group", NodeGroup = "group",
  #   sinksRight = FALSE)
  ##
  
  # 4. Modify analysis 3 ----
  # For analysis 2, latent clusters need to be renamed to names from analysis 3:
  ## 4.1 Get new and og latent cluster names ----
  names_clusters_3 <- tibble(
    name_og = c("Latent Cluster1", "Latent Cluster2"),
    name_new = c("<b>Proteins\nProfile 0</b>", "<b>Proteins\nProfile 1</b>")) 
  
  
  ## 4.2 Change cluster names ----
  lnks3_RNA_new <- sankey_dat3[["links"]] %>% 
    mutate(
      analysis = "3_RNA", 
      source = case_when(
        source == names_clusters_3$name_og[1] ~ names_clusters_3$name_new[1], 
        source == names_clusters_3$name_og[2] ~ names_clusters_3$name_new[2], 
        TRUE ~ source), 
      target = case_when(
        target == names_clusters_3$name_og[1] ~ names_clusters_3$name_new[1], 
        target == names_clusters_3$name_og[2] ~ names_clusters_3$name_new[2], 
        TRUE ~ target))
  
  lnks3_RNA_new[1,1] = "<b>Urine metabolites\nProfile 0</b>"
  lnks3_RNA_new[2,1] = "<b>Urine metabolites\nProfile 1</b>"
  ## 4.3 Change node names ----
  nodes3_RNA_new <- sankey_dat3[["nodes"]] %>% 
    mutate(
      name = case_when(
        name == names_clusters_3$name_og[1] ~ names_clusters_3$name_new[1], 
        name == names_clusters_3$name_og[2] ~ names_clusters_3$name_new[2], 
        TRUE ~ name), 
      group = case_when(group == "exposure" ~ "lc", 
                        group == "biomarker" ~ "RNA",
                        TRUE ~ as.character(group)))
  
  # Test/Visualize
  # sankeyNetwork(
  #   Links = lnks3_protein_new, 
  #   Nodes = nodes3_protein_new,
  #   Source = "IDsource", Target = "IDtarget",
  #   Value = "value", NodeID = "name", LinkGroup = "group", NodeGroup = "group",
  #   sinksRight = FALSE)
  
  
  
  # 5. Combine analysis 1-3 ----
  
  ## 5.1 Final Links ----
  links_all_1 <- bind_rows(lnks1_methylation_new, 
                           lnks2_miRNA_new,
                           lnks3_RNA_new) %>%
    dplyr::select(-IDsource, -IDtarget)
  
  
  ### 5.1.1 Arrange by magnitude ----
  omics_priority <- links_all_1 %>% 
    filter(str_detect(source, "Profile 0"), 
           str_detect(target, "Profile 0", negate = TRUE), 
           str_detect(target, "Profile 1", negate = TRUE), 
           str_detect(target, "outcome", negate = TRUE)) %>%
    group_by(source) %>%
    arrange(desc(group), desc(value), .by_group = TRUE) %>%
    mutate(omics_order = row_number()) %>%
    ungroup() %>%
    dplyr::select(target, omics_order)
  
  
  
  links_all <- links_all_1 
  
  
  
  
  ### 5.1.2 Get new source and target IDs ----
  # First, combine all layers, get unique identifier
  node_ids <- tibble(name = unique(c(unique(links_all$source), 
                                     unique(links_all$target)))) %>%
    mutate(ID = row_number()-1)
  
  # Then combine with original data 
  links_new <- links_all %>%
    left_join(node_ids, by = c("source" = "name")) %>%
    dplyr::rename(IDsource = ID) %>%
    left_join(node_ids, by = c("target" = "name")) %>%
    dplyr::rename(IDtarget = ID)
  
  
  ## 5.2 Final Nodes ----
  nodes_new <- node_ids %>%
    dplyr::select(name) %>%
    left_join(bind_rows(nodes1_methylation_new, 
                        nodes2_miRNA_new,
                        nodes3_RNA_new))
  # remove duplicates 
  nodes_new_nodup <- nodes_new[!base::duplicated(nodes_new),] %>%
    base::as.data.frame()
  
  
  # 6. Plotly Version ----
  library(plotly)
  
  # Add color scheme to nodes
  nodes_new_plotly <- nodes_new_nodup %>% 
    left_join(color_pal_sankey) %>%
    mutate(
      x = case_when(
        group == "exposure" ~ 0,
        str_detect(name, "Methylation") ~ 1/5, 
        str_detect(name, "Methylation") | 
          str_detect(group, "miRNA") ~ 2/5, 
        str_detect(name, "miRNA") | 
          str_detect(group, "RNA") ~ 3/5, 
        str_detect(group, "RNA") ~  4/5, 
        str_detect(group, "outcome") ~ 4.5/5, 
      ))
  
  
  ## 6.2 Get links for Plotly, set color ----
  links_new <- links_new  %>%
    mutate(
      link_color = case_when(
        str_detect(target, "outcome") &  group == TRUE  ~  "red",
        str_detect(target, "outcome") &  group == FALSE  ~  "#d9d2e9",
        
        str_detect(target, "<b>Serum metabolites\nProfile 1</b>") &  group == TRUE  ~  "red",
        str_detect(target, "<b>Serum metabolites\nProfile 1</b>") &  group == FALSE  ~  "#d9d2e9",
        
        str_detect(target, "<b>Urine metabolites\nProfile 1</b>") &  group == TRUE  ~  "red",
        str_detect(target, "<b>Urine metabolites\nProfile 1</b>") &  group == FALSE  ~  "#d9d2e9",
        
        str_detect(target, "<b>Proteins\nProfile 1</b>") &  group == TRUE  ~  "red",
        str_detect(target, "<b>Proteins\nProfile 1</b>") &  group == FALSE  ~  "#d9d2e9",
        
        # Ref link color
        value == 0 ~     "#f3f6f4", 
        # Methylation 
        str_detect(source, "Methylation") &  group == TRUE  ~  "#bf9000",
        str_detect(source, "Methylation") &  group == FALSE ~  "#ffd966",
        # Transcriptome
        str_detect(source, "RNA") &  group == TRUE  ~  "#38761d",
        str_detect(source, "RNA") &  group == FALSE ~  "#b6d7a8",
        # proteome
        str_detect(source, "miRNA") &  group == TRUE  ~  "#a64d79",
        str_detect(source, "miRNA") &  group == FALSE ~  "#ead1dc",
        
        
        
        links_new$group == FALSE ~ "#d9d2e9", # Negative association
        links_new$group == TRUE ~  "red")) # Positive association
  
  plotly_link <- list(
    source = links_new$IDsource,
    target = links_new$IDtarget,
    value = links_new$value+.00000000000000000000001, 
    color = links_new$link_color)  
  
  
  # Get list of nodes for Plotly
  plotly_node <- list(
    label = nodes_new_plotly$name, 
    color = nodes_new_plotly$color,
    pad = 15,
    thickness = 20,
    line = list(color = "black",width = 0.5), 
    x = nodes_new_plotly$x, 
    y = c(0.01, 
          0.1, 0.3, # Methylation clusters
          .45, .55, # Transcriptome clusters
          .80, .95, # Proteome clusters
          seq(from = .01, to = 1, by = 0.035)[1:n_omics_1], # Cpgs (10 total)
          seq(from = 0.35, to = 1, by = 0.025)[1:n_omics_2], # miRNA (8 total)
          seq(from = 0.75, to = 1, by = 0.03)[1:n_omics_3], # Transcript (10 total)
          .95
    ))
  
  
  ## 6.3 Plot Figure ----
  fig <- plot_ly(
    type = "sankey",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)),
    orientation = "h",
    node = plotly_node,
    link = plotly_link)
  
  (fig <- fig %>% layout(
    # title = "Basic Sankey Diagram",
    font = list(
      size = text_size
    ))
  )
  
  return(fig)
}

col_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
color_pal_sankey_serial <- matrix(c("exposure", col_pal[6],
                                    "lc",       "#b3d8ff",
                                    "RNA",      col_pal[1],
                                    "Methylation",       col_pal[2],
                                    "miRNA",  col_pal[3],
                                    "outcome",  "grey"), 
                                  ncol = 2, byrow = TRUE) %>%
  as_tibble(.name_repair = "unique") %>% 
  janitor::clean_names() %>%
  dplyr::rename(group = x1, color = x2)
