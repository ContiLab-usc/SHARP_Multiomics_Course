library(purrr)
library(dplyr)

#Now we make it more generalized, still for K = 2 cluster of each omics layer, 
# but we can have any number of omics layer (2-5) for parallel & serial.

plot_omics_profiles <- function(fit, integration_type, pattern_list, omics_list) {
  list_num = lapply(1:length(omics_list), function(x) paste(x))
  color_p = c("#bf9000", "#a64d79", "#38761d","yellow", "green")
  if(integration_type == "early"){
    M_mean = as.data.frame(fit$res_Mu)
    M_mean$cluster = as.factor(1:2)
    # Reshape the data
    M_mean_melt <- M_mean %>% 
      pivot_longer(cols = -cluster, names_to = "variable", values_to = "value")
    
    M_mean_melt <- M_mean_melt %>% 
      mutate(cluster = ifelse(cluster == 2, "High Risk", "Low Risk"))
    
    # add color label for omics layer
    M_mean_melt = M_mean_melt %>%
      mutate(color_label <- map2(pattern_list, omics_list, 
                  ~ case_when(
                    str_detect(variable, .x) ~ .y,
                    TRUE ~ NA_character_  # Default condition if no pattern matches
                  )
      ) %>%
        
        reduce(coalesce) %>%
        tibble(color_label = .)
        
      ) 
      
    # Filter only the top ## differential expressed features 
    M_mean_top <- M_mean_melt %>% 
      group_by(variable) %>% 
      filter(abs(value) == max(abs(value))) %>% 
      ungroup() %>% 
      arrange(max(abs(value))) %>% 
      group_by(color_label) %>% 
      slice_head(n=12) %>%
      ungroup()

    
    fig <- ggplot(M_mean_melt  %>% filter(variable %in% M_mean_top$variable), 
                  aes(fill = color_label, y = value, x = variable)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle("Omics profiles for the two latent clusters") +
      facet_grid(rows = vars(cluster), cols = vars(color_label),  scales = "free_x", space = "free") +
      theme(legend.position="none") +
      geom_hline(yintercept = 0) +
      xlab("") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle = 90, vjust = 1,
                                       hjust = 1),
            plot.margin = margin(10, 10, 10, 80),
            panel.background = element_rect(fill="white"), 
            strip.background = element_rect(fill = "white"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),) +
      scale_fill_manual(values = color_p[1:length(omics_list)])
    
  } else if(integration_type == "intermediate"){
    rownames(fit$res_Mu[[1]]) = fit$var.names$Znames[[1]]
    rownames(fit$res_Mu[[2]]) = fit$var.names$Znames[[2]]
    rownames(fit$res_Mu[[3]]) = fit$var.names$Znames[[3]]
    
    color_p_p <- c("#bf9000","#38761d","#a64d79","yellow", "green")
    M_mean = map_dfr(fit$res_Mu, ~ as_tibble(.x, rownames = "variable")) %>%
      mutate(omic_layer = map2(pattern_list, omics_list, 
                                 ~ case_when(
                                   str_detect(variable, .x) ~ .y,
                                   TRUE ~ NA_character_  # Default condition if no pattern matches
                                 )
      ) %>%
        reduce(coalesce))
    
    # Reorder results because mirna order is reversed
    M_mean1 <- M_mean %>% 
     #left_join(meta_df, by = c("variable" = "ftr_name")) %>%
     # mutate(`Low Risk`  =  if_else(omic_layer == "RNA", V2, V1), 
             #`High Risk` =  if_else(omic_layer == "RNA", V1, V2)) %>%
      mutate(`Low Risk`  =  V1, 
             `High Risk` =  V2) %>%
      dplyr::select(-c("V1", "V2"))
    
    # Pivot longer for figure 
    M_mean_l <- M_mean1 %>% 
      pivot_longer(cols = c(`Low Risk`, `High Risk`),
                   names_to = "cluster",
                   values_to = "value")
    
    # add color label for omics layer
    M_mean2 = M_mean_l %>%
      mutate(color_label = as.character(as.numeric(factor(omic_layer, levels = unique(omic_layer)))), 
             low_high = if_else(str_detect(cluster, "Low"), 0,1),
             omic = if_else(omic_layer == "miRNA", 
                            "mir",
                            str_sub(omic_layer, end = 1) %>% toupper()),
             omic_cluster = str_c(omic, low_high))
    
    # Filter only the top ## differential expressed features 
    M_mean2_top <- M_mean2 %>% 
      group_by(variable) %>% 
      filter(abs(value) == max(abs(value))) %>% 
      ungroup() %>% 
      arrange(max(abs(value))) %>% 
      group_by(omic_layer) %>% 
      slice_head(n=12) %>%
      ungroup()
    
    # Plots top 12 features
    fig <- ggplot(M_mean2  %>% filter(variable %in% M_mean2_top$variable),
                  aes(fill = color_label, y = value, x = variable)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle("Omics profiles for 2 latent clusters - Lucid in Parallel") +
      facet_grid(rows = vars(cluster),
                 cols = vars(omic_layer), scales = "free_x", space = "free") +
      theme(legend.position="none") +
      geom_hline(yintercept = 0) +
      xlab("") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle = 90, vjust = 1,
                                       hjust = 1),
            plot.margin = margin(10, 10, 10, 80),
            panel.background = element_rect(fill="white"),
            strip.background = element_rect(fill = "white"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),) +
      scale_fill_manual(values = color_p[1:length(omics_list)])
  }
  
  if(integration_type == "late"){
    
    colnames(fit$res_Mu[[1]]) = fit$var.names$Znames[[1]]
    colnames(fit$res_Mu[[2]]) = fit$var.names$Znames[[2]]
    colnames(fit$res_Mu[[3]]) = fit$var.names$Znames[[3]]
    
    fit$submodel[[3]]$var.names$Znames = paste("protein", fit$submodel[[3]]$var.names$Znames, sep = "_")
    colnames(fit$submodel[[3]]$res_Mu) = fit$submodel[[3]]$var.names$Znames
    M_mean = map_dfr(fit$submodel, ~ as_tibble(t(.x$res_Mu), rownames = "variable")) %>%
      mutate(omic_layer = map2(pattern_list, omics_list, 
                               ~ case_when(
                                 str_detect(variable, .x) ~ .y
                               )
      ) %>%
        reduce(coalesce))
    
    # Reorder results because mirna order is reversed
    M_mean1 <- M_mean %>% 
      #left_join(meta_df, by = c("variable" = "ftr_name")) %>%
      mutate(`Low Risk`  =  V1, 
             `High Risk` =  V2) %>%
      dplyr::select(-c("V1", "V2"))
    
    # Pivot longer for figure 
    M_mean_l <- M_mean1 %>% 
      pivot_longer(cols = c(`Low Risk`, `High Risk`),
                   names_to = "cluster",
                   values_to = "value")
    
    # add color label for omics layer
    M_mean2 = M_mean_l %>%
      mutate(color_label = as.character(as.numeric(factor(omic_layer, levels = unique(omic_layer)))),
             low_high = if_else(str_detect(cluster, "Low"), 0,1),
             omic = if_else(omic_layer == "miRNA", 
                            "mir",
                            str_sub(omic_layer, end = 1) %>% toupper()),
             omic_cluster = str_c(omic, low_high))
    
    # Filter only the top ## differential expressed features 
    M_mean2_top <- M_mean2 %>% 
      group_by(variable) %>% 
      filter(abs(value) == max(abs(value))) %>% 
      ungroup() %>% 
      arrange(max(abs(value))) %>% 
      group_by(omic_layer) %>% 
      slice_head(n=12) %>%
      ungroup()
    
    # Plots top 12 features
    fig <- ggplot(M_mean2  %>% filter(variable %in% M_mean2_top$variable),
                  aes(fill = color_label, y = value, x = variable)) +
      geom_bar(position="dodge", stat="identity") +
      ggtitle("Omics profiles for 2 latent clusters - Lucid in Serial") +
      facet_grid(rows = vars(cluster),
                 cols = vars(omic_layer), scales = "free_x", space = "free") +
      theme(legend.position="none") +
      geom_hline(yintercept = 0) +
      xlab("") +
      theme(text = element_text(size=10),
            axis.text.x = element_text(angle = 90, vjust = 1,
                                       hjust = 1),
            plot.margin = margin(10, 10, 10, 80),
            panel.background = element_rect(fill="white"),
            strip.background = element_rect(fill = "white"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),) +
      scale_fill_manual(values = color_p[1:length(omics_list)])
  }
  
  return(fig)
}




