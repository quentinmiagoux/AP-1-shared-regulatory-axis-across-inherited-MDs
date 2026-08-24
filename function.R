
extract_column <- function(list, db, method){
  
  b <- a %>% map(~ .$enrichment[[method]]) %>% map(~ .[[db]] ) %>% map(~  .$Description)
  return(b)
}
extract_degs <- function(list,updown,logtresh){
  
  b <- list %>% map(~ filter(.$dataset, group == updown,  abs(log2FoldChange) > log2(logtresh), padj < 0.05))  %>% map(~  unique(.$Symbol))
  return(b)
}

extract_degs_no <- function(list,logtresh){
  
  b <- list %>% 
    map(~ filter(.$dataset,  abs(log2FoldChange) > log2(logtresh), padj < 0.05)) %>% 
    map(~ mutate(.,  id = paste0(Symbol, "_",group))) %>%
    map(~  unique(.$id))
  return(b)
}

extract_degs_no_up_down <- function(list,logtresh){
  b <- list %>% map(~ filter(.$dataset,  abs(log2FoldChange) > log2(logtresh), padj < 0.05)) %>% 
    map(~ mutate(., Symbol = if_else(log2FoldChange>0, paste0(Symbol,"_","up"), paste0(Symbol,"_","down")))) %>%
    map(~  unique(.$Symbol))
  return(b)
}

extract_degs_egi <- function(list,updown,logtresh){
  b <- list %>% 
    map(~ dplyr::filter(.$dataset, group == updown,  abs(log2FoldChange) > log2(logtresh), padj < 0.05))  %>% 
    map(~  unique(.$entrezgene_id))
  return(b)
}
extract_degs_egi_no_diff <- function(list,logtresh){
  
  b <- list %>% 
    map(~ filter(.,  abs(log2FoldChange) > log2(logtresh), padj < 0.05))  %>% 
    map(~  unique(.$Symbol))
  return(b)
}
split_labels <- function(labels, max.char = 30) {
  sapply(labels, function(label) {
    if (nchar(label) <= max.char) return(label)
    
    last_space_before_max <- max.char - str_length(str_extract(substr(label, 1, max.char), "\\s+\\S*$"))
    
    first_part <- substr(label, 1, last_space_before_max)
    second_part <- substr(label, last_space_before_max + 1, nchar(label))
    
    paste0(first_part, "\n", str_trim(second_part))
  }, USE.NAMES = FALSE)
}

heatmap_degs_pathway <- function(pathway){
  
  
  DEGS <- kegg %>% filter(ID == pathway) %>% .$geneID %>% str_split("/") %>% unlist
  
  
  datasets_to_filt <- datasets %>%
    map(~ filter(.$dataset, Symbol %in% DEGS, padj<0.05, abs(log2FoldChange)>log2(1.5))) %>%
    map(~ arrange(., padj)) %>%
    map(~ distinct(., Symbol, .keep_all = TRUE)) %>% 
    .[1:4]
  
  fc_list <- lapply(datasets_to_filt, function(df) {
    # Sélectionne uniquement les deux colonnes nécessaires
    fc <- df[, c("Symbol", "log2FoldChange")]
    # Crée un vecteur nommé (noms = ID, valeurs = log2FC)
    setNames(fc$log2FoldChange, fc$Symbol)
  })
  
  
  # Étape 1 : extraire tous les gènes uniques
  all_genes <- unique(unlist(lapply(fc_list, names)))
  
  # Étape 2 : réindexer chaque vecteur pour avoir les mêmes gènes dans le même ordre
  fc_list_filled <- lapply(fc_list, function(vec) {
    out <- setNames(rep(NA, length(all_genes)), all_genes)
    out[names(vec)] <- vec
    return(out)
  })
  
  # Étape 3 : combiner en data.frame (lignes = conditions, colonnes = gènes)
  fc_mat <- as.data.frame(do.call(rbind, fc_list_filled))
  
  # Étape 4 : (optionnel) transposer pour avoir gènes en lignes, conditions en colonnes
  fc_mat_t <- t(fc_mat)
  
  fc_mat_t_with_n <- as.data.frame(fc_mat_t)
  fc_mat_t_with_n$n <- rowSums(!(fc_mat_t_with_n==0 | is.na(fc_mat_t_with_n)))
  # Résultat : matrice avec gènes en ligne et maladies en colonne
  fc_mat_t_with_n3 <-fc_mat_t_with_n %>% filter(n>=3)
  
  
  fc_mat_t_with_n3[is.na(fc_mat_t_with_n3)] <- 0
  
  # Palette de couleur centrée sur 0 (log2FC)
  col_fun <- circlize::colorRamp2(c(-3, 0, 3), c("blue", "grey90", "red"))
  
  
  # Heatmap de base
  # Définir style des noms :
  row_font <- gpar(fontsize = 10)        # gènes en gras
  col_font <- gpar(fontface = "bold", fontsize = 10)     # maladies inclinées
  
  # Heatmap
  heatmap_lysosomes <- Heatmap(
    fc_mat_t_with_n3[1:4],
    name = "log2FC",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = row_font,        # gènes en gras
    column_names_gp = col_font,
    na_col = "grey90",
    column_names_rot = 45,  # ← rotation à 90° des noms de colonnes
    heatmap_legend_param = list(
      title = "log2FC",
      color_bar = "continuous",
      legend_height = unit(4, "cm")
    )
  )
  return(heatmap_lysosomes)
  
}

datatable_fun <- function(x){x %>% mutate_if(is.numeric, ~format(., scientific = T, digits = 2)) %>% datatable(
  options = list(
    pageLength = 10,
    scrollX = TRUE,
    autoWidth = TRUE,
    columnDefs = list(list(className = "dt-nowrap", targets = "_all"))
  ),
  class = "nowrap",
  filter = 'top',
  rownames = F,
  escape = F
)
}


arrange_it <- function(ck_up){
  
  shared_terms <- ck_up %>%
    group_by(Description) %>%
    summarise(n_clusters = n_distinct(Cluster),mRichFactor =  mean(RichFactor)) %>%
    filter(n_clusters >= 3) %>%
    arrange(desc(mRichFactor))
  
  
  ck_df_filtered <- ck_up %>%
    filter(Description %in% shared_terms$Description)
  
  ck_df_filtered@compareClusterResult$Description <- factor(
    ck_df_filtered@compareClusterResult$Description,
    levels = shared_terms$Description
  )
  
  ck_df_filtered@compareClusterResult <- ck_df_filtered@compareClusterResult %>%
    arrange(factor(Description, levels = shared_terms$Description))
  
  return(ck_df_filtered)
  
}

# Calculer les valeurs minimales et maximales de RichFactor pour harmoniser les tailles
get_richfactor_limits <- function(ck) {
  df <- as.data.frame(ck)
  rf <- mapply(function(count, bgratio) {
    count / as.numeric(strsplit(bgratio, "/")[[1]][1])
  }, df$Count, df$BgRatio)
  range(rf, na.rm = TRUE)
}


filter_DEGs <- function(df) {
  df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > log2(1.5)) %>% arrange(padj) %>% distinct(Symbol, .keep_all = T) %>%
    dplyr::select(Symbol, log2FoldChange)
}

enrichplot_point_shape <- ggfun:::enrichplot_point_shape


plot_gsea_dotplot <- function(
    object,
    x = "Cluster",
    colorBy = "p.adjust",
    showCategory = 5,
    by = "geneRatio",
    size = "geneRatio",
    split = NULL,
    includeAll = TRUE,
    font.size = 12,
    title = "",
    label_format = 30,
    group = FALSE,
    shape = FALSE,
    facet = NULL,
    strip_width = 15
) {
  color <- NULL
  if (is.null(size)) {
    size <- by
  } ## by may deprecated in future release
  
  if (!is.null(facet) && facet == "intersect") {
    object <- append_intersect(object)
  }
  
  df <- fortify(
    object,
    showCategory = showCategory,
    by = size,
    includeAll = includeAll,
    split = split
  )
  
  # if (by != "geneRatio")
  #    df$GeneRatio <- parse_ratio(df$GeneRatio)
  label_func <- default_labeller(label_format)
  if (is.function(label_format)) {
    label_func <- label_format
  }
  
  if (size %in% c("rowPercentage", "count", "geneRatio")) {
    by2 <- switch(
      size,
      rowPercentage = "Percentage",
      count = "Count",
      geneRatio = "GeneRatio"
    )
  } else {
    by2 <- size
  }
  
  p <- ggplot(
    df,
    aes(x = .data[[x]], y = .data[["Description"]], size = .data[[by2]])
  ) +
    scale_y_discrete(labels = label_func)
  
  if (group) {
    p <- p +
      geom_line(
        aes(color = .data$Cluster, group = .data$Cluster),
        size = .3
      ) +
      ggnewscale::new_scale_colour()
  }
  
  if (shape) {
    check_installed('ggstar', 'for `dotplot()` with `shape = TRUE`.')
    ggstar <- "ggstar"
    require(ggstar, character.only = TRUE)
    # p <- p + ggsymbol::geom_symbol(aes(symbol = .data$Cluster, fill = .data[[colorBy]])) +
    p <- p +
      ggstar::geom_star(aes(
        starshape = .data$Cluster,
        fill = .data[[colorBy]]
      )) +
      set_enrichplot_color(type = "fill", transform = 'log10')
  } else {
    p <- p +
      geom_point(aes(fill = .data[[colorBy]])) +
      aes(shape = I(enrichplot_point_shape))
  }
  
  p <- p +
    set_enrichplot_color(type = "fill", transform = 'log10') +
    ylab(NULL) +
    ggtitle(title) +
    DOSE::theme_dose(font.size) +
    scale_size_continuous(range = c(3, 8))
  
  if (!is.null(facet)) {
    p <- p +
      facet_grid(
        .data[[facet]] ~ .,
        scales = "free",
        space = 'free',
        switch = 'x',
        labeller = ggplot2::label_wrap_gen(strip_width)
      ) +
      theme(strip.text = element_text(size = 14))
  }
  
  class(p) <- c("enrichplotDot", class(p))
  
  return(p)
}

default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

ep_str_wrap <- function(string, width) {
  x <- gregexpr(' ', string)
  vapply(seq_along(x),
         FUN = function(i) {
           y <- x[[i]]
           n <- nchar(string[i])
           len <- (c(y,n) - c(0, y)) ## length + 1
           idx <- len > width
           j <- which(!idx)
           if (length(j) && max(j) == length(len)) {
             j <- j[-length(j)]
           }
           if (length(j)) {
             idx[j] <- len[j] + len[j+1] > width
           }
           idx <- idx[-length(idx)] ## length - 1
           start <- c(1, y[idx] + 1)
           end <- c(y[idx] - 1, n)
           words <- substring(string[i], start, end)
           paste0(words, collapse="\n")
         },
         FUN.VALUE = character(1)
  )
}


make_pathway_heatmap <- function(
    df_gsea,
    datasets,
    pathway = "Lysosome",
    dataset_names = NULL,
    mode = c("filtered", "all"),
    fc_filter = 0.5,
    use_padj = FALSE,
    padj_cutoff = 0.05,
    min_presence_datasets = 1,
    fc_limits = c(-2, 0, 2),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    na_col = "grey90"
) {
  mode <- match.arg(mode)
  
  # --- checks ----------------------------------------------------------------
  required_gsea_cols <- c("Description", "core_enrichment")
  if (!all(required_gsea_cols %in% colnames(df_gsea))) {
    stop(
      "df_gsea must contain columns: ",
      paste(required_gsea_cols, collapse = ", ")
    )
  }
  
  if (!is.list(datasets) || length(datasets) == 0) {
    stop("'datasets' must be a non-empty list.")
  }
  
  # --- 1. pathway genes ------------------------------------------------------
  pathway_genes <- df_gsea %>%
    dplyr::filter(.data$Description == pathway) %>%
    dplyr::pull(.data$core_enrichment) %>%
    stringr::str_split("/") %>%
    unlist() %>%
    unique()
  
  if (length(pathway_genes) == 0) {
    stop("No genes found for pathway: ", pathway)
  }
  
  pathway_genes <- unique(pathway_genes)
  
  # --- 2. subset datasets ----------------------------------------------------
  if (!is.null(dataset_names)) {
    missing_datasets <- setdiff(dataset_names, names(datasets))
    if (length(missing_datasets) > 0) {
      stop(
        "These dataset_names are missing from 'datasets': ",
        paste(missing_datasets, collapse = ", ")
      )
    }
    datasets <- datasets[dataset_names]
  }
  
  # --- 3. clean datasets -----------------------------------------------------
  datasets_clean <- purrr::map(
    datasets,
    function(x) {
      if (!"dataset" %in% names(x)) {
        stop("Each element of 'datasets' must contain a $dataset data frame.")
      }
      
      df <- x$dataset
      
      # important: remove grouping
      df <- dplyr::ungroup(df)
      
      # if your gene column is mgi_symbol instead of Symbol, rename it
      if (!"Symbol" %in% colnames(df) && "mgi_symbol" %in% colnames(df)) {
        df <- dplyr::rename(df, Symbol = .data$mgi_symbol)
      }
      
      required_dataset_cols <- c("Symbol", "log2FoldChange", "padj")
      if (!all(required_dataset_cols %in% colnames(df))) {
        stop(
          "Each datasets[[i]]$dataset must contain columns: ",
          paste(required_dataset_cols, collapse = ", ")
        )
      }
      
      df %>%
        dplyr::filter(.data$Symbol %in% pathway_genes) %>%
        dplyr::arrange(.data$padj) %>%
        dplyr::distinct(.data$Symbol, .keep_all = TRUE) %>%
        dplyr::select(.data$Symbol, .data$log2FoldChange, .data$padj)
    }
  )
  
  # --- 4. gene presence ------------------------------------------------------
  gene_counts <- table(unlist(purrr::map(datasets_clean, ~ .x$Symbol)))
  genes_keep <- names(gene_counts[gene_counts >= min_presence_datasets])
  
  selected_genes <- pathway_genes[pathway_genes %in% genes_keep]
  
  if (length(selected_genes) == 0) {
    stop(
      "No genes remain after applying min_presence_datasets = ",
      min_presence_datasets
    )
  }
  
  # --- 5. build matrix -------------------------------------------------------
  if (mode == "all") {
    fc_matrix <- sapply(
      datasets_clean,
      function(df) {
        as.numeric(df$log2FoldChange[match(selected_genes, df$Symbol)])
      },
      simplify = "matrix"
    )
    
    rownames(fc_matrix) <- selected_genes
    colnames(fc_matrix) <- names(datasets_clean)
  }
  
  if (mode == "filtered") {
    filtered_list <- purrr::map(
      datasets_clean,
      function(df) {
        keep <- !is.na(df$log2FoldChange) & abs(df$log2FoldChange) >= fc_filter
        
        if (use_padj) {
          keep <- keep & !is.na(df$padj) & df$padj < padj_cutoff
        }
        
        out <- rep(NA_real_, length(selected_genes))
        names(out) <- selected_genes
        
        idx <- match(df$Symbol, selected_genes)
        
        # keep only valid matches
        valid <- keep & !is.na(idx)
        
        if (any(valid)) {
          out[idx[valid]] <- as.numeric(df$log2FoldChange[valid])
        }
        
        out
      }
    )
    
    fc_matrix <- do.call(cbind, filtered_list)
    rownames(fc_matrix) <- selected_genes
    colnames(fc_matrix) <- names(filtered_list)
    
    fc_matrix <- fc_matrix[
      rowSums(!is.na(fc_matrix)) >= min_presence_datasets,
      ,
      drop = FALSE
    ]
    
    if (nrow(fc_matrix) == 0) {
      stop(
        "No genes remain after filtering. ",
        "Try lowering fc_filter, padj_cutoff, or min_presence_datasets."
      )
    }
  }
  
  # --- 6. cap values ---------------------------------------------------------
  fc_matrix[fc_matrix < min(fc_limits, na.rm = TRUE)] <- min(fc_limits, na.rm = TRUE)
  fc_matrix[fc_matrix > max(fc_limits, na.rm = TRUE)] <- max(fc_limits, na.rm = TRUE)
  
  fc_matrix_clust <- fc_matrix
  fc_matrix_clust[is.na(fc_matrix_clust)] <- 0
  
  # --- 7. colors -------------------------------------------------------------
  col_fun <- circlize::colorRamp2(
    fc_limits,
    c("blue", "#FFFFFF", "red")
  )
  
  # --- 8. title --------------------------------------------------------------
  title_suffix <- if (mode == "filtered") {
    paste0(
      " | abs(log2FC) >= ", fc_filter,
      if (use_padj) paste0(" & padj < ", padj_cutoff) else "",
      " | min datasets = ", min_presence_datasets
    )
  } else {
    paste0(" | all log2FC | min datasets = ", min_presence_datasets)
  }
  
  # --- 9. heatmap ------------------------------------------------------------
  heatmap_obj <- ComplexHeatmap::Heatmap(
    fc_matrix,
    name = "log2FC",
    col = col_fun,
    na_col = na_col,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    clustering_distance_rows = function(x) stats::dist(fc_matrix_clust),
    clustering_distance_columns = function(x) stats::dist(t(fc_matrix_clust)),
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontface = "bold", fontsize = 8),
    column_names_gp = grid::gpar(fontface = "bold", fontsize = 10),
    column_names_rot = 45,
    heatmap_legend_param = list(
      title = "log2FC",
      at = fc_limits,
      labels = fc_limits,
      color_bar = "continuous",
      legend_height = grid::unit(4, "cm")
    ),
    column_title = paste0(pathway, title_suffix)
  )
  
  return(list(
    pathway_genes = pathway_genes,
    selected_genes = rownames(fc_matrix),
    datasets_clean = datasets_clean,
    fc_matrix = fc_matrix,
    fc_matrix_clust = fc_matrix_clust,
    heatmap = heatmap_obj
  ))
}
