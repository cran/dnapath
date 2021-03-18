#' @importFrom utils globalVariables
utils::globalVariables(c("dnapath.mart_species", 
                         "dnapath.mart_obj", 
                         "entrezgene_id", 
                         "symbol"))

#' Obtain Reactome pathways
#' 
#' Connects to `reactome.db` \insertCite{reactome}{dnapath}  
#' to obtain a list of pathways for a given species.
#' The pathway list is processed by combining any two pathways that have 
#' substantial overlap (default is over 90% overlap). This output if this 
#' function can be used for the `pathway_list` argument in 
#' \code{\link{dnapath}}.
#'
#' @param species A string, for example "Homo sapiens" or "Mus musculus", 
#' indicating the species to use. 
#' @param overlap_limit (Optional) Any pathways that have an overlap
#' greater than overlap_limit are combined. Set to NULL to disable
#' this option.
#' @param min_size The minimum pathway size. Any Reactome pathways with fewer
#' than `min_size` genes are removed from the list. Defaults to 10. 
#' @param max_size The maximum pathway size. Any Reactome pathways with more
#' than `max_size` genes are removed from the list. Defaults to 50. 
#' @param verbose Set to FALSE to turn off messages.
#' @return A named list of vectors. Each vector corresponds to a Reactome pathway
#' and contains the entrezgene IDs of the genes in that pathway.
#' @references 
#' \insertRef{reactome}{dnapath}
#' @seealso 
#' The genes in the Reactome pathways use entrezgene IDs. These can be converted
#' to gene symbols, if desired, using the \code{\link{entrez_to_symbol}} and
#' \code{\link{rename_genes}} functions.
#' @export
#' @examples 
#' # Obtaining a pathway list for human (Homo sapiens).
#' # In this example, overlapping pathways are not combined (this is
#' # specified by setting overlap_limit to NULL).
#' pathway_list <- get_reactome_pathways("Homo sapiens", overlap_limit = NULL,
#'                                       min_size = 10, max_size = 20)
get_reactome_pathways <- function(species, overlap_limit = 0.9, min_size = 10, 
                                  max_size = 50, verbose = TRUE) {
  if(!requireNamespace("reactome.db", quietly = TRUE)) {
    message("Warning: The `reactome.db` package must be installed to use get_reactome_pathways().",
            "Returning the p53_pathway list by default.")
    return(dnapath::p53_pathways)
  }
  
  # If species is "human" or "mouse", change to appropriate name.
  if(tolower(species) == "human") 
    species <- "Homo sapiens"
  if(tolower(species) == "mouse" || tolower(species) == "mice") 
    species <- "Mus musculus"
  
  if(verbose)
    cat("Obtaining reactome pathway information for species:", species, "\n")
  # Get list of reactome pathway names
  # Map names to ID
  # Map ID to entrez ID
  # Save results as p by g matrix. (Optional; not performed.)
  
  # library(reactome.db)
  # Map reactome NAME to reactome ID. Subset on pathways for this species.
  reactome_to_id <- as.list(reactome.db::reactomePATHNAME2ID)
  index <- which(tolower(species) == 
                   tolower(gsub(":.*", "", 
                                names(as.list(reactome.db::reactomePATHNAME2ID)))))
  if(length(index) == 0) {
    available_species <- unique(
      gsub("(: .*)", "", names(as.list(reactome.db::reactomePATHNAME2ID))))
    stop(species, " is not an available species. Use one of the following:\n\t", 
         paste(available_species[-length(available_species)], collapse = ",\n\t"), 
         available_species[length(available_species)], ".")
  }
  reactome_to_id <- reactome_to_id[index]
  
  # Map reactome ID to entrez ID.
  id_to_entrez <- as.list(reactome.db::reactomePATHID2EXTID)
  
  # Subset on reactome IDs obtained in first step.
  index <- which(names(id_to_entrez) %in% sapply(reactome_to_id, dplyr::first))
  id_to_entrez <- id_to_entrez[index]
  
  # Remove any pathways containing fewer than min_size genes or more than 
  # max_size genes.
  sizes <- sapply(id_to_entrez, function(x) length(unique(x)))
  id_to_entrez <- id_to_entrez[sizes >= min_size & sizes <= max_size]
  
  # Map reactome ID back to reactome NAME.
  id_to_reactome <- as.list(reactome.db::reactomePATHID2NAME)
  
  # Subset on reactome IDs obtained in previous step.
  index <- which(names(id_to_reactome) %in% names(id_to_entrez))
  id_to_reactome <- id_to_reactome[index]
  
  # Finally, map reactome NAME to entrez ID
  pathway_list <- id_to_entrez
  names(pathway_list) <- sapply(names(id_to_entrez), function(x) {
    dplyr::first(id_to_reactome[[which(names(id_to_reactome) == x)[1]]])
  })
  names(pathway_list) <- gsub("^[ A-Za-z]*: ", "", names(pathway_list))
  
  if(!is.null(overlap_limit)) {
    if(((overlap_limit >= 1) || (overlap_limit <= 0))) {
      warning("`overlap_limit` is not between 0 and 1. Overlapping pathways",
              " are not combined.")
    } else {
      if(verbose)
        cat("Combining pathways that have greater than", overlap_limit * 100,
            "% overlap.\n")
      pathway_list <- combine_overlapping_pathways(pathway_list, overlap_limit)
    }
  }
  
  return(pathway_list)
}

#' Modify a pathway list to combine overlapping pathways.
#'
#' @param pathway_list A list of pathways obtained from
#' \code{\link{get_reactome_pathways}}.
#' @param overlap_limit A percentage between 0 and 1. If two pathways
#' overlap by more than this amount, they are combined into one pathway.
#' @return A modified list with overlapping pathways combined together.
#' @keywords internal
combine_overlapping_pathways <- function(pathway_list, overlap_limit = 0.9) {
  if((overlap_limit >= 1) || (overlap_limit <= 0)) {
    warning("`overlap_limit` is not between 0 and 1. Returning pathway list",
            " unchanged.")
    return(pathway_list)
  }
  
  # Determine which pathways are similar.
  n <- length(pathway_list)
  similar_to <- diag(FALSE, n)
  for(i in 2:n) {
    for(j in 1:(i - 1)) {
      len_intersect <- length(intersect(pathway_list[[i]], pathway_list[[j]]))
      if(len_intersect == 0) {
        similar_to[i, j] <- FALSE
      } else {
        len_union <- length(union(pathway_list[[i]], pathway_list[[j]]))
        if(len_union != 0)
          similar_to[i, j] <- (len_intersect / len_union) > overlap_limit
      }
    }
  }
  pathways_to_remove <- NULL
  for(i in 1:ncol(similar_to)) {
    if(sum(similar_to[, i]) > 0) {
      index <- which(similar_to[, i])
      names(pathway_list)[i] <- paste0(names(pathway_list)[i], " (See also: ",
                                       paste0(names(pathway_list)[index],
                                              collapse = "; "), ")")
      
      pathway_list[[i]] <- unique(unlist(pathway_list[c(i, index)]))
      pathways_to_remove <- c(pathways_to_remove, index)
    }
  }
  
  pathways_to_remove <- unique(pathways_to_remove)
  
  if(length(pathways_to_remove) > 0) {
    pathway_list <- pathway_list[-pathways_to_remove]
  }
  return(pathway_list)
}







#' Obtain gene symbols for entrezgene IDs
#'
#' Uses `biomaRt` \insertCite{biomart}{dnapath} to map entrezgene IDs to gene 
#' symbols for a given species. Obtains MGI symbols for mouse species and
#' HGNC symbols for other species.
#' (Note: this mapping may not work for all species.)
#' The output of this function can be used in \code{\link{rename_genes}}. 
#' 
#' If entrezgene IDs are used in a `dnapath_list` or `dnapath`
#' object, or a pathway list, then \code{\link{get_genes}} can be used to 
#' extract them and used for the `x` argument here.
#' 
#' @param x A vector of entrezgene IDs.
#' @param species The species used to obtain the entrezgene IDs. For example:
#' "Homo sapiens", "m musculus", "C. elegans", or "S cerevisiae".
#' "Human" and "mouse" can also be used and will be converted to the 
#' correct species name.
#' @param dir_save The directory to store annotation reference. Future
#' calls to this function will use the stored annotations. This speeds up the
#' operation and allows for reproducibility in the event that the `biomaRt`
#' database is updated. Set to NULL to disable. By default, it uses a
#' temporary directory to store files during the R session.
#' @param verbose Set to FALSE to avoid messages.
#' @return A data frame with two columns: the first contains the original
#' entrezgene IDs, and the second contains the corresponding gene symbols.
#' MGI symbols are returned when `species = "Mus musculus"` and HGNC symbols
#' are returned otherwise.
#' @note
#' Internet connection is required to connect to biomaRt. If unavailable, the
#' default biomart and default species contained in the package is used, but
#' this may not match the desired species.
#' @references 
#' \insertRef{biomart}{dnapath}
#' @seealso 
#' \code{\link{symbol_to_entrez}}, \code{\link{get_genes}}
#' @export
#' @examples 
#' \donttest{
#' data(meso)
#' # The meso gene expression data contains entrezgene IDs. 
#' # These can be converted to gene symbols.
#' gene_mat <- entrez_to_symbol(colnames(meso$gene_expression), species = "human")
#' }
entrez_to_symbol <- function(x, 
                             species,
                             dir_save = tempdir(),
                             verbose = TRUE) {
  
  if(class(x) %in% c("dnapath", "dnapath_list")) {
    stop("Input should be a vector of entrezgene IDs, not a dnapath object. ",
         'Use rename_genes() with `to = "symbol"`')
  }
  
  species <- format_species_name(species)
  
  attribute <- "hgnc_symbol"
  if(species == "mmusculus") {
    attribute <- "mgi_symbol"
  }
  
  load_file <- NULL
  if(!is.null(dir_save))
    load_file <- file.path(dir_save, paste0("entrez_to_", species, ".rds"))
    
  if(!is.null(dir_save) && file.exists(load_file)) {
    if(verbose) cat("\t- loading gene info from", load_file, "\n")
    gene_info <- readRDS(load_file)
  } else {
    if(!all(c("dnapath.mart_obj", "dnapath.mart_species") %in% ls()) || 
       is.null(dnapath.mart_obj) || is.null(dnapath.mart_species) ||
       (dnapath.mart_species != species)) 
      init_mart(species)
    
    gene_info <- biomaRt::getBM(attributes = c("entrezgene_id", attribute),
                                mart = dnapath.mart_obj)
    gene_info %>%
      dplyr::filter(!is.na(entrezgene_id)) %>%
      dplyr::group_by(entrezgene_id) %>%
      dplyr::summarise(across(everything(), dplyr::first)) ->
      gene_info
    
    # Store the annotations for future reference.
    if(!is.null(dir_save)) {
      save_file <- load_file
      if(verbose)
        cat("\t- saving gene info to", save_file, "\n")
      if(!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
      saveRDS(gene_info, save_file)
    }
  }

  # If the entrezgene IDs are a factor or character object, coerce into numeric.
  if(!is.numeric(x)) {
    if(is.factor(x)) {
      x <- as.numeric(as.character(x))
    } else {
      x <- as.numeric(x)
    }
  }
  
  df <- data.frame(entrezgene_id = x, 
                   stringsAsFactors = FALSE)
  df <- dplyr::left_join(df, gene_info, by = "entrezgene_id")
  df <- dplyr::distinct(df)
  
  # If any symbols were not found, replace NA value(s) with the entrezgene ID.
  index_na <- which(is.na(df[, 2]))
  if(length(index_na) > 0) {
    df[index_na, 2] <- df[index_na, 1]
  }
  
  return(df)
}




#' Obtain entrezgene IDs for gene symbols
#'
#' Uses `biomaRt` \insertCite{biomart}{dnapath}   
#' to map entrezgene IDs to gene symbols for a given species. The output of
#' this function can be used in \code{\link{rename_genes}}.
#' 
#' If entrezgene IDs are used in a `dnapath_list` or `dnapath`
#' object, or a pathway list, then \code{\link{get_genes}} can be used to 
#' extract them and used for the `x` argument here.
#' 
#' @param x A vector of gene symbols.
#' @param species The species used to obtain the entrezgene IDs. For example:
#' "Homo sapiens", "m musculus", "C. elegans", or "S cerevisiae".
#' "Human" and "mouse" can also be used and will be converted to the 
#' correct species name.
#' @param dir_save The directory to store annotation reference. Future
#' calls to this function will use the stored annotations. This speeds up the
#' operation and allows for reproducibility in the event that the `biomaRt`
#' database is updated. Set to NULL to disable. By default, it uses a
#' temporary directory to store files during the R session.
#' @param verbose Set to FALSE to avoid messages.
#' @return A data frame with two columns: the first contains the original
#' gene symbols, and the second contains a corresponding entrezgene ID. 
#' @note
#' Internet connection is required to connect to biomaRt. If unavailable, the
#' default biomart and default species contained in the package is used, but
#' this may not match the desired species.
#' 
#' It is assumed that `x` contains MGI symbols when the biomart species is
#' "Mus musculus" and HGNC symbols otherwise.
#' @references 
#' \insertRef{biomart}{dnapath}
#' @seealso 
#' \code{\link{entrez_to_symbol}}, \code{\link{get_genes}} 
#' @export
#' @examples
#' \donttest{
#' # Convert a set of gene symbols to entrezgene IDs.
#' # Note that not all may have mapping (such as "MSX" in this example).
#' gene_mat <- symbol_to_entrez(c("SOX2", "SEMA3E", "COL11A1", "UBB", "MSX"),
#'                              species = "human")
#' }
symbol_to_entrez <- function(x, 
                             species,
                             dir_save = tempdir(),
                             verbose = TRUE) {
  
  if(class(x) %in% c("dnapath", "dnapath_list")) {
    stop("Input should be a vector of entrezgene IDs, not a dnapath object.")
  }
  
  species <- format_species_name(species)
  
  attribute <- "hgnc_symbol"
  if(species == "mmusculus") {
    attribute <- "mgi_symbol"
  }
  
  load_file <- NULL
  if(!is.null(dir_save))
    load_file <- file.path(dir_save, paste0("entrez_to_", species, ".rds"))
  
  if(!is.null(dir_save) && file.exists(load_file)) {
    if(verbose) cat("\t- loading gene info from", load_file, "\n")
    gene_info <- readRDS(load_file)
  } else {
    if(!all(c("dnapath.mart_obj", "dnapath.mart_species") %in% ls()) || 
       is.null(dnapath.mart_obj) || is.null(dnapath.mart_species) ||
       (dnapath.mart_species != species)) 
      init_mart(species)
    
    gene_info <- biomaRt::getBM(attributes = c("entrezgene_id", attribute),
                                mart = dnapath.mart_obj)
    
    # The order of columns is not fixed (sometimes entrezgene IDs are second).
    # Identify which column contains the symbols.
    index_symbol_column <- which(colnames(gene_info) == attribute)[1]
    colnames(gene_info)[index_symbol_column] <- "symbol"
    gene_info %>%
      dplyr::filter(!is.na(entrezgene_id)) %>%
      dplyr::group_by(symbol) %>%
      dplyr::summarise(across(everything(), dplyr::first)) %>%
      dplyr::ungroup() ->
      gene_info
    
    # The grouping variable (gene symbol) will now be the first column
    colnames(gene_info)[1] <- attribute
    
    # Store the annotations for future reference.
    if(!is.null(dir_save)) {
      save_file <- load_file
      if(verbose)
        cat("\t- saving gene info to", save_file, "\n")
      if(!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
      saveRDS(gene_info, save_file)
    }
  }
  
  # If the gene symbols are a factor, coerce into character.
  if(is.factor(x)) {
    x <- as.character(x)
  }
  
  df <- data.frame(symbol = x, 
                   stringsAsFactors = FALSE)
  colnames(df) <- attribute
  df <- dplyr::left_join(df, gene_info, by = attribute)
  df <- dplyr::distinct(df)
  
  # If any symbols were not found, replace NA value(s) with the entrezgene ID.
  index_na <- which(is.na(df[, 2]))
  if(length(index_na) > 0) {
    df[index_na, 2] <- df[index_na, 1]
  }
  
  return(df)
}



#' Initialize biomaRt for a given species
#'
#' @param species The species to obtain a biomart dataset for. For example:
#' "Homo sapiens", "m musculus", "C. elegans", or "S cerevisiae".
#' "Human" and "mouse" can also be used and will be converted to the 
#' correct species name.
#' @return None. Creates/updates the global variables "dnapath.mart_obj" and 
#' "dnapath.mart_species" used by the 
#' \code{\link{symbol_to_entrez}} and \code{\link{entrez_to_symbol}} methods.
#' @note Requires internet connection. If unavailable, the default biomart
#' and default species contained in the package is used. 
#' @keywords internal
#' @export
init_mart <- function(species) {
  # listMarts(); listAttributes(mart) # Useful functions to obtain more info.
  # listDatasets(useMart('ensembl'));
  
  species <- format_species_name(species)
  
  # Obtain the biomaRt for this species.
  mart_temp <- tryCatch(
    biomaRt::useMart(biomart = "ensembl", 
                     dataset = paste0(species, "_gene_ensembl")),
    error = function(e) {
      if(!curl::has_internet()) {
        message("Error: Internet connection must be available.\n")
      } else if(!grepl("curl", e)) {
        message("Error: ", paste0('"', species, '"'), " is not an available species.\n",
                "Example species include: `H. sapiens`, `M. Musculus`, `C. elegans`, etc.\n",
                "Use the following command to see which species are available:\n",
                "> biomaRt::listDatasets(useMart('ensembl'))$dataset\n",
                "Note: biomaRt uses the format `hsapiens` for `Homo sapiens`.\n")
      }
      return(NULL)
    })
  
  if(is.null(mart_temp)) {
    message("Warning: Using default species `H. sapiens`.\n")
    mart_temp <- dnapath::mart
    assign("dnapath.mart_obj", mart_temp, inherits = TRUE)
    assign("dnapath.mart_species", "default", inherits = TRUE)
  } else {
    # Save the biomaRt in global environment.
    assign("dnapath.mart_obj", mart_temp, inherits = TRUE)
    assign("dnapath.mart_species", species, inherits = TRUE)
  }
  
  return(list(mart_obj = mart_temp, 
              mart_species = species))
}

#' Format sepcies name input.
#' 
#' Internal function used to format species names to be used with biomart.
#' @param species The species to obtain a biomart dataset for. For example:
#' "Homo sapiens", "m musculus", "C. elegans", or "S cerevisiae", and
#' "Human" and "mouse" will be converted to the format required in biomart.
#' @return A string containing the formatted species name.
#' @keywords internal
#' @export
format_species_name <- function(species) {
  # Convert to lowercase. If a vector of names is provided, use only the first.
  species <- tolower(species[1])
  
  if(species == "human") {
    species <- "h sapiens"
  } else if(species == "mouse") {
    species <- "m musculus"
  }
  
  # If species contains a space, like "homo sapiens", "h. sapiens", or "h sapiens", 
  # convert to format "hsapiens".
  if(grepl(" ", species)) {
    first <- substr(species, 1, 1)
    species <- paste0(first, gsub("(.* )", "", species))
  }
  if(grepl("_gene_ensembl", species)) {
    species <- gsub("_gene_ensembl", "", species)
  }
  
  return(species)
}