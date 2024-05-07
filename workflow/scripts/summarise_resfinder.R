library(magrittr)

logger::log_threshold(logger::INFO)


read_res <- function(res_file, columns){
  sample <- dirname(res_file) %>%
    basename
  
  res_raw <- readr::read_tsv(
    file = res_file,
    skip = 1,
    col_names = names(columns),
    col_types = paste(columns, collapse ="")
  )
  
  dplyr::mutate(
    res_raw,
    sample = sample,
    module = "ResFinder"
  )
}


wrangle_res <- function(res, id_thresh, cov_thresh){
  res_wrangled <- dplyr::mutate(
    res,
    phenotype = stringr::str_replace_all(string = phenotype, pattern = " ", replacement = "_"),
    phenotype = stringr::str_replace_all(string = phenotype, pattern = "-", replacement = "_"),
    partial = identity < id_thresh | coverage < cov_thresh,
    gene = dplyr::case_when(
      identity == 100 & coverage == 100 ~ res_gene,
      identity == 100 ~ paste(res_gene, paste0("(", coverage, "% COV)")),
      coverage == 100 ~ paste(res_gene, paste0("(", identity, "% ID)")),
      TRUE ~ paste(res_gene, paste0("(", identity, "% ID, ", coverage, "% COV)"))
    ),
    class = dplyr::case_when(
      stringr::str_detect(string = phenotype, pattern = "Warning") ~ "Others",
      partial ~ paste(phenotype, "Partial", sep = "__"),
      TRUE ~ phenotype
    )
  )
}


import_res <- function(res_file, id_thresh, cov_thresh, columns){
  res_import <- read_res(res_file, columns = columns)
  res_wrangled <- wrangle_res(res = res_import, id_thresh = id_thresh, cov_thresh = cov_thresh)
  
}


summarise_resfinder <- function(res){
  dplyr::group_by(res, sample, class) %>%
    dplyr::summarise(gene = paste(gene, collapse = ", "), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = class, values_from = gene) 
}


write_summary <- function(res, outfile){
  type_str <- stringr::str_extract(string = outfile, pattern = "\\.[A-Za-z]*$")
  
  logger::log_info("Writing ResFinder results to: ", outfile)
  if (tolower(type_str) == ".tsv"){
    readr::write_tsv(x = res, file = outfile)
  } else if (tolower(type_str) == ".csv"){
    readr::write_csv(x = res, file = outfile)
  } else {
    stop("Wrong file extension for output file: ", type_str, "\nMust be either .tsv or .csv")
  }
}

import_all_res <- function(res_files, id_thresh, cov_thresh, columns, outfile){
  res_import <- purrr::map_dfr(.x = res_files, .f = import_res, id_thresh = id_thresh, cov_thresh = cov_thresh, columns = columns)
  
  res_summarised <- dplyr::arrange(res_import, class, res_gene) %>%
    summarise_resfinder
  
  others_there <- "Ohters" %in% colnames(res_summarised)
  if (others_there)
    res_summarised <- dplyr::relocate(res_summarised, Others, .after = dplyr::last_col())

  
  write_summary(res = res_summarised, outfile)
    
}


save(snakemake, file = "snakemake.RData")


import_all_res(
  res_files = snakemake@input[["res_files"]],
  id_thresh = snakemake@params[["id_thresh"]],
  cov_thresh = snakemake@params[["cov_thresh"]],
  outfile = snakemake@output[["resfinder_results"]],
  columns = c(
    "res_gene" = "c",
    "identity" = "d",
    "reference_alignment" = "c",
    "coverage" = "d",
    "reference_position" = "c",
    "contig_id" = "c",
    "contig_position" = "c",
    "phenotype" = "c",
    "gene_accession" = "c"
  )
)

