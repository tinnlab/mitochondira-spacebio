library(jsonlite)

setwd("./gene-set-enrichment")

download_and_process_data <- function(accid) {
  message(accid)
  path <- file.path("./data", accid)
  if (dir.exists(path)) return()
  
  dir.create(path)
  id <- strsplit(accid, "-")[[1]][2]
  url <- paste0("https://genelab-data.ndc.nasa.gov/genelab/data/glds/files/", id)
  document <- fromJSON(txt = url)
  files <- document$studies[[1]]$study_files

  # Download and extract metadata
  download_metadata(files, path)

  # Download and save differential expression data
  download_de_data(files, path)
}

download_metadata <- function(files, path) {
  metadata <- files$file_name[grepl("ISA.zip", files$file_name) & grepl("metadata", files$file_name)]
  metadata_url <- files$remote_url[files$file_name == metadata]
  temp <- tempfile()
  download.file(paste0("https://genelab-data.ndc.nasa.gov", metadata_url), temp)
  
  files_list <- unzip(temp, list = TRUE)
  s_files <- grep("^s_", files_list$Name)
  unzip(temp, files = files_list$Name[s_files], exdir = path)
  unlink(temp)
}

download_de_data <- function(files, path) {
  de_data <- files$file_name[grepl("_differential_expression.*\\.csv", files$file_name)]
  de_urls <- files$remote_url[files$file_name %in% de_data]
  
  for (url in de_urls) {
    filename <- strsplit(url, "=")[[1]][3]
    filename <- sub("\\.csv$", "", filename)
    data <- read.csv(paste0("https://genelab-data.ndc.nasa.gov", url), stringsAsFactors = FALSE)
    saveRDS(as.data.frame(data), file.path(path, paste0(filename, ".rds")))
  }
}

# Main execution
for (accid in accessionIds) {
  download_and_process_data(accid)
}

# Check datasets with multiple DE analysis files
list_data <- list.dirs("./data")
multi_de <- Filter(function(dir) {
  length(list.files(dir, pattern = "\\.rds$")) > 1
}, list_data)

print(multi_de)