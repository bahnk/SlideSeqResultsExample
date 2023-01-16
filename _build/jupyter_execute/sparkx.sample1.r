library(SPARK)
library(Seurat)
library(dplyr)
library(future)
library(readr)
library(tibble)

dir.create("output", showWarnings=F)

counts <- Seurat::Read10X("sample1", gene.column=2)

spatial <-
    readr::read_csv("sample1.csv") %>%
    tibble::column_to_rownames("Barcode")

spatial <- spatial[colnames(counts),]

mt_idx <- grep("mt-", rownames(counts), ignore.case=T)
if ( length(mt_idx) != 0 )
{
   counts <- counts[-mt_idx,]
}

spark <- 
   SPARK::sparkx(
      count_in=counts,
      locus_in=as.matrix(spatial),
      numCores=12,
      option="mixture",
      verbose=T
   )

spark$res_mtest %>%
   tibble::rownames_to_column("gene") %>%
   readr::write_csv("output/results.csv")

sessionInfo()
