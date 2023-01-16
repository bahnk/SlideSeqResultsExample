library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(plotly)
library(readr)
library(spacexr)
library(tibble)

dir.create("output", showWarnings=F)

ref_counts <- Seurat::Read10X("reference", gene.column=2)

cell_types <-
   readr::read_tsv(file.path("reference", "types.tsv.gz"), col_names=F) %>%
   dplyr::pull(X1) %>%
   setNames(colnames(ref_counts)) %>%
   as.factor()

ref_nUMIs <- colSums(ref_counts)

ref <- spacexr::Reference(round(ref_counts), cell_types, ref_nUMIs)
saveRDS(ref, "output/reference.rds")

counts <- Seurat::Read10X("sample1", gene.column=2)

nUMIs <- colSums(counts)

spatial <-
    readr::read_csv("sample1.csv") %>%
    tibble::column_to_rownames("Barcode")

spatial <- spatial[colnames(counts),]

spots <- spacexr::SpatialRNA(spatial, counts, nUMIs)

plan("multicore", workers=12)
rctd <- spacexr::create.RCTD(spots, ref, max_cores=12)

rctd <- spacexr::run.RCTD(rctd, doublet_mode="doublet")

saveRDS(rctd, "output/rctd.rds")

readr::write_csv(
   rctd@results$results_df %>% tibble::rownames_to_column("cell"),
   "output/rctd.csv"
)

if ( "first_type" %in% colnames(rctd@results$results_df) )
{
   df <-
      rctd@results$results_df %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::inner_join( spatial %>% tibble::rownames_to_column("cell") )
   
   fig <- plotly::plot_ly(data=df, x=~x, y=~y, color=~first_type)
   
   plotly::save_image(fig, "output/first_type.png")
   plotly::save_image(fig, "output/first_type.pdf")
   
   #fig

   g <- ggplot2::ggplot(df, aes(x=x, y=y, color=first_type)) + geom_point()
   g
}

if ( "second_type" %in% colnames(rctd@results$results_df) )
{
   df <-
      rctd@results$results_df %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::inner_join( spatial %>% tibble::rownames_to_column("cell") )
   
   fig <- plotly::plot_ly(data=df, x=~x, y=~y, color=~second_type)
   
   plotly::save_image(fig, "output/second_type.png")
   plotly::save_image(fig, "output/second_type.pdf")
   
   #fig

   g <- ggplot2::ggplot(df, aes(x=x, y=y, color=second_type)) + geom_point()
   g
}

sessionInfo()
