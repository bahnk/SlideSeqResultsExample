library(BiocParallel)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(future)
library(ggplot2)
library(patchwork)
library(purrr)
library(sceasy)

dir.create("output", showWarnings=F)

counts <- Seurat::Read10X("sample2", gene.column=2)

spatial <-
   readr::read_csv("sample2.csv") %>%
   dplyr::filter( Barcode %in% colnames(counts) ) %>%
   dplyr::arrange(Barcode) %>%
   tibble::column_to_rownames("Barcode")

stopifnot( sum( ! spatial$Barcode == colnames(counts) ) == 0 )

args <- list("project"="sample2", "assay"="Spatial", "meta.data"=spatial)
obj <- do.call(SeuratObject::CreateSeuratObject, c(counts, args))

img <- new(
   "SlideSeq",
   coordinates=obj@meta.data[,c("x","y")],
   assay="Spatial",
   key="_images"
   )

obj@images <- list(image=img)

obj$log_nCount_Spatial <- log(obj$nCount_Spatial)

obj %>%
   Seurat::VlnPlot(., features="nCount_Spatial", pt.size=0, log=TRUE) +
   Seurat::NoLegend() -> vln_count_plot

ggplot2::ggsave("output/violin_count.png", vln_count_plot, width=9, height=6)
ggplot2::ggsave("output/violin_count.pdf", vln_count_plot, width=9, height=6)

vln_count_plot

obj %>%
   Seurat::SpatialFeaturePlot(., features="log_nCount_Spatial") +
   ggplot2::theme(legend.position="right") -> spatial_umi_plot

ggplot2::ggsave("output/spatial_umi.png", spatial_umi_plot, width=9, height=6)
ggplot2::ggsave("output/spatial_umi.pdf", spatial_umi_plot, width=9, height=6)

spatial_umi_plot

obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern="^mt-")

obj %>%
   Seurat::VlnPlot(., features="percent.mt", pt.size=0, log=TRUE) +
   Seurat::NoLegend() -> vln_mitoch_plot

ggplot2::ggsave("output/violin_mitoch.png", vln_mitoch_plot, width=9, height=6)
ggplot2::ggsave("output/violin_mitoch.pdf", vln_mitoch_plot, width=9, height=6)

vln_mitoch_plot

scatter_plot <-
   Seurat::FeatureScatter(
      obj,
      feature1="nCount_Spatial",
      feature2="percent.mt"
   )

ggplot2::ggsave("output/scatter_count_mitoch.png", scatter_plot, width=9, height=6)
ggplot2::ggsave("output/scatter_count_mitoch.pdf", scatter_plot, width=9, height=6)

scatter_plot

cat("Original\n")
print(obj)

cat("Filter mitochondria\n")
obj <- base::subset(obj, percent.mt < 30)
print(obj)

cat("Filter UMIs\n")
obj <- base::subset(obj, nCount_Spatial >= 10)
print(obj)

obj <- Seurat::SCTransform(obj, assay="Spatial", ncells=3000, verbose=T)

obj <- Seurat::RunPCA(obj)

obj <- Seurat::RunUMAP(obj, dims=1:10)

obj <- Seurat::FindNeighbors(obj, dims=1:10)

BiocParallel::register(MulticoreParam(12))
base::options(future.globals.maxSize = 5000 * 1024^2)
future::plan("multiprocess", workers=12)

obj <- Seurat::FindClusters(obj, resolution=0.1, verbose=T)

umap_plot <- Seurat::DimPlot(obj, reduction="umap", label=TRUE)

ggplot2::ggsave("output/umap.png", umap_plot, width=9, height=6)
ggplot2::ggsave("output/umap.pdf", umap_plot, width=9, height=6)

umap_plot

spatial_cluster_plot <- Seurat::SpatialDimPlot(obj, stroke=0)

ggplot2::ggsave("output/spatial.png", spatial_cluster_plot, width=9, height=6)
ggplot2::ggsave("output/spatial.pdf", spatial_cluster_plot, width=9, height=6)

spatial_cluster_plot

obj %>%
   slot(., "meta.data") %>%
   dplyr::pull(paste0("SCT_snn_res.", 0.1)) %>%
   levels() %>%
   as.numeric(.) %>%
   sort() -> clusters


ncol <- 3
1:ceiling( length(clusters) / ncol) %>%
   purrr::map(function(x) ncol * (x-1) + 1:ncol) %>%
   purrr::map(function(x) clusters[x]) %>%
   purrr::map(function(x) {

      plot <- Seurat::SpatialDimPlot(
         obj,
         cells.highlight=Seurat::CellsByIdentities(object=obj, idents=x),
         facet.highlight=T
      )

      basename <- paste0("output/individual_cluster_", paste(na.omit(x), collapse="-"))

      ggplot2::ggsave(
         filename=paste0(basename, ".png"),
         plot=plot,
         widt=9, height=6, dpi=300
         )

      ggplot2::ggsave(
         filename=paste0(basename, ".pdf"),
         plot=plot,
         widt=9, height=6, dpi=300
         )

      print(plot)

      list("clusters"=x, "plot"=plot)

   }) -> clusters_plots

markers <- Seurat::FindAllMarkers(obj, assay="SCT", only.pos=T)
write.csv(markers, "output/markers.csv", row.names=F)
markers

if ( nrow(markers) > 0 )
{
   markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::group_map(function(x, n) {
   
         df <-
            x %>%
            dplyr::arrange(dplyr::desc(avg_log2FC) ) %>%
            dplyr::top_n( max(6, nrow(x)) )
   
         genes <- df$gene
   
         for ( gene in genes ) {
   
            g <- Seurat::FeaturePlot(obj, features=gene, ncol=1)
            g <- g + patchwork::plot_annotation(title=paste0("Cluster ", n))
   
            ggplot2::ggsave(
               filename=sprintf("output/markers_cluster_%s_gene_%s.png", n, gene),
               plot=g,
               height=3, width=3.5, dpi=200
               )
   
            ggplot2::ggsave(
               filename=sprintf("output/markers_cluster_%s_gene_%s.pdf", n, gene),
               plot=g,
               height=3, width=3.5, dpi=200
               )
   
            print(g)
         }
      })
}

saveRDS(obj, "output/seurat_object.rds")

sceasy::convertFormat(
   obj,
   from="seurat",
   to="anndata",
   assay="Spatial",
   outFile="output/seurat_object.h5ad"
   )

sessionInfo()
