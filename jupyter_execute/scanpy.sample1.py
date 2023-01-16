#!/usr/bin/env python
# coding: utf-8

# # sample1

# This is a very a very basic analysis of this sample with scanpy.

# Here are the libraries we need.

# In[1]:


from IPython.core.display import display, Image
from numpy.random import choice
from os import makedirs
from os.path import join
from plotly.subplots import make_subplots

import ipywidgets as widgets
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
import seaborn as sns
import session_info


# Some useful functions.

# In[2]:


##################
def in_jupyter():#
##################
   try:
      __IPYTHON__
   except NameError:
      return False
   else:
      return "IPKernelApp" in get_ipython().config.keys()
   ############################################################################


# We create the output directory for this noteboook.
# Every outputs will save there.

# In[3]:


try:
   makedirs("output")
except OSError as e:
   pass


# Configuration of scanpy.

# In[4]:


sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3


# We load the digital expression matrix and the spatial information and create an AnnData object.

# In[5]:


spatial = pd.read_csv("sample1.csv").set_index("Barcode")

adata = sc.read_10x_mtx("sample1", var_names="gene_symbols")
adata.obs = spatial.loc[ adata.obs.index ]
adata = adata[ : , ~ np.all(adata.X.toarray() == 0, axis=0) ]

# removes genes without a name (in Drosophila melanogaster for example)
adata = adata[ : , ~ adata.var.index.isna() ]


# ## Quality control
# 
# Scanpy will compute basic QC metrics for us.

# In[6]:


adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


# ### Gene counts
# 
# Genes detected for each bead.

# In[7]:


plots = [
   ("total_counts", "Total counts"),
   ("n_genes_by_counts", "Detected genes")
]

tab_content = []
tab_titles = []

n = min(adata.obs.shape[0], 3000)

for column, title in plots:

   fig = px.violin(adata.obs.sample(n), y=column, box=True, points="all")
   fig.write_image( join("output", f"qc_beads.{column}.pdf") )
   fig.write_image( join("output", f"qc_beads.{column}.png") )
   tab_content.append( go.FigureWidget(fig) )
   tab_titles.append(title)

   fig = px.scatter(adata.obs, x="x", y="y", color=column)
   fig.write_image( join("output", f"qc_beads.{column}.spatial.pdf") )
   fig.write_image( join("output", f"qc_beads.{column}.spatial.png") )
   tab_content.append( go.FigureWidget(fig) )
   tab_titles.append(f"{title} (xy)")

fig = px.scatter(adata.obs.sample(n), x="total_counts", y="n_genes_by_counts")
fig.write_image( join("output", "total_counts_VS_n_genes_by_counts.png") )
fig.write_image( join("output", "total_counts_VS_n_genes_by_counts.pdf") )

tab_content += [go.FigureWidget(fig)]
tab_titles += ["Counts VS Genes"]

# display plots
if in_jupyter():
   tab = widgets.Tab(tab_content)
   for i, title in enumerate(tab_titles):
      tab.set_title(i, title)
   display(tab)


# ### Genes
# 
# The dropouts.

# In[8]:


adata.var["detected"] = np.where(adata.var.n_cells_by_counts > 0,
   "Detected", "Undetected")

df = adata   .var   .detected   .value_counts()   .to_frame()   .reset_index()   .rename(columns={"index": "Gene", "detected": "Number"})

fig = px.bar(df, x="Gene", y="Number", text_auto=True)

tab_content = [go.FigureWidget(fig)]
tab_titles = ["Detected genes"]

plots = [
   ("total_counts", "Total counts"),
   ("mean_counts", "Mean counts"),
   ("n_cells_by_counts", "Beads"),
   ("pct_dropout_by_counts", "Dropout")
]

n = min(adata.obs.shape[0], 3000)

for column, title in plots:
   fig = px.violin(adata.var.sample(n), y=column, box=True, points="all")
   fig.write_image( join("output", f"qc_genes.{column}.pdf") )
   fig.write_image( join("output", f"qc_genes.{column}.png") )
   tab_content.append( go.FigureWidget(fig) )
   tab_titles.append(title)

# display plots
if in_jupyter():
   tab = widgets.Tab(tab_content)
   for i, title in enumerate(tab_titles):
      tab.set_title(i, title)
   display(tab)


# ### Mitochondrial activity
# 
# Percentage of mitochondrial genes detected for each bead.

# In[9]:


#("pct_counts_mt", "Percent Mitoch")

n = min(adata.obs.shape[0], 3000)

fig1 = px.violin(adata.obs.sample(n), y="pct_counts_mt", box=True, points="all")
fig1.write_image( join("output", "pct_counts_mt.pdf") )
fig1.write_image( join("output", "pct_counts_mt.png") )

fig2 = px.scatter(adata.obs.sample(n), x="total_counts", y="pct_counts_mt")
fig2.write_image( join("output", "total_counts_VS_pct_counts_mt.png") )
fig2.write_image( join("output", "total_counts_VS_pct_counts_mt.pdf") )

# display plots
if in_jupyter():
   tab = widgets.Tab([go.FigureWidget(fig1), go.FigureWidget(fig2)])
   tab.set_title(0, "Percent Mitoch")
   tab.set_title(1, "Counts VS Mitoch")
   display(tab)


# ## Filtering
# 
# We filter anormal beads.

# In[10]:


print("Before filtering:")
print(adata)

adata = adata[ adata.obs.n_genes_by_counts > 5 , :]
adata = adata[ adata.obs.n_genes_by_counts < np.infty , :]
adata = adata[ adata.obs.pct_counts_mt < 30 , :]

sc.pp.filter_genes(adata, min_cells=10)

print("\n\n")
print("After filtering:")
print(adata)


# ## Data transformation
# 
# Normalization by sequencing depth to remove technical variability.
# And nonlinear transformation to stabilize the variance across genes with different expression levels.

# In[11]:


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata)


# ## Principal components analysis
# 
# We compute the PCA.

# In[12]:


sc.pp.pca(adata, use_highly_variable=False)


# We plot the variance explained by the different principal components.

# In[13]:


var_ratio = adata.uns["pca"]["variance_ratio"]

fig = px.scatter(x=range(1, len(var_ratio)+1), y=var_ratio)
fig.write_image( join("output", "pca_variance_ratio.png") )
fig.write_image( join("output", "pca_variance_ratio.pdf") )

fig


# We plot the 5 first principal components.

# In[14]:


df = pd.DataFrame(
   np.hstack([adata.obs.index[:,None], adata.obsm["X_pca"][:,:5]]),
   columns=["Bead"] + [f"PC{i+1}" for i in range(5)]
   )\
   .melt(id_vars=["Bead", "PC1"], var_name="PC", value_name="Value")

fig = px.scatter(df, x="PC1", y="Value", facet_col="PC", facet_col_wrap=2)
fig.write_image( join("output", "pca.png") )
fig.write_image( join("output", "pca.pdf") )

fig


# ## Clustering & UMAP
# 
# We cluster the beads based on their expression.

# In[15]:


sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10, use_rep="X_pca")
sc.tl.leiden(adata, key_added="clusters", resolution=0.2)


# We compute the UMAP.

# In[16]:


sc.tl.umap(adata)


# We plot the clusters on top of the UMAP.

# In[17]:


args = {
      "x": np.apply_along_axis(lambda x: x[0], arr=adata.obsm["X_umap"], axis=1),
      "y": np.apply_along_axis(lambda x: x[1], arr=adata.obsm["X_umap"], axis=1),
      "color": adata.obs.clusters
}
fig = px.scatter(**args)
fig.write_image( join("output", "clusters_umap.png") )
fig


# We plot the clusters on top of the spatial coordinates.

# In[18]:


fig = px.scatter(adata.obs, x="x", y="y", color="clusters")
fig.write_image( join("output", "clusters_spatial.png") )
fig


# ## Save

# In[19]:


adata.write_h5ad( join("output", "anndata.h5ad") )


# ## Session info

# In[20]:


session_info.show()

