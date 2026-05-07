"""
editor: yangdong
date: 260506
"""

atac4cistopic_dir="/data/work/scenic/Os_rds2cistopic"
mallet_path="/data/work/scenic/Mallet-202108/bin/mallet"
n_cpu=16
mallet_mem='128G'
numTopics = 45
atac_key='sctype'

import os 
os.chdir('/data/work/scenic/cistopic2/')

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
# sys.path.append('/home/sussmanj/miniconda3/envs/scenicplus/lib/')
_stderr = sys.stderr
null = open(os.devnull,'wb')


import importlib
import sys
import os 
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from scenicplus.scenicplus_class import SCENICPLUS, create_SCENICPLUS_object
from scenicplus.loom import *
from scenicplus.preprocessing.filtering import *
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import evaluate_models
from pycisTopic.clust_vis import run_umap
from loomxpy.loomxpy import SCopeLoom
from pycisTopic.loom import *
import itertools
import anndata
import pickle
import dill
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
from pycisTopic.lda_models import run_cgs_models_mallet
import scipy.sparse as sp
import scipy.io as spio
import pyranges as pr
import numpy as  np
from pycistarget.utils import region_names_to_coordinates
from pycisTopic.lda_models import run_cgs_models
from scenicplus.wrappers.run_pycistarget import run_pycistarget
import pybiomart as pbm
from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
from scenicplus.eregulon_enrichment import score_eRegulons
import seaborn as sns
from scenicplus.plotting.dotplot import heatmap_dotplot
from scenicplus.networks import create_nx_tables, create_nx_graph, plot_networkx, export_to_cytoscape
from scenicplus.RSS import *
from scenicplus.plotting.correlation_plot import *
from pycisTopic.diff_features import find_highly_variable_features
from scenicplus.networks import create_nx_tables, create_nx_graph, plot_networkx, export_to_cytoscape
import seaborn as sns
import sklearn
from scenicplus.differentiation_potential import *


#Create cisTopic object from the counts matrix 
sparse_csr_matrix=spio.mmread(os.path.join(atac4cistopic_dir, 'atac_peaks_sparse.mtx')).tocsr()
cell_list = []
region_list = [] 

with open(os.path.join(atac4cistopic_dir, 'atac_cell_names.txt'), "r") as file:
    cell_list = [line.strip() for line in file]


with open(os.path.join(atac4cistopic_dir, 'atac_region_names.txt'), "r") as file:
    region_list = [line.strip() for line in file]


region_list = [s.replace('-', ':', 1) for s in region_list]

# path_to_blacklist='/mnt/isilon/tan_lab/sussmanj/Single_Cell_Tools/ScenicPlus/Genome_Files/hg38-blacklist.v2.bed'
path_to_blacklist=None

cistopic_obj = create_cistopic_object(fragment_matrix=sparse_csr_matrix, cell_names=cell_list, region_names=region_list,
                                      path_to_blacklist=path_to_blacklist)

#Add cell metadata
cell_data =  pd.read_csv(os.path.join(atac4cistopic_dir, 'atac_metadata.txt'), sep='\t')
cell_data.index = cell_list
cistopic_obj.add_cell_data(cell_data)

if not os.path.exists('Tmp_Files'):
    os.makedirs('Tmp_Files')

# !wget https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz
# !tar -xf Mallet-202108-bin.tar.gz

os.environ['MALLET_MEMORY'] = mallet_mem # 700G'
# models=run_cgs_models_mallet(cistopic_obj,
#                     mallet_path = mallet_path, 
#                     n_topics=[45],
#                     n_cpu=n_cpu, # 12
#                     n_iter=500,
#                     random_state=555,
#                     alpha=50,
#                     alpha_by_topic=True,
#                     eta=0.1,
#                     eta_by_topic=False,
#                     save_path = 'Tmp_Files',
#                     tmp_path = 'Tmp_Files')
models=run_cgs_models_mallet(cistopic_obj,
                    mallet_path = mallet_path, 
                    n_topics=[2,5,10,15,20,25,30,35,40,45,50],
                    n_cpu=n_cpu, # 12
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path = 'Tmp_Files',
                    tmp_path = 'Tmp_Files')

pickle.dump(cistopic_obj, open('ATAC_cistopic_obj.pkl', 'wb'))
pickle.dump(models, open('ATAC_Models_500_iter_LDA.pkl', 'wb'))

# #Load models 
# cistopic_obj = pickle.load(open('ATAC_cistopic_obj.pkl', 'rb'))
# models = pickle.load(open('ATAC_Models_500_iter_LDA.pkl', 'rb'))


model = evaluate_models(models,
                     select_model = numTopics,
                     return_model = True,
                     metrics = ['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics = True,
                     save = "model_metrics.pdf")

#Add model to cisTopic object and save again 
cistopic_obj.add_LDA_model(model)
cistopic_obj.cell_data[atac_key] = cistopic_obj.cell_data[atac_key].str.replace(' ', '_')
pickle.dump(cistopic_obj,
            open('ATAC_cistopic_obj_with_model.pkl', 'wb'))

#Run UMAP
run_umap(cistopic_obj, target = "cell", scale = True)

print(cistopic_obj.cell_data.columns.tolist())

plot_metadata(
    cistopic_obj,
    reduction_name = 'UMAP',
    variables = [atac_key],
    figsize = (10, 10), 
    save = "UMAP_CisTopics.pdf")

plot_topic(cistopic_obj, reduction_name = 'UMAP', num_columns = 5, 
           save = "FeaturePlot_CisTopics.pdf")


#Binarize topics, using the 'otsu' method
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

#Calculating DARs per cell type 
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable=atac_key, var_features=variable_regions, split_pattern = '-')

os.makedirs("region_sets", exist_ok = True)
os.makedirs(os.path.join("region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join("region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join("region_sets", "DARs_cell_type"), exist_ok = True)

for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join("region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )


for topic in region_bin_topics_top3k:
    region_names_to_coordinates(
        region_bin_topics_top3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join("region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )


for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join("region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )