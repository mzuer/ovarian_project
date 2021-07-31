#
import datetime
start_time = str(datetime.datetime.now().time())
print('> START: tcga_met_vae.py \t' + start_time)

# python tcga_met_vae.py

from tensorflow.keras import backend as K
from tensorflow.keras import optimizers
from tensorflow.keras.layers import BatchNormalization as BN, Concatenate, Dense, Input, Lambda,Dropout
from tensorflow.keras.models import Model
from tensorflow.keras.losses import mean_squared_error,binary_crossentropy
from keras.utils.vis_utils import plot_model

import tensorflow as tf
import pandas as pd
import numpy as np
import os
import sys
import pickle
import seaborn as sns
import argparse
import matplotlib.pyplot as plt
import pandas as pd

import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from scipy.stats import spearmanr

wd = os.path.join('/home','marie','Documents','FREITAS_LAB','ovarian_project','vae')
os.chdir(wd)

pkg_dir = os.path.join('/home','marie','Documents','FREITAS_LAB','VAE_tutos','CancerAI-IntegrativeVAEs')
#module_path = r'C:\Users\d07321ow\Google Drive\SAFE_AI\CCE_DART\code\IntegrativeVAEs\code'
module_path = os.path.join(pkg_dir, 'code')

if module_path not in sys.path:
    sys.path.append(module_path)

from models.common import sse, bce, mmd, sampling, kl_regu
from misc.dataset import Dataset, DatasetWhole
from misc.helpers import normalizeRNA,save_embedding

from my_vae import *

outfolder = os.path.join('TCGA_MET_VAE')
os.makedirs(outfolder, exist_ok=True)

parser = argparse.ArgumentParser()
args = parser.parse_args()

###### PARAMETERS TO CHANGE BY THE USER #############

latent_dims = 64
args.ls = latent_dims # latent dimension size
args.ds = 256 # The intermediate dense layers size
args.distance = 'mmd'
args.beta = 1
args.act = 'elu'
args.epochs= 150
args.bs= 128  # Batch size
args.dropout = 0.2
args.save_model = True

outsuffix = "_" + str(args.epochs) + "epochs_"  + str(args.bs) + "bs"
args.out_model_file = os.path.join(outfolder, "metVAE" + outsuffix + ".h5")
args.out_history_file = os.path.join(outfolder, 'history'+ outsuffix +'.sav')

###### load methylation data

df = pd.read_csv(os.path.join('..','tcga_data','DOWNLOAD_TCGADATA_TCGABIOLINKS','TCGAbiolinks_OV_DNAmet27_hg38_metDT.txt'),
                 sep="\t")

# if there is na -> the model will not work  (visible as nan loss)
# i dont know what to do with the nan
df = df.fillna(0) 
#df.replace(np.nan,0)

keepFeatures = 1000
nskip = 4
ncolfull = df.shape[1]
meth_data = df.iloc[:,nskip:ncolfull].copy()
assert meth_data.shape[1] + nskip == ncolfull
# beta values -> no need to scale

meth_data = meth_data.iloc[:,0:keepFeatures]
assert meth_data.shape[1] == keepFeatures
meth_data_df = meth_data.copy()
meth_data = meth_data.values

print("... # samples = " + str(meth_data.shape[0]))
print("... # probes = " + str(meth_data.shape[1]))

args.input_size = meth_data.shape[1] # number of features (genes/probes)

assert args.input_size == keepFeatures

########## run the model

cncvae = CNCVAE(args)
cncvae.build_model()

cncvae.train(meth_data, meth_data)
emb_train = cncvae.predict(meth_data) # this it the latent space representation !

filename = os.path.join(outfolder,'cncvae_vae' + outsuffix + '.sav')
cncvae.vae.save(filename)
print("... written: " + filename )

filename = os.path.join(outfolder,'cncvae_decoder' + outsuffix + '.sav')
cncvae.decoder.save(filename)
print("... written: " + filename )

filename = os.path.join(outfolder,'cncvae_encoder' + outsuffix + '.sav')
cncvae.encoder.save(filename)
print("... written: " + filename )

filename = os.path.join(outfolder,'emb_train' + outsuffix + '.sav')
pickle.dump(emb_train, open(filename, 'wb'))
print("... written: " + filename )

outfile = os.path.join(outfolder, "methData_ls64_hs256_mmd_beta1_scaled" + outsuffix + ".csv")
np.savetxt(outfile, emb_train, delimiter = ',')
print("... written: " + outfile )

#################### plot the training performance

file = open(args.out_history_file, 'rb')
history = pickle.load(file)

loss_train = history['loss']
loss_val = history['val_loss']
epochs = range(1,args.epochs+1)
plt.plot(epochs, loss_train, 'g', label='Training loss')
plt.plot(epochs, loss_val, 'b', label='validation loss')
plt.title('Training and Validation loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.show()

filename = os.path.join(outfolder, 'train_valid_loss'+ outsuffix +'.png')
plt.savefig(filename, dpi=300) 
print('... written: ' + filename)


sys.exit(0)

#################### plot representations

latent_repr = emb_train
assert latent_repr.shape[0] == meth_data.shape[0]
assert latent_repr.shape[1] == latent_dims

df_plot = df.copy()
df_plot['figo_stages'] = df_plot['figo_stages'].replace(0,"-")

df_plot['figo_stages_short'] =  df_plot['figo_stages']
df_plot['figo_stages_short'] = df_plot['figo_stages_short'].str.replace('A',"")
df_plot['figo_stages_short'] = df_plot['figo_stages_short'].str.replace('B',"")
df_plot['figo_stages_short'] = df_plot['figo_stages_short'].str.replace('C',"")

all_label_cols = ['figo_stages_short', 'sample_types']

for labcol in all_label_cols:
    
    # PCA
    outfile = os.path.join(outfolder, "latent_vs_all_repr_" + labcol + "_PCA" + outsuffix + ".png")
    pca = PCA(n_components=2)
    pca.fit(latent_repr)
    latent_repr_pca = pca.transform(latent_repr)
    assert latent_repr_pca.shape == (meth_data.shape[0], 2)
    pca2 = PCA(n_components=2)
    pca2.fit(meth_data)
    fullData_pca = pca2.transform(meth_data)
    assert fullData_pca.shape == (meth_data.shape[0], 2)
    plot_2plots(ld_toplot_dt=latent_repr_pca, 
                raw_toplot_dt=fullData_pca, 
                labels=df_plot[labcol], 
                ld_pca_evr = pca.explained_variance_ratio_.sum(),
                raw_pca_evr = pca2.explained_variance_ratio_.sum(),
                dr_type="PCA", file_name=None)
    
    # PLOT UMAP
    outfile = os.path.join(outfolder, "latent_vs_all_repr_" + labcol + "_UMAP" + outsuffix + ".png")
    mapper = umap.UMAP(n_neighbors=15, n_components=2).fit(latent_repr)
    latent_repr_umap = mapper.transform(latent_repr)
    assert latent_repr_umap.shape == (meth_data.shape[0], 2)
    mapper2 = umap.UMAP(n_neighbors=15, n_components=2).fit(meth_data)
    fullData_umap = mapper2.transform(meth_data)
    assert fullData_umap.shape == (meth_data.shape[0], 2)
    plot_2plots(ld_toplot_dt=latent_repr_umap, 
                    raw_toplot_dt=fullData_umap, 
                    labels=df_plot[labcol], 
                    dr_type="UMAP", file_name=None)
    
    
    # PLOT TSNE
    outfile = os.path.join(outfolder, "latent_vs_all_repr_" + labcol + "_TSNE" + outsuffix + ".png")
    latent_repr_tsne = TSNE(n_components=2, perplexity=30 ).fit_transform(latent_repr)
    assert latent_repr_tsne.shape == (meth_data.shape[0], 2)
    fullData_tsne = TSNE(n_components=2, perplexity=30 ).fit_transform(meth_data)
    assert fullData_tsne.shape == (meth_data.shape[0], 2)
    plot_2plots(ld_toplot_dt=latent_repr_tsne, 
                    raw_toplot_dt=fullData_tsne, 
                    labels=df_plot[labcol], 
                    dr_type="t-SNE", file_name=None)
    






#####################################################################################################################


correlations_all=[]
p_values_all=[]
for gene_i in range(meth_data.shape[1]):
    correlations=[]
    p_values=[]
    for latent_dim_i in range(latent_dims):
        corr_, p_value = spearmanr(meth_data[:,gene_i], latent_repr[:,latent_dim_i])
        correlations.append(corr_)
        p_values.append(p_value)
    correlations_all.append(correlations)
    p_values_all.append(p_values)

correlations_all = np.array(correlations_all)
correlations_all_df = pd.DataFrame(correlations_all.T, columns = meth_data_df.columns)
p_values_all = np.array(p_values_all)
p_values_all_df  = pd.DataFrame(p_values_all.T, columns = meth_data_df.columns)



######### clustering
labels = df_plot['figo_stages_short'].values

lut = dict(zip(set(labels), sns.hls_palette(len(set(labels)))))
col_colors = pd.DataFrame(labels)[0].map(lut)

sns.clustermap(correlations_all_df, col_colors=col_colors)
out_file_name = os.path.join(outfolder, 'correlations_clustermap.png')
plt.savefig(out_file_name, dpi=300) 
print('... written: ' + out_file_name)

sns.clustermap(p_values_all_df)
out_file_name = os.path.join(outfolder, 'pvalues_clustermap.png')
plt.savefig(out_file_name, dpi=300) 
print('... written: ' + out_file_name)




all_corrs = correlations_all_df.values
all_corrs = all_corrs.flatten()

#plt.hist(all_corrs, color = 'blue', edgecolor = 'black',
#         bins = int(180/5))
# seaborn histogram
#sns.distplot(all_corrs, hist=True, kde=False, 
#             bins=int(180/5), color = 'blue',
#             hist_kws={'edgecolor':'black'})
# Density Plot and Histogram of all arrival delays
sns.distplot(all_corrs, hist=True, kde=True, 
             bins=int(180/5), color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4})
# Add labels
plt.title('Histogram of corr values')
plt.xlabel('Corr. coeff.')
plt.ylabel('Probes')

out_file_name = os.path.join(outfolder, 'correlations_distplot.png')
plt.savefig(out_file_name, dpi=300) 
print('... written: ' + out_file_name)



########## most correlated probes

sorted_min = correlations_all_df.apply(min, axis=0).sort_values()
sorted_max = correlations_all_df.apply(max, axis=0).sort_values(ascending=False)


probe_dt = pd.read_csv(os.path.join('..','tcga_data','DOWNLOAD_TCGADATA_TCGABIOLINKS','TCGAbiolinks_OV_DNAmet27_hg38_probesDT.txt'),
                 sep="\t")



probe2genes = dict(zip(probe_dt['Composite.Element.REF'],probe_dt['Gene_Symbol']))

sorted_min.index.map(probe2genes)
# Index(['.', 'CYP2W1;CYP2W1', 'EFNB1', 'PDILT;PDILT', 'TAAR5', 'SLN;SLN;SLN',
#        'INS;INS;INS;INS;INS;INS-IGF2;INS-IGF2', 'CEACAM8;PSG3;PSG3;PSG3;PSG3',
#        'MNDA;MNDA', 'CD1B;CD1B',

sorted_max.index.map(probe2genes)       
# Index(['KRTAP19-7;KRTAP19-7', 'LCE2D', 'MBL2', 'SPRR3;SPRR3;SPRR3;SPRR3',
#        'KRTAP19-3;KRTAP19-3', 'CLDN17', 'LCE4A;LCE4A',
#        'MS4A2;MS4A2;MS4A2;MS4A2', '.', 'LYZL4;LYZL4',
#        .


# =============================================================================
# for latent_dim_i in range(latent_dims):
#     fig, ax = plt.subplots(figsize=(15,6))
#     corrs = correlations_all_df.iloc[latent_dim_i,:]
#     corrs.sort_values(ascending=False)[:30].plot.bar(ax=ax)
# 
# out_file_name = os.path.join(outfolder, 'correlations_barplot.png')
# plt.savefig(out_file_name, dpi=300) 
# print('... written: ' + out_file_name)
# 
# for latent_dim_i in range(latent_dims):
#     fig, ax = plt.subplots(figsize=(15,6))
#     p_values = p_values_all_df.iloc[latent_dim_i,:]
#     p_values.sort_values(ascending=True)[:30].plot.bar(ax=ax)
# out_file_name = os.path.join(outfolder, 'pvalues_barplot.png')
# plt.savefig(out_file_name, dpi=300) 
# print('... written: ' + out_file_name)
# 
# =============================================================================
    
print('***** DONE\n' + start_time + " - " +  str(datetime.datetime.now().time()))