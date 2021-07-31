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

from my_vae import CNCVAE

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

# PCA
pca = PCA(n_components=2)
pca.fit(latent_repr)
latent_repr_pca = pca.transform(latent_repr)
assert latent_repr_pca.shape == (meth_data.shape[0], 2)

fullData_pca = PCA(n_components=2)
fullData_pca.fit(meth_data)
fullData_pca = fullData_pca.transform(meth_data)
assert fullData_pca.shape == (meth_data.shape[0], 2)


outfile = os.path.join(outfolder, "latent_repr_pca" + outsuffix + ".png")

fig, axs = plt.subplots(1,2,figsize = (12,6))
palette = 'tab10'
g = sns.scatterplot(latent_repr_pca[:,0], latent_repr_pca[:,1],
                hue = list(df_plot['figo_stages_short']), 
                ax=axs[0],
                linewidth=0, s=15, alpha=0.7, palette = palette)
g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=1)
g.set(title="A")

g = sns.scatterplot(fullData_pca[:,0], fullData_pca[:,1],
                hue = list(df_plot['figo_stages_short']), 
                ax=axs[1],
                linewidth=0, s=15, alpha=0.7, palette = palette)
g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=len(set(labels)))


def plot_2plots(ld_toplot_dt, raw_toplot_dt, labels, dr_type, file_name=None):
    
    assert ld_toplot_dt.shape[0] == raw_toplot_dt.shape[0]
    assert ld_toplot_dt.shape[0] == len(labels)
    
    fig, axs = plt.subplots(1,2,figsize = (12,6))
    palette = 'tab10'
    g = sns.scatterplot(ld_toplot_dt[:,0], ld_toplot_dt[:,1],
                    hue = list(labels),
                    ax=axs[0],
                    linewidth=0, s=15, alpha=0.7, palette = palette)
    g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=len(set(labels)))
    g.set(title='latent dims')
    
    g = sns.scatterplot(raw_toplot_dt[:,0], raw_toplot_dt[:,1],
                    hue = list(labels),
                    ax=axs[1],
                    linewidth=0, s=15, alpha=0.7, palette = palette)
    #g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=1)
    g.set(title='full data')
    g.legend_.remove()
    
    for ax in axs:
        ax.set_xlabel('{} 1'.format(dr_type))
        ax.set_ylabel('{} 2'.format(dr_type))

    
    fig.suptitle('{}'.format(dr_type), x=0.5, y=0.99)


    if file_name:
        plt.savefig(file_name, dpi=300) 


plot_2plots(ld_toplot_dt=latent_repr_pca, 
            raw_toplot_dt=fullData_pca, 
            labels=df_plot['figo_stages_short'], 
            dr_type="PCA", file_name=None)



plot_3plots(data_to_plot=latent_repr_pca, data_with_labels=df, type_='PCA', pca=pca, file_name=outfile)




###########################################################################################
    
def plot_3plots(data_to_plot, data_with_labels,file_name='', type_ = 'PCA', pca=None):
    
    fig, axs = plt.subplots(1,3,figsize = (15,6))
    palette = 'tab10'
    ### ! depending on the version of matplotlib -> should pass a list to hue !!!!
    g = sns.scatterplot(data_to_plot[:,0], data_to_plot[:,1],
                        #hue = data_with_labels['ER_Expr'], 
                        hue = list(data_with_labels['ER_Expr']), 
                        ax=axs[0],linewidth=0, s=15, alpha=0.7, palette = palette)
    g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=1)

    g = sns.scatterplot(data_to_plot[:,0], data_to_plot[:,1],
                        # hue = data_with_labels['Pam50Subtype'], 
                        hue = list(data_with_labels['Pam50Subtype']), 
                        ax=axs[1],linewidth=0, s=15, alpha=0.7, palette = palette)
    g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3)
    g = sns.scatterplot(data_to_plot[:,0], data_to_plot[:,1],
                        # hue = data_with_labels['iC10'], 
                        hue = list(data_with_labels['iC10']), 
                        ax=axs[2],linewidth=0, s=15, alpha=0.7, palette = palette)
    g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=5)
    #plt.
    # ax[0].plot(latent_repr_pca[:,0], latent_repr_pca[:,1], '.')
    # ax[0].plot(latent_repr_pca[:,0], latent_repr_pca[:,1], '.')
    
    for ax in axs:
        ax.set_xlabel('{} 1'.format(type_))
        ax.set_ylabel('{} 2'.format(type_))
    
    if type_ =='PCA':
        fig.suptitle('{}\nPCA - explained variance ratio: {:.2f}'.format(file_name,pca.explained_variance_ratio_.sum()), x=0.5, y=0.99)
    else:
        fig.suptitle('{}\n{}'.format(file_name,type_), x=0.5, y=0.99)
        
    plt.tight_layout()
    
    if file_name != '':
        plot_file_name = str.replace(file_name, '\\','_').split('.')[0]
        #out_file_name = os.path.join(outfolder,'downstream_results/{}_{}.png'.format(plot_file_name, type_)) # r -> treated as raw string
        out_file_name = os.path.join('{}_{}.png'.format(plot_file_name, type_)) # r -> treated as raw string
        plt.savefig(out_file_name, dpi=300) 
        print('... written: ' + out_file_name)
    return
    



# PLOT PCA
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
pca.fit(latent_repr)
latent_repr_pca = pca.transform(latent_repr)
#plot_3plots(data_to_plot=latent_repr_pca, data_with_labels=df, type_='PCA', pca=pca)
outfile = os.path.join(outfolder, "latent_repr_pca")
plot_3plots(data_to_plot=latent_repr_pca, data_with_labels=df, type_='PCA', pca=pca, file_name=outfile)


# PLOT UMAP
data_to_umap = latent_repr

mapper = umap.UMAP(n_neighbors=15, n_components=2).fit(data_to_umap)
latent_repr_umap = mapper.transform(data_to_umap)
plot_3plots(latent_repr_umap, df, type_='UMAP')
outfile = os.path.join(outfolder, "latent_repr_umap")
plot_3plots(data_to_plot=latent_repr_umap, data_with_labels=df, type_='UMAP', file_name=outfile)

# PLOT TSNE

latent_repr_tsne = TSNE(n_components=2, perplexity=30 ).fit_transform(latent_repr)
plot_3plots(latent_repr_tsne, df, type_='tSNE')
outfile = os.path.join(outfolder, "latent_repr_tsne")
plot_3plots(data_to_plot=latent_repr_tsne, data_with_labels=df, type_='tSNE', file_name=outfile)
    

# PLOT UMAP for RAW MRNA
mapper = umap.UMAP(n_neighbors=15, n_components=2).fit(mrna_data)
latent_repr_umap = mapper.transform(mrna_data)
outfile = os.path.join(outfolder, "latent_repr_umap_raw")
plot_3plots(latent_repr_umap, df, type_='UMAP')

#####################################################################################################################

from scipy.stats import spearmanr
correlations_all=[]
p_values_all=[]
for gene_i in range(mrna_data.shape[1]):
    correlations=[]
    p_values=[]
    for latent_dim_i in range(latent_dims):
        corr_, p_value = spearmanr(mrna_data[:,gene_i], latent_repr[:,latent_dim_i])
        correlations.append(corr_)
        p_values.append(p_value)
    correlations_all.append(correlations)
    p_values_all.append(p_values)

correlations_all = np.array(correlations_all)
correlations_all_df = pd.DataFrame(correlations_all.T, columns = df.iloc[:,34:1034].columns)
p_values_all = np.array(p_values_all)
p_values_all_df  = pd.DataFrame(p_values_all.T, columns = df.iloc[:,34:1034].columns)

labels = df['Pam50Subtype'].values

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


for latent_dim_i in range(latent_dims):
    fig, ax = plt.subplots(figsize=(15,6))
    corrs = correlations_all_df.iloc[latent_dim_i,:]
    corrs.sort_values(ascending=False)[:30].plot.bar(ax=ax)

out_file_name = os.path.join(outfolder, 'correlations_barplot.png')
plt.savefig(out_file_name, dpi=300) 
print('... written: ' + out_file_name)

for latent_dim_i in range(latent_dims):
    fig, ax = plt.subplots(figsize=(15,6))
    p_values = p_values_all_df.iloc[latent_dim_i,:]
    p_values.sort_values(ascending=True)[:30].plot.bar(ax=ax)
out_file_name = os.path.join(outfolder, 'pvalues_barplot.png')
plt.savefig(out_file_name, dpi=300) 
print('... written: ' + out_file_name)

    
print('***** DONE\n' + start_time + " - " +  str(datetime.datetime.now().time()))