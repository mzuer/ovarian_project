
# python ov_tcga_norm_tum_VAE.py

### !!!!! il faudrait avoir des STICs plutôt que des norm !!!!

from tensorflow.keras import backend as K
from tensorflow.keras import optimizers
from tensorflow.keras.layers import BatchNormalization as BN, Concatenate, Dense, Input, Lambda,Dropout
from tensorflow.keras.models import Model

import tensorflow as tf
import numpy as np

import os
import sys

import pandas as pd

wd = os.path.join('/home','marie','Documents','FREITAS_LAB','ovarian_project','CancerAI_VAE')
os.chdir(wd)


module_path = os.path.join('/home','marie','Documents','FREITAS_LAB','VAE_tutos','CancerAI-IntegrativeVAEs', 'code')
if module_path not in sys.path:
    sys.path.append(module_path)


from models.common import sse, bce, mmd, sampling, kl_regu
from tensorflow.keras.losses import mean_squared_error,binary_crossentropy
import numpy as np

from misc.dataset import Dataset, DatasetWhole
from misc.helpers import normalizeRNA,save_embedding

outfolder = os.path.join('OV_TCGA_NORM_TUM_MOSTVAR_VAE')
os.makedirs(outfolder, exist_ok=True)

import vae_ow
from vae_ow import *


import argparse
parser = argparse.ArgumentParser()
args = parser.parse_args()

###### PARAMETERS TO CHANGE BY THE USER #############

latent_dims = 100 # was 64
args.ls = latent_dims # latent dimension size
args.ds = 256 # The intermediate dense layers size
args.distance = 'mmd'
args.beta = 1
        
args.act = 'elu'
args.epochs= 150
args.bs= 128  # Batch size
args.dropout = 0.2
args.save_model = False

outsuffix = "_" + str(latent_dims) + "LD_" + str(args.ds) + "DS_" + str(args.epochs) + "_"+ str(args.bs) + "_" + str(args.dropout) + "_" + str(args.beta) + "_" + str(args.distance) + "."

###### END #############

input_data = os.path.join('/home', 'marie' ,'Documents', 'FREITAS_LAB', 'VAE_tutos/NeuralDec/ovarian_data/log_filtNormCount_mostVar_dt.txt')

# training data
df=pd.read_csv(input_data)
del df['geneID_symb']
del df['geneID_ensembl']

args.input_size = df.shape[0]  ## !!! will need to try with most Var genes !!!

in_data = df.copy().values
in_data_scaled = (in_data - in_data.min(axis=1).reshape(-1,1))/ (in_data.max(axis=1)-in_data.min(axis=1)).reshape(-1,1)

in_data_scaled = in_data_scaled.T

cncvae = CNCVAE(args)
cncvae.build_model()

cncvae.train(in_data_scaled, in_data_scaled)
emb_train = cncvae.predict(in_data_scaled) # this it the latent space representation !



latent_repr = emb_train

outfile = os.path.join(outfolder, "latent_repr" + outsuffix +"csv")
np.savetxt(outfile, latent_repr, delimiter = ',')
print("... written: " + outfile )


#####################################################################################################################

in_data = in_data.T

from scipy.stats import spearmanr
correlations_all=[]
p_values_all=[]
for gene_i in range(in_data.shape[1]):
    correlations=[]
    p_values=[]
    for latent_dim_i in range(latent_dims):
        
        corr_, p_value = spearmanr(in_data[:,gene_i], latent_repr[:,latent_dim_i])
        
        correlations.append(corr_)
        p_values.append(p_value)
        
    correlations_all.append(correlations)
    p_values_all.append(p_values)

correlations_all = np.array(correlations_all)

p_values_all = np.array(p_values_all)


outfile = os.path.join(outfolder, "correlations_all" + outsuffix + "csv")
np.savetxt(outfile, correlations_all, delimiter = ',')
print("... written: " + outfile )

outfile = os.path.join(outfolder, "p_values_all" + outsuffix + "csv")
np.savetxt(outfile, p_values_all, delimiter = ',')
print("... written: " + outfile )

