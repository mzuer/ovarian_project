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

class CNCVAE:
    def __init__(self, args):
        self.args = args
        self.vae = None
        self.encoder = None
        
        self.decoder=None

    def build_model(self):
        np.random.seed(42)
        # tf.random.set_random_seed(42)   # moved to (not able to install tensorflow 1.11)
        tf.random.set_seed(42)
        # Build the encoder network
        # ------------ Input -----------------
        inputs = Input(shape=(self.args.input_size,), name='concat_input')
        #inputs = [concat_inputs]

        # ------------ Encoding Layer -----------------
        x = Dense(self.args.ds, activation=self.args.act)(inputs)
        x = BN()(x)      

        # ------------ Embedding Layer --------------
        z_mean = Dense(self.args.ls, name='z_mean')(x)
        z_log_sigma = Dense(self.args.ls, name='z_log_sigma', kernel_initializer='zeros')(x)

        z = Lambda(sampling, output_shape=(self.args.ls,), name='z')([z_mean, z_log_sigma])

        self.encoder = Model(inputs, [z_mean, z_log_sigma, z], name='encoder')
        self.encoder.summary()

        # Build the decoder network
        # ------------ Dense out -----------------
        latent_inputs = Input(shape=(self.args.ls,), name='z_sampling')
        x = latent_inputs
        x = Dense(self.args.ds, activation=self.args.act)(x)
        x = BN()(x)
        
        x=Dropout(self.args.dropout)(x)

        # ------------ Out -----------------------
        
        #if self.args.integration == 'Clin+CNA':
        #    concat_out = Dense(self.args.input_size,activation='sigmoid')(x)
        #else:
        concat_out = Dense(self.args.input_size)(x)
        
        decoder = Model(latent_inputs, concat_out, name='decoder')
        
        decoder.summary()
        
        self.decoder=decoder

        outputs = decoder(self.encoder(inputs)[2])
        self.vae = Model(inputs, outputs, name='vae_mlp')
        

        # Define the loss
        if self.args.distance == "mmd":
            true_samples = K.random_normal(K.stack([self.args.bs, self.args.ls]))
            distance = mmd(true_samples, z)
        if self.args.distance == "kl":
            distance = kl_regu(z_mean,z_log_sigma)


        #if self.args.integration == 'Clin+CNA':
        #    reconstruction_loss = binary_crossentropy(inputs, outputs)
        #else:
        reconstruction_loss = mean_squared_error(inputs, outputs)
        vae_loss = K.mean(reconstruction_loss + self.args.beta * distance)
        self.vae.add_loss(vae_loss)

        adam = optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.001, amsgrad=False)
        self.vae.compile(optimizer=adam)
        self.vae.summary()
        
            
    def train(self, s_train, s_test):
        train = s_train#np.concatenate((s1_train,s2_train), axis=-1)
        test = s_test#np.concatenate((s1_test,s2_test), axis=-1)
        history = self.vae.fit(train, epochs=self.args.epochs, batch_size=self.args.bs, shuffle=True, validation_data=(test, None))
        if self.args.save_model:
            #self.vae.save_weights('./models/vae_cncvae.h5')
            self.vae.save_weights(self.args.out_model_file)
            pickle.dump(history.history, open(self.args.out_history_file, 'wb'))


    def predict(self, s_data):

        return self.encoder.predict(s_data, batch_size=self.args.bs)[0]



def plot_2plots(ld_toplot_dt, raw_toplot_dt, labels, dr_type, 
                ld_pca_evr=None, raw_pca_evr=None, file_name=None):
    
    assert ld_toplot_dt.shape[0] == raw_toplot_dt.shape[0]
    assert ld_toplot_dt.shape[0] == len(labels)
    
    fig, axs = plt.subplots(1,2,figsize = (12,6))
    palette = 'tab10'
    g = sns.scatterplot(ld_toplot_dt[:,0], ld_toplot_dt[:,1],
                    hue = list(labels),
                    ax=axs[0],
                    linewidth=0, s=25, alpha=0.9, palette = palette)
    g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=len(set(labels)))
    
    if ld_pca_evr:
                g.set(title='latent dims (explained variance ratio: {:.2f}'.format(ld_pca_evr) +")")
    else:
        g.set(title='latent dims')
    
    
    g = sns.scatterplot(raw_toplot_dt[:,0], raw_toplot_dt[:,1],
                    hue = list(labels),
                    ax=axs[1],
                    linewidth=0, s=25, alpha=0.9, palette = palette)
    #g.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=1)
    
    if raw_pca_evr:
                g.set(title='full data (explained variance ratio: {:.2f}'.format(raw_pca_evr) +")")
    else:
        g.set(title='full data')
    
    g.legend_.remove()
    
    for ax in axs:
        ax.set_xlabel('{} 1'.format(dr_type))
        ax.set_ylabel('{} 2'.format(dr_type))
    
    fig.suptitle('{}'.format(dr_type), x=0.5, y=0.99, weight="bold")

    if file_name:
        plt.savefig(file_name, dpi=300) 
        print("written : " + file_name)
