# Methylation data were then quantile-normalized across samples (Yang et al. 2019)

# To investigate the feasibility of using variational autoencoder (VAE) models to extract
# biologically meaningful features from DNA methylation data, we trained a VAE model on 
# the top 100,000 most variable CpGs across 1, 230 samples from three publicly available
#  data sets, as defined by median absolute deviation (MAD) of methylation beta values (Table 1).
#  The VAE model learned a 100-dimensional latent representation of the original 
#  100,000-dimensional input data, representing each of the 100,000 original data inputs 
#  as a non-linear combination of the 100 intermediate dimensions. This 100D representation 
#  of the input data is refined by sampling from the 100 dimensions to attempt to reconstruct 
#  the original 100, 000 CpG methylation values. These 100 intermediate dimensions are referred 
#  to here as the latent dimensions of the model, and the trained VAE model is optimized to
#  generate a lower-dimensional representation of the target DNA methylation data.
# Titus et al. 2018


# TCGA-02-0001-01C-01D-0192-01
 ### 0001 is the participant then 01 indicates that it is tumor 
 # (tumor from 01-09; 10-19 normal and control 20-29) (C = the third vial)

# yearstobirth Age (at first diagnosis) Integer 
# vitalstatus Vital status Binary; 1 - dead, 0 - alive/censored
# daystodeath Number of days to death (overall survival time)Integer or NA if the vitalstatus is 0
# daystolastfollowup Number of days to the last follow-up (last known survival time)Integer or NA (sometimes) if the vitalstatus is 1
# gender Gender Categorical; "female" or "male"
# race Race Categorical; "asian", "white", "black or african american", "american indian or alaska native", or "native hawaiian or other pacific islander"
# ethnicity Ethnicity Categorical; "hispanic or latino" or "not hispanic or latino"
# pathologicstage Pathologic stage Categorical; vary slightly based on cancer types but in general could range from Stage I to IV,  and sub-stages such as Stage Ia, Ib, etc.
# pathologyTstage Tumor stage (in TNM staging system) describing the size and location of the tumorCategorical; vary based on cancer types but in general could be "TX", "T0", "Tis",  "T1", "T2", "T3", "T4", and substages like "T2a","T2b" etc.
# pathologyNstage Lymph nodes status (in TNM staging system) describing if the cancer has spread into nearby lymph nodesCategorical; vary based on cancer types but in general could be "NX", "N0", "N1", N2", "N3", and substages such as "N1a", "N2a", etc
# pathologyMstage Metastasis status (in TNM staging system) describing if the cancer has spread to other parts of the bodyCategorical; vary based on cancer types but in general could be "MX", "M0", "M1", and substages "M1a", "M1b", etc.
# 

import sys
sys.path.append('/home/marie/Documents/FREITAS_LAB/ovarian_project/tcga_data/tcga2stat-python-v2.0/')

import getTCGA
#from getTCGA import getRS
from getTCGA.getTCGA import getTCGA
from getTCGA.getRS import getRS
# !!! need to comment the line: 
# from getRS import getRS 
# in getTCGA

import pandas as pd

# Get 450k methylation profiles for ovarian cancer patients
methyl45_ov = getTCGA(disease = "OV", 
                              dataType = 'Methylation',
                              type = "450K",
                              cvars='pathologicstage',
                              clinical = True)
# this downloads a dictionary
# dict_keys(['cpgs', 'dat', 'clinical', 'merged_dat'])

len(methyl45_ov['dat']) # 1
full_dt = methyl45_ov['dat'][0]
full_dt = pd.DataFrame(full_dt)
full_dt.shape
#482421x10

all_samples = full_dt.columns
all_sampleTypes = [x[13:15] for x in all_samples]
np.unique(all_sampleTypes, return_counts=True)
# Out[173]: (array(['01', '02', '11'], dtype='<U2'), array([582,  18,  12]))


full_dt.dtypes
full_dt = full_dt.astype(float)

from scipy.stats import median_abs_deviation

full_dt.apply(median_abs_deviation, axis=1, nan_policy="omit")




methyl27_ov = getTCGA(disease = "OV", 
                              dataType = 'Methylation',
                              type = "27K",
                              clinical = True)

len(methyl27_ov['dat']) # 1
full27_dt = methyl27_ov['dat'][0]
full27_dt.shape
# 27578 x 612
all_samples = full27_dt.columns
all_sampleTypes = [x[13:15] for x in all_samples]
np.unique(all_sampleTypes, return_counts=True)
# Out[173]: (array(['01', '02', '11'], dtype='<U2'), array([582,  18,  12]))

# get the most variant
full27_dt = full27_dt.astype(float)
mad_probes = full27_dt.apply(median_abs_deviation, axis=1, nan_policy="omit")
mad_probes_sorted = mad_probes.sort_values(ascending=False)



len(methyl27_ov['clinical'])
full27_clinical = methyl27_ov['clinical'][0]


set(full27_clinical['tumortissuesite'])
#Out[156]: {nan, 'omentum', 'ovary', 'peritoneum ovary'}

set(full27_clinical['pathologicstage'])
#Out[157]: {nan}

set(full27_clinical['histologicaltype'])
#Out[158]: {nan, 'serous cystadenocarcinoma'}
set(full27_clinical['tumorgrade'])
# Out[160]: {nan}


len(methyl27_ov['merged_dat']) # 1
full27_mdt = methyl27_ov['merged_dat'][0]
full27_dt.shape




