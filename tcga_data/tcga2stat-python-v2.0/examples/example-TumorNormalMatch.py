import os
import sys

#
# ATTENTATION: 
# (1) uncomment the following line
# (2) add the PATH of directory, getTCGA, in sys.path.append()
#sys.path.append('directory_of_getTCGA')

try:
    import getTCGA
except ModuleNotFoundError as err:
    print (err)
    print ('Please provide path for getTCGA by adding a line \"sys.path.append(\'directory_of_getTCGA\')\"')
    print ('    at the top of this script')
    exit()

############################################
# example: Tumor Normal Match 
############################################

# Get the RNA-SeqV2 data for LUSC patients
lusc_ranseq2 = getTCGA.getTCGA(disease='LUSC', dataType='RNASeq2')
#Tumor-normal matched profiles
lusc_ranseq2_tum_norm = getTCGA.TumorNormalMatch( lusc_ranseq2['dat'][0] )

# Look at the keys of all the contents in m1 and m2
print ('lusc_ranseq2 keys: %s' %lusc_ranseq2.keys() )
print ('lusc_ranseq2_tum_norm: %s' %lusc_ranseq2_tum_norm.keys() )

# Look at the dimentions of merged_data of m1 and m2
print ('dim dat: ', lusc_ranseq2['dat'][0].shape)
print ('dim primary tumor:', lusc_ranseq2_tum_norm['primary_tumor'].shape)
print ('dim normal:', lusc_ranseq2_tum_norm['normal'].shape)

