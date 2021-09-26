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
# example: Combine multiple OMICs profiles 
############################################

# Get the RNA-Seq, methylation, and mutation profiles for OV cancer patients
seq = getTCGA.getTCGA(disease='OV', dataType='RNASeq2')
meth = getTCGA.getTCGA(disease='OV', dataType='Methylation', type='27K')
mut = getTCGA.getTCGA(disease='OV', dataType='Mutation', type='all')

# Look at the dat matrix
print ( len(seq['dat']) )
print ('seq dat shape:', seq['dat'][0].shape )

print ( len(meth['dat']) )
print ( 'meth dat shape:', meth['dat'][0].shape )

print ( len(mut['dat']) )
print ( 'mut dat shape:', mut['dat'][0].shape )

# Merge the three OMICs-data into one object
m1 = getTCGA.OMICSBind(dat1=seq['dat'][0], dat2=mut['dat'][0])
m2 = getTCGA.OMICSBind(dat1=m1['merged_data'], dat2=meth['dat'][0])

# Look at the keys of all the contents in m1 and m2
print ('m1 keys: %s' %m1.keys() )
print ('m2 keys: %s' %m2.keys() )

# Look at the dimentions of merged_data of m1 and m2
print ('merged_data shape:', m1['merged_data'].shape)
print ('merged dahe shape:', m2['merged_data'].shape)

# Look at the dimension of X and Y of m1 and m2
print ('m1 X shape:', m1['X'].shape)
print ('m1 Y shape:', m1['Y'].shape)

print ('m2 X shape:', m2['X'].shape)
print ('m2 Y shape:', m2['Y'].shape)


