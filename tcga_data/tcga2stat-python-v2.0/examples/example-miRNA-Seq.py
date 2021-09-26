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

# Get miRNA expression data via microarray for ovarian cancer patients
gdats =getTCGA.getTCGA( disease='OV', dataType='miRNASeq', type='count')

# Look at the dat matrix
print ('keys: %s' %gdats.keys() )

print ('# of arrays in dat:', len(gdats['dat']) )
print ('dim of 1st array in dat:',  gdats['dat'][0].shape )

# Look at the first 6 rows and 3 columns
print ('1st 6 rows and 3 columns')
print ( gdats['dat'][0].iloc[0:6,0:3] )

