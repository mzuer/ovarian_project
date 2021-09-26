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

# Get somatic non-silent mutations for ovarian cancer patients
gdats = getTCGA.getTCGA(disease='OV', dataType='Mutation', type='somatic' )

print ('keys: %s' %gdats.keys() )

print ('# arrays in \'dat\':', len(gdats['dat']) )
print ('# arrays in \'clinical\':', len(gdats['clinical']) )
print ('# arrays in \'merged_dat\':',len(gdats['merged_dat']) )

print ('shape of 1st array in \'dat\':', gdats['dat'][0].shape )
# Look at the 1st 6 rows and 6 columns
print ( gdats['dat'][0].iloc[0:6,0:6] )

