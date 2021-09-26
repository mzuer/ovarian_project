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


gdats = getTCGA.getTCGA(disease='OV', dataType='RNASeq', type='count' )
#gdats = getTCGA.getTCGA(disease='OV', dataType='RNASeq', type='RPKM' )
#gdats = getTCGA.getTCGA(disease='OV', dataType='RNASeq2' )

print ('keys: %s' %gdats.keys() )

print ('# arrays in \'dat\':', len(gdats['dat']) )
print ('# arrays in \'clinical\':', len(gdats['clinical']) )
print ('# arrays in \'merged_dat\':',len(gdats['merged_dat']) )

print ('shape of 1st array in \'dat\':', gdats['dat'][0].shape )
print ( gdats['dat'][0].iloc[0:5,0:3] )

