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

# Get 27K methylation profiles for ovarian cancer patients
gdats = getTCGA.getTCGA(disease='OV', dataType='Methylation', type='27K' )

print ('keys: %s' %gdats.keys() )

print ('# arrays in \'dat\':', len(gdats['dat']) )
print ('# arrays in \'cpgs\':', len(gdats['cpgs']) )
print ('# arrays in \'clinical\':', len(gdats['clinical']) )
print ('# arrays in \'merged_dat\':',len(gdats['merged_dat']) )

print ('shape of 1st array in \'dat\':', gdats['dat'][0].shape )
print ('shape of 1st array in \'cpgs\':', gdats['cpgs'][0].shape )

# Look at the 1st 6 rows and 3 columns
print ('data in dat[0]')
print ( gdats['dat'][0].iloc[0:6,0:3] )

print ('data in cpgs[0]')
print (gdats['cpgs'][0].iloc[0:6,0:3] )

