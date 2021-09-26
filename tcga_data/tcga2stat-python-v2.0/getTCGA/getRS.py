import os
import pandas as pd
import numpy as np

from ctypes import *
import ctypes

def getRS(geneInfo, dat):

    # find path to getTCGA
    len1 = len(os.path.realpath(__file__))
    len2 = len('getRS.py')
    tcgaPath = os.path.realpath(__file__)[:len1-len2]

    # load the getmap func
    getmap = CDLL( tcgaPath + 'getmap_1.so' )

    gStart = geneInfo['start'].to_list()
    gEnd   = geneInfo['end'].to_list()
    gChrom = geneInfo['chrom'].to_list()

    dat_grps = dat.groupby('ID')


    gChrom = [ord(item) if item in ['X', 'Y'] else int(item) for item in gChrom]
    chromArray1=(c_char * len(gChrom))( *gChrom )
    startArray1=(c_ulong * len(gStart))( *gStart)
    endArray1=(c_ulong * len(gEnd))( *gEnd )
    length1 = c_int( len(gStart) )

    rs_dict={}
    for grpIdx, grp in dat_grps:

        chrom2 = grp['chrom'].to_list()
        chrom2= [ord(item) if item in ['X', 'Y'] else int(item) for item in chrom2]
        start2 = grp['loc.start'].to_list()
        end2   = grp['loc.end'].to_list()
        values = grp['seg.mean'].to_list()

        chromArray2= (c_char*len(chrom2))( *chrom2 )
        startArray2= (c_ulong * len(start2))( *start2 )
        endArray2=(c_ulong * len(end2))( *end2 )
        valuesArray = (c_double * len(values))( *values )
        length2 = c_int( len(start2) )

        merged = np.ones( len(gStart), dtype=float )
        mergedArray = merged.ctypes.data_as(c_void_p)

        getmap.getmap( chromArray1, startArray1, endArray1, length1,
                       chromArray2, startArray2, endArray2, length2, valuesArray, mergedArray) 

        rs_dict[grpIdx] = merged

    rs_df = pd.DataFrame(rs_dict)

    print ('')
    print ('rs_df shape:', rs_df.shape)
    print ('geneInfo shape:', geneInfo.shape)

    RS_df = pd.concat([geneInfo, rs_df], axis=1, sort=False)
    return RS_df #RS_df

