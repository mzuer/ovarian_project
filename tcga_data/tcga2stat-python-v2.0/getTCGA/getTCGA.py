import requests
import urllib.request
import time
#from beautifulsoup4 import BeautifulSoup
from bs4 import BeautifulSoup

import tarfile
import re
import pandas as pd
import numpy as np
import os
import time
import collections
import math


def get_all_values_TagA( url = 'http://gdac.broadinstitute.org/runs/stddata__latest/' ):

    r = requests.get(url, allow_redirects=True)
    soup = BeautifulSoup(r.text, "html.parser")
    valuesTagA = soup.findAll('a')

    return valuesTagA

def extract_dlinks( valuesTagA, dataset ):
    return [msgTagA['href'] for msgTagA in valuesTagA if '/data/' + dataset +'/' in msgTagA['href'] ]


def extract_plinks ( valuesTagA, keyWord ):
    return [msgTagA.string for msgTagA in valuesTagA if keyWord in msgTagA['href'] and msgTagA.string[-7:]=='.tar.gz' ]


def file_download_cleanup( dlink, plink, dataset, fileDest, newNameStr, dlFileNamePattern ):

    download_link = dlink + '/' + plink

#   print ('download_link')
#   print (download_link)

    # download the tar data file
    r = requests.get(download_link, allow_redirects=True)
    gzFile = fileDest + '/' + dataset + newNameStr + '.tar.gz'
    open(gzFile, 'wb').write(r.content)
 
    # get the list of desired files from the tar file
    myTar = tarfile.open(gzFile, 'r')
    myTarList = myTar.getmembers()
    grepSearchStr = "." + dataset + dlFileNamePattern

    fileList = [member.name for member in myTarList if re.search(grepSearchStr, member.name)]

    # extract file
    #    fileList[0] for now. It looks like the tar.gz file only contains one .txt file
    fname = fileDest + '/' + fileList[0]
    myTar.extract(fileList[0], path = fileDest)
    myTar.close()

    # time_stamp 
    time_stamp = dlink.split('/')[-1]

    # read in the extracted file
    newName =  fileDest + '/' + dataset + '_' + time_stamp + newNameStr + '.txt'
    os.rename(fname, newName)
    df = pd.read_csv(newName, sep='\t', low_memory=False)

    # cleanup
    os.remove(gzFile)
    os.remove(newName)
    dirToDelete = fileDest + '/' + fileList[0].split('/')[0]
    os.rmdir(dirToDelete)

    return df


def file_download_cleanup_mutation( dlink, plink, dataset, fileDest, newNameStr, dlFileNamePattern ):

    download_link = dlink + '/' + plink

    # download data file
    r = requests.get(download_link, allow_redirects=True)
    gzFile = fileDest + '/' + dataset + newNameStr + '.tar.gz'
    open(gzFile, 'wb').write(r.content)

    myTar = tarfile.open(gzFile, 'r')
    myTar_list = myTar.getmembers()

    grepSearchStr = dlFileNamePattern

    fileList=[]
    ldIndex = 0
    for member in myTar_list:
        match = re.search( grepSearchStr, member.name)
        if match:
            fileList.append( member.name )
            myTar.extract( member, path = fileDest )
            fname = fileDest + '/' + member.name
            if ldIndex == 0:
                df= pd.read_csv(fname, sep='\t', low_memory=False)
            else:
                df= pd.concat( [df, pd.read_csv(fname, sep='\t', low_memory=False)] )
            ldIndex += 1
            os.remove(fname)
    myTar.close()

    # cleanup
    dir_to_delete = fileDest + '/' + fileList[0].split('/')[0]
    os.rmdir(dir_to_delete)
    os.remove(gzFile)

    return df



def processDataFrame_RNASeqV2( df ):

    colNames = list(df.columns.values)
    colsToKeep = []
    for k in range( len(df.iloc[0].tolist()) ):
        if df.iloc[0].tolist()[k] == "normalized_count":
            colsToKeep.append( df.columns.values[k] )

    # drop 1st row which contains 'normalized_count' 
    df = df.iloc[1:].reset_index(drop=True)

    col_one = df[ colNames[0] ].to_list()
    df =df[colsToKeep]

    # set new index
    df.index = col_one 

    # extract gen names from 1st column
    genNames = [ s.split('|')[0] for s in col_one ]

    # replace df index with gen names
    df.index = genNames

    # eliminate rows with row names either be "?" or duplicate
    #    find rows with index values == '?'
    test_list_0 = [s=='?' for s in genNames]
    #     find row indexes that duplicate
    test_list_1 = df.index.duplicated(keep='first')
    #     indexes either with "?" or duplicate
    test_list_1 = test_list_0 | test_list_1
    #     eliminate rows with row names either be "?" or duplicate
    df = df[~test_list_1]

    return df


def processDataFrame_RNASeq( df, type ):

    colNames = list(df.columns.values)
    colsToKeep = []

    for k in range( len(df.iloc[0].tolist()) ):
        if df.iloc[0].tolist()[k] == type:
            colsToKeep.append( df.columns.values[k] )

    # drop 1st row which contains 'normalized_count' 
    df = df.iloc[1:].reset_index(drop=True)

    col_one = df[ colNames[0] ].to_list()
    df =df[colsToKeep]

    # set new index
    df.index = col_one 

    # extract gen names from 1st column
    genNames = [ s.split('|')[0] for s in col_one ]

    # replace df index with gen names
    df.index = genNames

    # eliminate rows with row names either be "?" or duplicate
    #    find rows with index values == '?'
    test_list_0 = [s=='?' for s in genNames]
    #     find row indexes that duplicate
    test_list_1 = df.index.duplicated(keep='first')
    #     indexes either with "?" or duplicate
    test_list_1 = test_list_0 | test_list_1
    #     eliminate rows with row names either be "?" or duplicate
    df = df[~test_list_1]

    # replace '.' with '-' in column names
    colNames = list(df.columns.values)
    colNames = [s.replace('.', '-') for s in colNames]
    df.columns = colNames
    return df


def processDataFrame_miRNASeq( df, type ):

    colNames = list(df.columns.values)
    colsToKeep = []
    for k in range( len(df.iloc[0].tolist()) ):
        if df.iloc[0].tolist()[k] == type:
            colsToKeep.append( df.columns.values[k] )

    # drop 1st row which contains 'normalized_count' 
    df = df.iloc[1:].reset_index(drop=True)

    col_one = df[ colNames[0] ].to_list()
    df =df[colsToKeep]

    # set new index
    df.index = col_one 

    # find row indexes that duplicate
    test_list_1 = df.index.duplicated(keep='first')
    # eliminate duplicate rows
    df = df[~test_list_1]

    # replace '.' with '-' in column names
    colNames = list(df.columns.values)
    colNames = [s.replace('.', '-') for s in colNames]
    df.columns = colNames
    return df


def processDataFrame_Methylation( df, type ):

    colNames = list(df.columns.values)
    row_1st = df.iloc[0].tolist()
 
    # keep columns with the type == type
    colsToKeep = [i for (i,j) in zip(colNames, row_1st) if j==type]
    colsToKeep = [df.columns.values[2],df.columns.values[3],df.columns.values[4]] + colsToKeep
    newNames_3cols = [row_1st[2], row_1st[3], row_1st[4]]

    # drop 1st row which contains type info 
    df = df.iloc[1:].reset_index(drop=True)

    # use 1st column as row names
    col_1st = df[ colNames[0] ].to_list()

    # eliminates columns that type!=type
    df =df[colsToKeep]

    # set new row index
    df.index = col_1st 

    # rename the first 3 column names
    colNames = list(df.columns.values)
    colNames[0:3] = newNames_3cols

    # replace '.' with '-' in column names
    colNames = [s.replace('.', '-') for s in colNames]
    df.columns = colNames

    # keep only cpg probes
    df = df[df.index.str.contains('^cg')]

    return df


def processDataFrame_mRNA_Array( df, type ):

    colNames = list(df.columns.values)
    colsToKeep = colNames[1:]

    # drop 1st row which contains 'normalized_count' 
    df = df.iloc[1:].reset_index(drop=True)

    col_one = df[ colNames[0] ].to_list()
    df =df[colsToKeep]

    # set new index
    df.index = col_one 

    # replace '.' with '-' in column names
    colNames = list(df.columns.values)
    colNames = [s.replace('.', '-') for s in colNames]
    df.columns = colNames
    return df

def processDataFrame_miRNA_Array( df, type ):

    colNames = list(df.columns.values)
    colsToKeep = colNames[1:]

    # drop 1st row which contains 'normalized_count' 
    df = df.iloc[1:].reset_index(drop=True)

    col_one = df[ colNames[0] ].to_list()
    df =df[colsToKeep]

    # set row index
    df.index = col_one 

    # replace '.' with '-' in column names
    colNames = list(df.columns.values)
    colNames = [s.replace('.', '-') for s in colNames]
    df.columns = colNames

    # check and remove control probes
    findStr = '(^dmr)|(^NC)|(Corner)|(Control)|(^hur)|(^mr)'
    rowNames = list(df.index.values)

    rowNamesToDrop = [name for name in rowNames if re.search(findStr,name)]
    df = df.drop(rowNamesToDrop)

    return df

def processDataFrame_Clinical( df ):

    colNames = list(df.columns.values)
    colsToKeep = colNames[1:]

    col_one = df[ colNames[0] ].to_list()
    df =df[colsToKeep]

    # transpose dataframe
    df = df.T

    colNames = [s.replace('_', '') for s in col_one]
    rowNames = [s.upper() for s in colsToKeep]
    rowNames = [s.replace('.', '-') for s in rowNames]

    df.columns = colNames
    df.index = rowNames

    return df





def RNASeqV2(dlinks, dataset, dataType):

    keyWord = 'Level_3__RSEM_genes_normalized__data.Level_3'
    fileDest = '.'
    newNameStr = '-RNAseq2GeneNorm'
    dlFileNamePattern = "[.]rnaseqv2__.*.__Level_3__RSEM_genes_normalized__data.data.txt$"

    # get values of all Tag \a from download link for dataset
    valuesTagA = get_all_values_TagA( dlinks[0] )

    # get the download link for the actual data set
    plinks = extract_plinks( valuesTagA, keyWord )

    plink_str = dataset + '.Merge_rnaseqv2__illuminahiseq*._rnaseqv2__.*.tar[.]gz$'
    plinks = [plink for plink in plinks if re.search( plink_str, plink) ]

    gdats = []
    for k in range(len(plinks)):
        # formulate fownload link
        dlink = dlinks[0]
        plink = plinks[k]
        df = file_download_cleanup( dlink, plink, dataset, fileDest, newNameStr, dlFileNamePattern )

        # post process data frame
        gdat = processDataFrame_RNASeqV2( df )

        gdats.append( gdat )

    return gdats

def RNASeq(dlinks, dataset, dataType, type):

    if type=='count':
        type='raw_counts'

    keyWord = 'Level_3__gene_expression__data.Level_3'
    fileDest = '.'  
    newNameStr = '-RNAseqGene' # .tar.gz'
    dlFileNamePattern = '[.]rnaseq__.*.__Level_3__gene_expression__data.data.txt$'
 
    # get values for the Tag a from download link for dataset
    valuesTagA = get_all_values_TagA( dlinks[0] )

    # get the download link for the actual data set
    plinks = extract_plinks( valuesTagA, keyWord )

    plink_str = dataset + '*.Merge_rnaseq__.*._rnaseq__.*.tar[.]gz$'
    plinks = [plink for plink in plinks if re.search( plink_str, plink) ]

    if len( plinks ) > 1:
        [plink for plink in plinks if "*_illuminahiseq_*" in plink ]


    gdats = []
    for k in range(len(plinks)):
        # formulate fownload link
        dlink = dlinks[0]
        plink = plinks[k]
        df = file_download_cleanup( dlink, plink, dataset, fileDest, newNameStr, dlFileNamePattern )

        # post process data frame
        gdat = processDataFrame_RNASeq( df, type )

        gdats.append( gdat )

    return gdats

# type could be "count" and "rpmmm"
def miRNASeq(dlinks, dataset, dataType, type='count'):

    if(type not in ["count", "rpmmm"]):
        print ('Error: invalid error type')

    if type=='count':
        type = 'read_count'

    if type=='rpmmm':
        type = 'reads_per_million_miRNA_mapped'

    keyWord = 'Level_3__miR_gene_expression__data.Level_3'
    fileDest = '.'
    newNameStr = '-miRNAseqGene' # .tar.gz'
    dlFileNamePattern = '[.]mirnaseq__.*.__Level_3__miR_gene_expression__data.data.txt$'

    # get values for the Tag a from download link for dataset
    valuesTagA = get_all_values_TagA( dlinks[0] )

    # get the download link for the actual data set
    plinks = extract_plinks( valuesTagA, keyWord )

    plink_str = '.' + dataset + '[.]Merge_mirnaseq__.*.hiseq_mirnaseq__.*.tar[.]gz$'
    plinks = [plink for plink in plinks if re.search( plink_str, plink) ]

    if len( plinks ) > 1:
        [plink for plink in plinks if "*_illuminahiseq_*" in plink ]

    print ('download files')
    gdats = []
    for k in range(len(plinks)):
        # formulate fownload link
        dlink = dlinks[0]
        plink = plinks[k]
        df = file_download_cleanup( dlink, plink, dataset, fileDest, newNameStr, dlFileNamePattern )

        # post process data frame
        gdat = processDataFrame_miRNASeq( df, type )

        gdats.append( gdat )

    return gdats


# type: "27K", "450K"
def Methylation(dlinks, dataset, dataType, type='read_count'):

    keyWord = '__Level_3__within_bioassay_data_set_function__data.Level_3'
    fileDest = '.' 
    newNameStr = '-Methylation.tar'
    dlFileNamePattern = '[.]methylation__.*.__Level_3__within_bioassay_data_set_function__data.data.txt$'

    # get values for the Tag a from download link for dataset
    valuesTagA = get_all_values_TagA( dlinks[0] )

    # get the download link for the actual data set
    plinks = extract_plinks( valuesTagA, keyWord )


    if type == '27K':
        methy_str = '[.]Merge_methylation__.*.methylation27.*.__Level_3__within_bioassay_data_set_function__data.Level_3.*.tar[.]gz$'
    if type == '450K':
        methy_str = '[.]Merge_methylation__.*.methylation450.*.__Level_3__within_bioassay_data_set_function__data.Level_3.*.tar[.]gz$'

    plink_str = '.' + dataset + methy_str
    plinks = [plink for plink in plinks if re.search( plink_str, plink) ]

    if len( plinks ) > 1:
        [plink for plink in plinks if "*_illuminahiseq_*" in plink ]


    print ('download files')
    gdats = []
    for k in range(len(plinks)):
        # formulate fownload link
        dlink = dlinks[0]
        plink = plinks[k]
        df = file_download_cleanup( dlink, plink, dataset, fileDest, newNameStr, dlFileNamePattern )

        # post process data frame
        gdat = processDataFrame_Methylation( df, 'Beta_value' )

        gdats.append( gdat )

    return gdats


def Proc_Mutation( dat ):

    dat.drop_duplicates()
    dat['Tumor_Sample_Barcode'] = dat.apply(lambda x: x['Tumor_Sample_Barcode'][:12],axis=1)
    dat.drop_duplicates()

    genes = list(dat['Hugo_Symbol'].drop_duplicates())
    samples = list(dat['Tumor_Sample_Barcode'].drop_duplicates())

    dat_dict={}
    for sample in samples:
        gene_list =  list(dat[dat['Tumor_Sample_Barcode'] == sample]['Hugo_Symbol'])  
        dat_dict[sample] = [1 if gene in gene_list else 0 for gene in genes ]

    dat = pd.DataFrame(dat_dict)
    dat.index = genes

    return dat



# Mutations type allows "all" or "somatic"
def Mutation(dlinks, dataset, dataType, type="somatic", method="auto" ):

    keyWord = 'Mutation_Packager_Calls'
    fileDest = '.'
    newNameStr = '-Mutation'
    dlFileNamePattern = '.maf.txt$'

    # get values for the Tag a from download link for dataset
    valuesTagA = get_all_values_TagA( dlinks[0] )

    # get the download link for the actual data set
    plinks = extract_plinks( valuesTagA, keyWord )

    plink_str = '.' + dataset + '[.]Mutation_Packager_Calls[.]Level_3[.].*.tar[.]gz$'
    plinks = [plink for plink in plinks if re.search( plink_str, plink) ]

    print ('download files')
    gdats = []
    mafs = []
    for k in range(len(plinks)):

        # download files and clean up
        gdat = file_download_cleanup_mutation( dlinks[0], plinks[k], dataset, fileDest, newNameStr, dlFileNamePattern )

        if type == 'somatic':
            gdat = gdat[ (gdat['Mutation_Status']=='Somatic') & (gdat['Variant_Classification']!='Silent') ]

        gdat = Proc_Mutation( gdat )

        gdats.append( gdat )

    return gdats


def mRNA_Array(dlinks, dataset, dataType, type='G450'):

    if type not in ['G450', 'U133', 'Huex']:
        print ( 'Error: Invalid miRNA-array platform.' )
        return []

    fileDest = '.'
    newNameStr = '-mRNAArray'
    dlFileNamePattern = '.*__data.data.txt$'

    if type == 'G450':
        keyWord1 = 'Merge_transcriptome__agilentg4502a_07'
        keyWord2 = '[.]Merge_transcriptome__agilentg4502a_.*.__Level_3__unc_lowess_normalization_gene_level__data.Level_3.*.tar[.]gz$' 
    if type == 'U133':
        keyWord1 = 'Merge_transcriptome__ht_hg_u133a'
        keyWord2 = '[.]Merge_transcriptome__ht_hg_u133a__.*.__Level_3__gene_rma__data.Level_3.*.tar[.]gz$'
    if type == 'Huex':
        keyWord1 = 'Merge_exon__huex_1_0_st_v2'
        keyWord2 = '[.]Merge_exon__huex_1_0_st_v2__.*.__Level_3__quantile_normalization_gene__data.Level_3.*.tar[.]gz$'

    # get values for the Tag a from download link for dataset
    valuesTagA = get_all_values_TagA( dlinks[0] )
    
    # get the download link for the actual data set
    plinks = extract_plinks( valuesTagA, keyWord1 )

    plink_str = '.' + dataset + keyWord2 
    plinks = [plink for plink in plinks if re.search( plink_str, plink) ]

    if len( plinks ) == 0:
        print('Error: No data available for download. Please ensure this data is available from TCGA.')
        return []

    print ('download files')
    dats = []
    for k in range(len(plinks)):
        # download file
        df = file_download_cleanup( dlinks[0], plinks[k], dataset, fileDest, newNameStr, dlFileNamePattern )

        # post process data frame
        gdat = processDataFrame_mRNA_Array( df, type )

        dats.append( gdat )

    if len(dats) > 1:
        for k in range(1,len(dats)):
            dats[0] = pd.concat( [dats[0], dats[k]], axis=1 )

    gdats = []
    gdats.append( dats[0] )

    return gdats



def miRNA_Array(dlinks, dataset, dataType, type=''):

    fileDest = '.'
    newNameStr = '-miRNAArray'
    dlFileNamePattern = '.*__data.data.txt$'

    keyWord1 = 'h_mirna_8x15k'
    keyWord2 = '.Merge_mirna__h_mirna_8x15k.*.data.Level_3.*.tar[.]gz$'

    # get values for the Tag a from download link for dataset
    valuesTagA = get_all_values_TagA( dlinks[0] )
    
    # get the download link for the actual data set
    plinks = extract_plinks( valuesTagA, keyWord1 )

    plink_str = keyWord2 
    plinks = [plink for plink in plinks if re.search( plink_str, plink) ]

    if len( plinks ) == 0:
        print('Error: No data available for download. Please ensure this data is available from TCGA.')
        return []

    print ('download files')
    dats = []
    for k in range(len(plinks)):
        # download file
        df = file_download_cleanup( dlinks[0], plinks[k], dataset, fileDest, newNameStr, dlFileNamePattern )

        # post process data frame
        gdat = processDataFrame_miRNA_Array( df, type )
        dats.append( gdat )

    return dats


#from getRS import getRS

def getSegmentsByGene( dat, geneInfo, filter=['Y']):

    if isinstance(filter,list)==False:
        print ('Expect filter is a list')
        exit (1)

    dat = dat[~dat['chrom'].isin( filter )]
    rs = getRS(geneInfo, dat)
    rs =  rs[ ~rs['chrom'].isin( filter ) ]

    return rs


def CNA_SNP(dlinks, dataset, dataType, geneInfo, filter ):

    fileDest = '.'
    newNameStr = '-CNASNPHg19'
    dlFileNamePattern = '[.]snp__.*.__Level_3__segmented_scna_hg19__seg.seg.txt$'

    keyWord1 = 'Level_3__segmented_scna_hg19__seg.Level_3'
    keyWord2 = '[.]Merge_snp__.*.__Level_3__segmented_scna_hg19__seg.Level_3.*.tar[.]gz$'

    # get values for the Tag a from download link for dataset
    valuesTagA = get_all_values_TagA( dlinks[0] )
    
    # get the download link for the actual data set
    plinks = extract_plinks( valuesTagA, keyWord1 )

    plink_str = '.' + dataset + keyWord2 
    plinks = [plink for plink in plinks if re.search( plink_str, plink) ]

    if len( plinks ) == 0:
        print('Error: No data available for download. Please ensure this data is available from TCGA.')
        exit(1) 

    print ('download files')
    dats = []
    for k in range(len(plinks)):
        # download file
        dat = file_download_cleanup( dlinks[0], plinks[k], dataset, fileDest, newNameStr, dlFileNamePattern )
        dat.columns = ['ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']

        rs = getSegmentsByGene ( dat, geneInfo, ['Y'] )

        geneSymbols = rs['genename'].to_list()
        chroms = rs['chrom'].to_list()

        dupGeneSymbols = [x for x, y in collections.Counter( geneSymbols ).items() if y > 1]
        rchr = ['_' + chrom if symbol in dupGeneSymbols else '' for symbol, chrom in zip(geneSymbols, chroms)]
        geneSymbols = [s + c for s, c in zip(geneSymbols,rchr)]

        rs = rs[list(rs.columns.values)[5:]]
        rs.index = geneSymbols

        dats.append( rs )

    return dats



def CNV_SNP(dlinks, dataset, dataType, geneInfo, filter=['Y']):


    fileDest = '.' 
    newNameStr = '-CNVSNPHg19'
    dlFileNamePattern = '[.]snp__.*.Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt$'

    keyWord1 = 'Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3'
    keyWord2 = '[.]Merge_snp__.*.Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.*.tar[.]gz$'

    # get values for the Tag a from download link for dataset
    valuesTagA = get_all_values_TagA( dlinks[0] )
    
    # get the download link for the actual data set
    plinks = extract_plinks( valuesTagA, keyWord1 )

    plink_str = '.' + dataset + keyWord2 
    plinks = [plink for plink in plinks if re.search( plink_str, plink) ]

    if len( plinks ) == 0:
        print('Error: No data available for download. Please ensure this data is available from TCGA.')
        exit(1) 

    print ('download files')
    dats = []
    for k in range(len(plinks)):
        # download file
        dat = file_download_cleanup( dlinks[0], plinks[k], dataset, fileDest, newNameStr, dlFileNamePattern )
        dat.columns = ['ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']

        rs = getSegmentsByGene ( dat, geneInfo, ['Y'] )

        geneSymbols = rs['genename'].to_list()
        chroms = rs['chrom'].to_list()

        dupGeneSymbols = [x for x, y in collections.Counter( geneSymbols ).items() if y > 1]
        rchr = ['_'+chrom if symbol in dupGeneSymbols else '' for symbol, chrom in zip(geneSymbols, chroms)]
        geneSymbols = [s + c for s, c in zip(geneSymbols,rchr)]

        rs = rs[list(rs.columns.values)[5:]]
        rs.index = geneSymbols

        dats.append( rs )

    return dats


def CNA_CGH(dlinks, dataset, dataType, type, geneInfo, filter=['Y']):

    if type not in ['415K', '244A']:
        print ('Error: Invalid type for CGH data.')
        exit (1)


    fileDest = '.'
    newNameStr = '-CNACGH'
    dlFileNamePattern = '[.]cna__.*.__Level_3__segmentation__seg.seg.txt$'

    keyWord1 = '__Level_3__segmentation__seg.Level_3'
    if type == '415K':
        keyWord2 = '[.]Merge_cna__.*.cgh_415k.*.__Level_3__segmentation__seg.Level_3.*.tar[.]gz$'
    if type == '244A':
        keyWord2 = '[.]Merge_cna__.*.cgh_244a.*.__Level_3__segmentation__seg.Level_3.*.tar[.]gz$'

    # get values for the Tag a from download link for dataset
    valuesTagA = get_all_values_TagA( dlinks[0] )
    
    # get the download link for the actual data set
    plinks = extract_plinks( valuesTagA, keyWord1 )

    plink_str = '.' + dataset + keyWord2 
    plinks = [plink for plink in plinks if re.search( plink_str, plink) ]

    if len( plinks ) == 0:
        print('Error: No data available for download. Please ensure this data is available from TCGA.')
        exit(1) 

    print ('download files')
    dats = []
    for k in range(len(plinks)):
        # download file
        dat = file_download_cleanup( dlinks[0], plinks[k], dataset, fileDest, newNameStr, dlFileNamePattern )
        dat.columns = ['ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean']

        rs = getSegmentsByGene ( dat, geneInfo, ['Y'] )

        geneSymbols = rs['genename'].to_list()
        chroms = rs['chrom'].to_list()

        dupGeneSymbols = [x for x, y in collections.Counter( geneSymbols ).items() if y > 1]
        rchr = ['_' + chrom if symbol in dupGeneSymbols else '' for symbol, chrom in zip(geneSymbols, chroms)]
        geneSymbols = [s + c for s, c in zip(geneSymbols,rchr)]

        rs = rs[list(rs.columns.values)[5:]]
        rs.index = geneSymbols

        dats.append( rs )

    return dats




def Clinical(dlinks, dataset):

    fileDest = '.'
    newNameStr = '-Clinical.tar.gz'
    dlFileNamePattern = '.clin.merged.picked.txt$'

    keyWord1 = '.Clinical_Pick_Tier1.Level_4'
    keyWord2 = '.tar[.]gz$'

    # get values for the Tag a from download link for dataset
    valuesTagA = get_all_values_TagA( dlinks[0] )
    
    # get the download link for the actual data set
    plinks = extract_plinks( valuesTagA, keyWord1 )

    plink_str = keyWord2 
    plinks = [plink for plink in plinks if re.search( plink_str, plink) ]

    if len( plinks ) == 0:
        print('Error: No Clinical data available for download. Please ensure this data is available from TCGA.')
        return []

    gdats = []
    for k in range(len(plinks)):
        # download file
        df = file_download_cleanup( dlinks[0], plinks[k], dataset, fileDest, newNameStr, dlFileNamePattern )

        # post process data frame
        cli = processDataFrame_Clinical( df )

        gdats.append( cli )

    return gdats

#def overallSurvial (ctbl ):
def survival(vitalstatus, daystodeath, daystolastfollowup):

    if vitalstatus=='1':
        x = daystodeath
    else:
        x = daystolastfollowup
    return int(x)

def MatrixMerge( dats, cli_list, cvars ):

    print ('cvars:', cvars)

    if (len(cli_list)==0):
        print ('Error: no clinical data available for this disease type.')
        return []


    cli = cli_list[0]
    colNames = list(cli.columns.values)
    colNamesUpcase = [elm.upper() for elm in colNames]
    colNamesUpcasePlusOS = colNamesUpcase + ['OS']

    if cvars.upper() not in colNamesUpcasePlusOS:
        print ('Error: invalid clinical covariates.')
        return []

    if cvars.upper() == 'OS':
        ctbl = cli[ ['vitalstatus', 'daystodeath', 'daystolastfollowup'] ].copy()
    else:
        cli.columns = colNamesUpcase
        ctbl = cli[cvars.upper()].copy().to_frame()

    dat_clis = []
    for dat in dats:
        colNames = list( dat.columns.values )

        temp = dat

        # for methylation
        if colNames[:3] == ['Gene_Symbol', 'Chromosome', 'Genomic_Coordinate']:
            temp = dat[colNames[3:]]

        sum_tl = sum( [True if len(name) ==12 else False for name in colNames ] ) 
        if sum_tl == len(colNames):
            tum_pat = temp
        else:
            tum = [ s.split('-')[3] for s in list(temp.columns.values) ]
            tum_pat = temp[ [col for (s, col) in zip(tum, colNames) if re.search('^01',s)] ]
            tum_pat.columns = [ s[:12] for s in list(tum_pat.columns.values)]

        if cvars.upper() != 'OS':
            temp1 = ctbl.copy()
        else:
            lstatus = list(ctbl['vitalstatus'])

            ldtd = list( ctbl['daystodeath'] )
            ldtd = [str(l) for l in ldtd]
            ldtd = [int(elem) if elem!='nan' else elem for elem in ldtd]

            ldtlf = list( ctbl['daystolastfollowup'] )
            ldtlf = [str(l) for l in ldtlf]
            ldtlf = [int(elem) if elem!='nan' else elem for elem in ldtlf]

            ctbl['OS'] = [ dtd if status=='1' else dtlf for status, dtd, dtlf in zip(lstatus, ldtd, ldtlf) ]
            ctbl['status'] = [int(status) for status in lstatus]

            temp1 = ctbl[ ['status', 'OS'] ].copy()
            temp1.index = list(ctbl.index.values)

        temp2 = tum_pat.T

        temp1.insert(0, 'bcr', list(temp1.index))
        temp2.insert(0, 'bcr', list(temp2.index))

        dat_cli = pd.merge(temp1, temp2, on='bcr', how='inner')

    dat_clis.append(dat_cli) 
    return dat_clis



def getTCGA( disease='GBM', dataType='RNASeq2', type='', clinical=False, cvars='OS', filter=['Y'] ): 

    dataset = disease

    print ('getTCGA - ', dataType)
    print ('dataset:', dataset)
    print ('dataType:', dataType)
    print ('type:', type)


    # find path to getTCGA
    len1 = len(os.path.realpath(__file__))
    len2 = len('getTCGA.py')
    geneInfoDir = os.path.realpath(__file__)[:len1-len2]
    print (geneInfoDir)

    # get values for all the Tag a
    valuesTagA = get_all_values_TagA()

    # extract download links
    dlinks = extract_dlinks( valuesTagA, dataset )

    gdats = {} 

    start = time.time()
    if dataType == 'RNASeq2':
        dat = RNASeqV2( dlinks, dataset, dataType)
        gdats['dat'] = dat
    end = time.time()

    start = time.time()
    if dataType == 'RNASeq':
        dat = RNASeq(dlinks, dataset, dataType, type )
        gdats['dat'] = dat
    end = time.time()

    if dataType == 'miRNASeq':
        dat = miRNASeq(dlinks, dataset, dataType, type)
        gdats['dat'] = dat

    if dataType == 'Methylation':
        dats = Methylation(dlinks, dataset, dataType, type)
        dat = dats[0]
        colNames = list(dat.columns.values)
        gdats={}
        gdats['cpgs'] = []
        gdats['cpgs'].append( dat[colNames[0:3]] )
        gdats['dat'] = []
        gdats['dat'].append( dat[colNames[3:]] )
        print ('dat shape:', gdats['dat'][0].shape)
        print ('cpgs shape:', gdats['cpgs'][0].shape)

        gdats['clinical'] = []
        gdats['merged_dat'] = []
        if clinical==True:
            cli = Clinical( dlinks, dataset )
            mdat = MatrixMerge( gdats['dat'], cli, cvars )
            gdats['clinical'] = cli
            gdats['merged_dat'] = mdat
        return gdats

    if dataType == 'Mutation':
        dat = Mutation(dlinks, dataset, dataType, type)
        gdats['dat'] = dat

    if dataType == 'mRNA_Array':
        dat = mRNA_Array(dlinks, dataset, dataType, type)
        gdats['dat'] = dat

    if dataType == 'CNV_SNP':
        geneInfo = pd.read_csv( geneInfoDir + '/geneinfo.csv')
        dat = CNV_SNP(dlinks, dataset, dataType, geneInfo, filter )
        gdats['dat'] = dat

    if dataType == 'CNA_SNP':
        geneInfo = pd.read_csv( geneInfoDir + '/geneinfo.csv')
        dat = CNA_SNP(dlinks, dataset, dataType, geneInfo, filter )
        gdats['dat'] = dat

    if dataType == 'CNA_CGH':
        geneInfo = pd.read_csv( geneInfoDir + '/geneinfo.csv')
        dat = CNA_CGH(dlinks, dataset, dataType, type, geneInfo, filter )
        gdats['dat'] = dat

    if dataType == 'miRNA_Array':
        dat = miRNA_Array(dlinks, dataset, dataType, type)
        gdats['dat'] = dat

    gdats['clinical'] = []
    gdats['merged_dat'] = []
    if clinical==True:
        cli = Clinical( dlinks, dataset )
        mdat = MatrixMerge( dat, cli, cvars )
        gdats['clinical'] = cli
        gdats['merged_dat'] = mdat

    return gdats

def SampleSplit( dat ):

    tumors={}
    colnames = list(dat.columns.values)
    if colnames[0]=='Gene_Symbol' and colnames[1]=='Chromosome' and colnames[2]=='Genomic_Coordinate':
        tempCols = dat[colnames[3:]]

#       colnames_temp = list( temp.columns.values )
        tempTypes = [ s.split('-')[3] for s in tempCols ]

        cols_annot = colnames[:3]
        cols_pTumor = [ colname for (colname, coltype) in zip(tempCols, tempTemp) if re.search('^01', coltype) ]
        primary_tumor = dat[ cols_annot + cols_pTumor ]
 
        cols_rcurrTumor = [ colname for (colname, coltype) in zip(tempCols, tempTemp) if re.search('^02', coltype) ]
        recurrent_tumor = dat[ cols_annot + cols_recurrTumor ]

        cols_10 = [ colname for (colname, coltype) in zip(tempCols, tempTemp) if re.search('^10', coltype) ]
        cols_11 = [ colname for (colname, coltype) in zip(tempCols, tempTemp) if re.search('^11', coltype) ]
        cols_12 = [ colname for (colname, coltype) in zip(tempCols, tempTemp) if re.search('^12', coltype) ]
        normal = dat[ cols_annot + cols_10 + cols_11 + cols_12]
    else:
        tempTypes = [ s.split('-')[3] for s in colnames ]

        cols_pTumor = [ colname for (colname, coltype) in zip(colnames, tempTypes) if re.search('^01', coltype) ]
        primary_tumor = dat[ cols_pTumor ]
 
        cols_rcurrTumor = [ colname for (colname, coltype) in zip(colnames, tempTypes) if re.search('^02', coltype) ]
        recurrent_tumor = dat[ cols_rcurrTumor ]

        cols_10 = [ colname for (colname, coltype) in zip(colnames, tempTypes) if re.search('^10', coltype) ]
        cols_11 = [ colname for (colname, coltype) in zip(colnames, tempTypes) if re.search('^11', coltype) ]
        cols_12 = [ colname for (colname, coltype) in zip(colnames, tempTypes) if re.search('^12', coltype) ]
        normal = dat[ cols_10 + cols_11 + cols_12]


    tumors['primary_tumor'] = primary_tumor
    tumors['recurrent_tumor'] = recurrent_tumor
    tumors['normal'] = normal
    return tumors

def TumorNormalMatch( dat ):

    colNames = list(dat.columns.values)
    tempTypes = [ s.split('-')[3] for s in colNames ]

    cols_pTumor = [ colname for (colname, coltype) in zip(colNames, tempTypes) if re.search('^01', coltype) ]
    primary_tumor = dat[ cols_pTumor ]
 
    cols_10 = [ colname for (colname, coltype) in zip(colNames, tempTypes) if re.search('^10', coltype) ]
    cols_11 = [ colname for (colname, coltype) in zip(colNames, tempTypes) if re.search('^11', coltype) ]
    cols_12 = [ colname for (colname, coltype) in zip(colNames, tempTypes) if re.search('^12', coltype) ]
    normal = dat[ cols_10 + cols_11 + cols_12]

    tum_bcr = [ s[:12] for s in list(primary_tumor.columns.values) ]
    norm_bcr = [ s[:12] for s in list(normal.columns.values) ]
    matching = list( set(tum_bcr).intersection( set(norm_bcr) ) )

    tumors = {}
    if len(matching)==0:
        tumors['primary_tumor'] = []
        tumors['normal'] = []
        print ('No matching samples on tumor and normal.')
        return tumors
  
    if len(matching)>0:
        primary_tumor.columns = tum_bcr
        normal.columns = norm_bcr
        tumors['primary_tumor'] = primary_tumor[matching]
        tumors['normal'] = normal[matching]
        return tumors

def OMICSBind(dat1, dat2):

    colnames1 = list(dat1.columns.values)
    colnames2 = list(dat2.columns.values)

    dat1_tumor = dat1
    dat2_tumor = dat2

    sum1 = sum( [True if len(name) >12 else False for name in colnames1 ] )
    sum2 = sum( [True if len(name) >12 else False for name in colnames2 ] )

    dat1_tumor = dat1
    if sum1== len(colnames1):
        dat1Types = [ s.split('-')[3] for s in colnames1 ]
        colsTumor1 = [ colname for (colname, coltype) in zip(colnames1, dat1Types) if re.search('^01', coltype) ]
        dat1_tumor = dat1[colsTumor1] 
        dat1_tumor.columns = [ s[:12] for s in colsTumor1 ]

    dat2_tumor = dat2
    if sum2== len(colnames2):
        dat2Types = [ s.split('-')[3] for s in colnames2 ]
        colsTumor2 = [ colname for (colname, coltype) in zip(colnames2, dat2Types) if re.search('^01', coltype) ]
        dat2_tumor = dat2[colsTumor2] 
        dat2_tumor.columns = [ s[:12] for s in colsTumor2 ]

    matching = list( set(dat1_tumor).intersection( set(dat2_tumor) ) )

    tumors={}
    if len(matching)==0:
        tumros['mergedData'] = []
        tumors['dat1Good'] = []
        tumors['dat2Good'] = []
        return tumors

    dat1_good = dat1_tumor[matching]
    dat2_good = dat2_tumor[matching]

    rownames1 = ['d1.'+name for name in list(dat1_good.index.values) ]
    rownames2 = ['d2.'+name for name in list(dat2_good.index.values) ]

    dat1_good.index = rownames1
    dat2_good.index = rownames2

    mdata = pd.concat([dat1_good, dat2_good])
    tumors['merged_data'] = mdata
    tumors['X'] = dat1_good
    tumors['Y'] = dat2_good
    return tumors
