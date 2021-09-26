import os
import sys
sys.path.append('./getTCGA')

import getTCGA


def main():

#   test = "RNASeq2"
#   test = 'RNASeq'
#   test = 'miRNASeq'
#   test = 'Methylation'
#   test = 'Mutation'
#   test = 'mRNA_Array'
#   test = 'miRNA_Array'
#   test = 'CNV_SNP'
#   test = 'CNA_SNP'
#   test = 'CNA_CGH'
#   test = 'SampleSplit'
#   test = 'TumorNormalMatch'
    test = 'OMICSBind'
    SAVE_FILE = 'No'
#   SAVE_FILE = 'Yes'


    gdats=[]
    if test=='RNASeq2':
        dataset = 'GBM'
        dataType = 'RNASeq2'
        type = ''
        filter = ['Y']
        cvars = 'OS'
        filter = ['Y']
        gdats = getTCGA.getTCGA(dataset, dataType, type, True, cvars, filter)


    if test == 'RNASeq':
        dataset = 'BRCA'
        dataType = 'RNASeq'
        type = 'RPKM'
        filter = ['Y']
        cvars = 'OS'
        filter = ['Y']
        gdats = getTCGA.getTCGA(dataset, dataType, type, True, cvars, filter)


    if test == 'miRNASeq':
        dataset = 'ACC'
        dataType = 'miRNASeq'
        type = 'read_count'
        filter = ['Y']
        cvars = 'OS'
        filter = ['Y']
        gdats = getTCGA.getTCGA(dataset, dataType, type, True, cvars, filter)


    if test == 'Methylation':
        dataset = 'KIRP'
        dataType = 'Methylation'
        type = '27K'
        filter = ['Y']
        cvars = 'OS'
        filter = ['Y']
        gdats = getTCGA.getTCGA(dataset, dataType, type, True, cvars, filter)

        print ('shape of gdats[\'dat\'][0]:', gdats['dat'][0].shape)
        print ('shape of gdats[\'cpgs\'][0]:', gdats['cpgs'][0].shape)
        print ('shape of gdats[\'clinical\'][0]:', gdats['clinical'][0].shape)
        print ('shape of gdats[\'merged_dat\'][0]:', gdats['merged_dat'][0].shape)
        print ('')

        name = './results/' + dataType + '_' + dataset + '_' + type + '_dat'        + '_p' + '.csv'
        print ( name )
        gdats['dat'][0].to_csv(name, index=True)
        name = './results/' + dataType + '_' + dataset + '_' + type + '_cpgs'       + '_p' + '.csv'
        print ( name )
        gdats['cpgs'][0].to_csv(name, index=True)
        name = './results/' + dataType + '_' + dataset + '_' + type + '_clinical'   + '_p' + '.csv'
        print ( name )
        gdats['clinical'][0].to_csv(name, index=True)
        name = './results/' + dataType + '_' + dataset + '_' + type + '_merged_dat' + '_p' + '.csv'
        print ( name )
        gdats['merged_dat'][0].to_csv(name, index=True)
        exit()


    if test == 'Mutation':
        dataset = 'ACC'
        dataType = 'Mutation'
        type = 'somatic'
        filter = ['Y']
        cvars = 'OS'
        filter = ['Y']
        gdats = getTCGA.getTCGA(dataset, dataType, type, True, cvars, filter)

    if test == 'mRNA_Array':
        dataset = 'BRCA'
        dataType = 'mRNA_Array'
        type = 'G450'
        filter = ['Y']
        cvars = 'OS'
        filter = ['Y']
        gdats = getTCGA.getTCGA(dataset, dataType, type, True, cvars, filter)

    if test == 'miRNA_Array':
        dataset = 'GBM'
        dataType = 'miRNA_Array'
        type = ''
        filter = ['Y']
        cvars = 'OS'
        filter = ['Y']
        gdats = getTCGA.getTCGA(dataset, dataType, type, True, cvars, filter)

    if test == 'CNV_SNP':
        dataset = 'ACC'
        dataType = 'CNV_SNP'
        type = ''
        filter = ['Y']
        cvars = 'OS'
        filter = ['Y']
        gdats = getTCGA.getTCGA(dataset, dataType, type, True, cvars, filter)

    if test == 'CNA_SNP':
        dataset = 'ACC'
        dataType = 'CNA_SNP'
        type = ''
        filter = ['Y']
        cvars = 'OS'
        filter = ['Y']
        gdats = getTCGA.getTCGA(dataset, dataType, type, True, cvars, filter)

    if test == 'CNA_CGH':
        dataset = 'GBMLGG'
        dataType = 'CNA_CGH'
        type = '244A'
        filter = ['Y']
        cvars = 'OS'
        filter = ['Y']
        gdats = getTCGA.getTCGA(dataset, dataType, type, True, cvars, filter)

    if test == 'SampleSplit':
        gdats = getTCGA.getTCGA(disease='LUSC', dataType='RNASeq2')
        tumors = getTCGA.SampleSplit( gdats['dat'][0] )
        print ('dim dat: ', gdats['dat'][0].shape)
        print ('dim primary tumor:', tumors['primary_tumor'].shape)
        print ('dim recurrent tumor:', tumors['recurrent_tumor'].shape )
        print ('dim normal:', tumors['normal'].shape)
        exit()

    if test == 'TumorNormalMatch':
        gdats = getTCGA.getTCGA(disease='LUSC', dataType='RNASeq2')
        tumors = getTCGA.TumorNormalMatch( gdats['dat'][0] )
        print ('dim dat: ', gdats['dat'][0].shape)
        print ('dim primary tumor:', tumors['primary_tumor'].shape)
        print ('dim normal:', tumors['normal'].shape)
        exit()

    if test == 'OMICSBind':
        seq = getTCGA.getTCGA(disease='OV', dataType='RNASeq2')
        meth = getTCGA.getTCGA(disease='OV', dataType='Methylation', type='27K')
        mut = getTCGA.getTCGA(disease='OV', dataType='Mutation', type='all')
        m1 = getTCGA.OMICSBind( dat1 = seq['dat'][0], dat2 = mut['dat'][0] )
        m2 = getTCGA.OMICSBind( dat1=m1['merged_data'], dat2=meth['dat'][0] )
        print ('seq dat shape:', seq['dat'][0].shape)
        print ('meth dat shape:', meth['dat'][0].shape)
        print ('mut data shape:', mut['dat'][0].shape)
        print ('ml dim:', m1['merged_data'].shape )
        print ('m2 dim:', m2['merged_data'].shape )
        exit()

    if SAVE_FILE == 'Yes':
        print ('')
        print ('len of gdats[\'dat\']:', len(gdats['dat']) )
        print ('shape of gdats[\'dat\'][0]:', gdats['dat'][0].shape)
        print ('shape of gdats[\'clinical\'][0]:', gdats['clinical'][0].shape)
        print ('shape of gdats[\'merged_dat\'][0]:', gdats['merged_dat'][0].shape)
        print ('')
        name = './results/' + dataType + '_' + dataset + '_' + type + '_dat' + '_p' + '.csv'
        print ( name )
        gdats['dat'][0].to_csv(name, index=True)

        name = './results/' + dataType + '_' + dataset + '_' + type + '_clinical' + '_p' + '.csv'
        print ( name )
        gdats['clinical'][0].to_csv(name, index=True)

        name = './results/' + dataType + '_' + dataset + '_' + type + '_merged_dat' + '_p' + '.csv'
        print ( name )
        gdats['merged_dat'][0].to_csv(name, index=True)



if __name__ == "__main__":
    main()

