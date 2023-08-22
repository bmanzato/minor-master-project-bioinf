# This script creates the Feature Sets files (Genes for scRNA and cCREs for scATAC)
# for both AD and Controls group that will be later used to make the annotation file for LDSC analysis.
# The gene sets are selected with two methods:
    # 1. for each cell type, filter out genes with specificity value 0 and retain the 10% 
    # 2. for each cell type, keep the genes with specificity higher or equal to 0.52 (in this way sets of unique features are created)
# Moreover, the script creates the ENSG_coord_file.txt = bed file of all the genes in the analysis


import pandas as pd
import numpy

########## Create Gene Sets for scRNA sequencing

# Loop over two diagnoses: 'AD' and 'CONTROLS'
for diag in ['AD','CONTROLS']:

    # Read the specificity matrix for the current diagnosis
    spec_matrix = pd.read_csv('/home/bmanzato/ad_data/scrna/specificity_files/specificity_avg_'+diag+'.csv',header=0,index_col=0)

    chrom = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
    # Read gene information from the BioMart file
    mart = pd.read_csv('../mart_export.txt',sep='\t')
    # Filter BioMart data to retain only specified chromosomes
    mart = mart[mart['Chromosome/scaffold name'].isin(chrom)]             

    # Filter Ensembl IDs in the specificity matrix to retain only the ones present in the BioMart file (build 37!)
    overlap = list(set(list(spec_matrix.index)) & set(list(mart['Gene stable ID'])))
    spec_matrix =  spec_matrix[spec_matrix.index.isin(overlap)]
    print('Length of overlap between AnnData object Ensembl IDs ('+diag+') and BioMart Ensembl IDs is ' + str(len(overlap)) + ' genes')


    # Method 1 is retain top 10% most specific genes for each cell type
    # First filter out genes with specificity zero, then take top 10% genes for each cell type.
    cell_types = list(spec_matrix.columns)
    top_genes = dict(zip(cell_types, [None]*len(cell_types)))

    for i in cell_types:
        ct = spec_matrix.sort_values(by=[i], ascending=False)
        ct = ct.loc[ct[i]!=0]
        top = ct.iloc[0:int(ct.shape[0]/10),]
        top = top[i]
        genes = list(top.index)
        top_genes[i] = genes
        save = pd.DataFrame(top.index)
        save.to_csv('/home/bmanzato/ad_data/scrna/genes_files_method1_'+diag+'/'+i+'.tsv',sep="\t",index=False,header=True)



    # Method 2 is select genes whose specificity value is equal or higher than 0.52
    size, avg = [],[]
    for i in cell_types:
        outlist,spec = [],[]
        for j in spec_matrix.index:
            if float(spec_matrix.loc[j,i]) > 0.52:
                outlist.append(j)
                spec.append(spec_matrix.loc[j,i])
        size.append(len(outlist))
        avg.append(numpy.average(spec))
        out = pd.DataFrame(outlist)
        out.to_csv('/home/bmanzato/ad_data/scrna/genes_files_method2_'+diag+'/'+i+'.tsv',sep='\t',header=True,index=False)


    # Print some stats about the number of genes identifies per cell type and their specificity range
    df = spec_matrix.mean()
    df = pd.DataFrame(df,columns=['Mean'])

    top_feat = dict(zip(cell_types, [None]*len(cell_types)))

    for i in cell_types:
        ct = spec_matrix.sort_values(by=[i], ascending=False)
        ct = ct.loc[ct[i]!=0]
        top = ct.iloc[0:int(ct.shape[0]/10),]
        top = top[i]
        feat = list(top.index)
        top_feat[i] = feat
        df.loc[i,'N Genes (Method1)'] = int(len(top))
        df.loc[i,'Mean Specificity (Method1)'] = top.mean()

    # Add Method 2 statistics to the DataFrame
    df['N Genes (Method 2)'] = size
    df['Mean specificity (Method 2)'] = avg
    # Add the number of cells for each cell type
    if diag =='CONTROLS': df['N cells'] = [1760,1644,1806,1384,14993,1038,171]
    else: df['N cells'] = [2996,4725,4156,2742,22059,1702,296]
    df.drop(['Mean'],axis=1,inplace=True)

    print('\nStats about number of genes selected with each method ('+diag+'):')
    print(df)

    # Create ENSG_coord.txt file, will be later needed to create an annotation file in LDSC step make_annot.py
    print('\nCreating gene coordinate file for '+diag)
    gene_coord =  mart[mart['Gene stable ID'].isin(overlap)]
    gene_coord = gene_coord[['Gene stable ID','Chromosome/scaffold name','Gene start (bp)','Gene end (bp)']]
    gene_coord.rename(columns={'Gene stable ID': 'GENE', 'Chromosome/scaffold name': 'CHR',
                               'Gene start (bp)': 'START', 'Gene end (bp)': 'END'}, inplace=True)
    gene_coord.iloc[0,:] = ['GENE','CHR','START','END']
    print('...')
    gene_coord.to_csv('/home/bmanzato/ad_data/scrna/coord_files/ENSG_coord_'+diag+'.tsv',sep="\t",index=False,header=False)
    print('Done creating ENSG coordinate file\n')



########## Create cCREs Sets for scATAC sequencing

for diag in ['AD', 'Controls']:
    spec_matrix = pd.read_csv('/home/bamnzato/ad_data/scatac/specificity_files/specificity_avg_' + diag + '.csv',
                              header=0, index_col=0)

    df = spec_matrix.mean()
    df = pd.DataFrame(df, columns=['Mean Specificity (incl 0)'])

    cell_types = ['ASC', 'EX', 'INH', 'MG', 'ODC', 'OPC', 'PER.END']
    top_feat = dict(zip(cell_types, [None] * len(cell_types)))

    # Method 1
    for i in cell_types:
        ct = spec_matrix.sort_values(by=[i], ascending=False)
        ct = ct.loc[ct[i] != 0]
        top = ct.iloc[0:int(ct.shape[0] / 10), ]
        top = top[i]
        feat = list(top.index)
        top_feat[i] = feat
        df.loc[i, 'Mean Specificity (excl 0)'] = top.mean()
        df.loc[i, 'Top 10% Max'] = top.max()
        df.loc[i, 'Top 10% Min'] = top.min()
        df.loc[i, 'N Top 10% Feat'] = len(top)
        top = pd.DataFrame(top.index)
        top.to_csv('/home/bmanzato/ad_data/scatac/peaks_files_method1_' + diag + '/' + i + '.tsv', sep="\t",
                   index=False, header=True)

    # Method 2
    size, avg = [], []
    for i in cell_types:
        outlist, spec = [], []
        for j in spec_matrix.index:
            if float(spec_matrix.loc[j, i]) > 0.52:
                outlist.append(j)
                spec.append(spec_matrix.loc[j, i])
        size.append(len(outlist))
        avg.append(numpy.average(spec))
        out = pd.DataFrame(outlist)
        out.to_csv('/home/bmanzato/ad_data/scatac/peaks_files_method2_' + diag + '/' + i + '.tsv', sep='\t',
                   header=True, index=False)
