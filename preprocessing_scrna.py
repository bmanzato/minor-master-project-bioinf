# pwd : morabitoetal/scrna/scripts

# This scripts:
#    - reads the count matrix in the H5 file
#    - separates AD and Controls
#    - filters out genes which are detected in less than 5 cells
#    - averages all the cells belonging to the same cell type (size=N cell types x N genes) 
#    - writes the new count matrix, the new metadata and the averaged count matrix

import h5py
import pandas as pd
import collections
import scipy.sparse as sp_sparse
import tables
import numpy as np
import scanpy as sc
import numbers


CountMatrix = collections.namedtuple('CountMatrix', ['feature_ref', 'barcodes', 'matrix'])

def get_matrix_from_h5(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)

        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read()
        feature_names = getattr(feature_group, 'name').read()
        feature_types = getattr(feature_group, 'feature_type').read()
        feature_ref['id'] = feature_ids
        feature_ref['name'] = feature_names
        feature_ref['feature_type'] = feature_types
        tag_keys = getattr(feature_group, '_all_tag_keys').read()
        for key in tag_keys:
            key = key.decode("utf-8")
            feature_ref[key] = getattr(feature_group, key).read()

        return CountMatrix(feature_ref, barcodes, matrix)


filtered_matrix_h5_rna = "/home/bmanzato/ad_data/scrna/data/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5"
filtered_feature_bc_matrix_rna = get_matrix_from_h5(filtered_matrix_h5_rna)

matrix = filtered_feature_bc_matrix_rna[2]
dense_matrix = matrix.todense()

# Assign feature name and barcodes
pd_matrix = pd.DataFrame(dense_matrix)

features = []
barcodes = []

for i in list(filtered_feature_bc_matrix_rna[0]['id']):
    j = str(i)[2:-1]
    features.append(j)

for i in list(filtered_feature_bc_matrix_rna[1]):
    j = str(i)[2:-1]
    barcodes.append(j)

pd_matrix.index = features
pd_matrix.columns = barcodes

print('Size of the Whole Count Matrix is '+str(len(features))+' genes and '+str(len(barcodes))+' cells')

rna_metadata = pd.read_csv('/home/bmanzato/ad_data/scrna/data/GSE174367_snRNA-seq_cell_meta.csv',header=0)
print(rna_metadata)
rna_metadata = rna_metadata.sort_values(by = 'Barcode')

for i in ['AD','Control']:

    # Filter Diagnosis groups and keep only genes detected in at least 5 cells 
    metagroup = rna_metadata.loc[rna_metadata['Diagnosis']==i]
    counts = pd_matrix[list(metagroup['Barcode'])]
    counts2 = counts.where(counts < 1, 1)
    counts['sum'] = counts2.sum(axis=1)
    counts3 = counts.loc[counts['sum'] > 4]
    counts3.drop(['sum'],axis=1,inplace=True)
    print('Size of the '+i+' Count Matrix is '+str(counts.shape[0])+' genes and '+str(counts.shape[1])+' cells')

    # remove version from ENSG gene name: e.g. ENSG00398492.2 -> remove .2
    counts3['genes'] = spec_matrix.index
    counts3[['ensembl', 'none']] = counts3['genes'].str.split(".", expand = True)
    counts3.index = spec_matrix['ensembl']
    counts3.drop(['genes','ensembl','none'],axis=1,inplace=True)

    # Save the averaged count matrix and the metadata
    counts3.to_csv('/home/bmanzato/ad_data/scrna/data/'+i+'.csv')
    metagroup.to_csv('/home/bmanzato/ad_data/scrna/data/METADATA_'+i+'.csv')
    
    # Average all cells belonging to the same Cell Type
    counts.columns = metagroup['Cell.Type']
    df2 = counts3.transpose()
    df2 = df2.groupby(by=df2.index, axis=0).apply(lambda g: g.mean() if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[0])
    df = df2.transpose()
    df.to_csv('/home/bmanzato/ad_data/scrna/data/'+i+'_avg.csv')


