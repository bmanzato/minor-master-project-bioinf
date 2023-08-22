# pwd : morabitoetal/scatac/scripts

# This script:
#    - reads the original Cases-Controls H5 Matrix
#    - separates AD and Controls
#    - filters out regions with less than 5 cumulative peak intensities in all cells
#    - write the new control and metadata matrices to csv files
#    - average the peak intensities per cell-type
#    - performs LIFTOVER creating 'COORD_file_ALL.tsv' = file with index the 'old' regions (build 38) and
#      the corresponding build 37 coordinates

import pandas as pd
import collections
import scipy.sparse as sp_sparse
import tables
from collections import Counter
import numpy as np
import scanpy as sc
import numbers
import liftover


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


# Read the .h5 file
filtered_matrix_h5_atac = "/home/bmanzato/ad_data/scatac/data/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5"
filtered_feature_bc_matrix_atac = get_matrix_from_h5(filtered_matrix_h5_atac)

# Convert the Sparse matrix to Dense
matrix = filtered_feature_bc_matrix_atac[2]
dense_matrix = matrix.todense()

pd_matrix = pd.DataFrame(dense_matrix)

features = []
barcodes = []

for i in list(filtered_feature_bc_matrix_atac[0]['name']):
    j = str(i)[2:-1]
    features.append(j)

for i in list(filtered_feature_bc_matrix_atac[1]):
    j = str(i)[2:-1]
    barcodes.append(j)

pd_matrix.index = features
pd_matrix.columns = barcodes

# Filter to keep Control individuals only 
for i in ['AD','Controls']:
    atac_metadata = pd.read_csv('/home/bmanzato/ad_data/scatac/data/GSE174367_snATAC-seq_cell_meta.csv')
    meta = atac_metadata.loc[atac_metadata['Diagnosis']==i]
    counts = pd_matrix[list(meta['Barcode'])]

    print('size of AD scATAC before filtering the regions:')
    print(counts.shape)

    # Filter out regions that have a peak intensity (sum) count lower than 5
    counts['sum'] = counts.where(counts < 1, 1).sum(axis=1)
    counts = counts.loc[counts['sum'] > 5]

    counts.drop(['sum'],axis=1,inplace=True)
    print('Shape of the preprocessed Peak Intensity Count Matrix:')
    print(counts.shape)
    # Save
    counts.to_csv('/home/bmanzato/ad_data/scatac/data/'+i+'.csv')
    meta.to_csv('/home/bmanzato/ad_data/scatac/data/METADATA_'+i+'.csv')

    # Average per cell-type
    counts.columns = meta['Cell.Type']
    df2 = controls.transpose()
    df2 = df2.groupby(by=df2.index, axis=0).apply(lambda g: g.mean() if isinstance(g.iloc[0,0], numbers.Number) else g.iloc[0])
    df = df2.transpose()

    df.to_csv('/home/bmanzato/ad_data/scatac/data/'+i+'_avg.csv')


# LiftOver step

converter = liftover.get_lifter('hg38', 'hg19')
skipped = []
keeps = []
features.rename(columns={0: 'GRCh38'},inplace=True)
features.index = features['GRCh38']

# Create Chr Start and End columns from the Specificity Matrix row names (chr:start-end)
features[['Chr', 'Pos']] = features['GRCh38'].str.split(":", expand = True)
#features['Chr'] = features['Chr'].str.replace('chr', '')
features[['Start', 'End']] = features['Pos'].str.split("-", expand = True)
features.drop(['Pos'],axis=1,inplace=True)
chrom, start, end = features['Chr'], features['Start'], features['End']

chromkeep = ['chr1', 'chr2', 'chr3', 'chr4','chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
             'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']
features = features[features['Chr'].isin(chromkeep)]
outlist = features


for i in features.index:
    x = converter[features.loc[i,'Chr']][int(features.loc[i,'Start'])]
    if len(x) == 0:
        n += 1
        skipped.append(features.loc[i,'GRCh38'])
    else:
        chrom = x[0][0]
        start = x[0][1]
        keeps.append(features.loc[i,'GRCh38'])


    x = converter[features.loc[i,'Chr']][int(features.loc[i,'End'])]
    if len(x) == 0:
        n += 1
        skipped.append(features.loc[i,'GRCh38'])
    else:
        end = x[0][1]
        keeps.append(features.loc[i,'GRCh38'])

    # if the chomosome changes in the conversion remove the region
    if chrom != features.loc[i,'Chr']:
        skipped.append(features.loc[i,'GRCh38'])

    length = abs(int(features.loc[i,'Start']) - int(features.loc[i,'End']))
    new_length = abs(int(start) - int(end))
    # check that the length of the GRCh37 region is not longer than the GRCh38 region +/- 20%
    if length-0.2*length < new_length and length+0.2*length > new_length:
        if int(start) < int(end):
            outlist.loc[i,'Chr'] = chrom
            outlist.loc[i,'Start'] = start
            outlist.loc[i,'End'] = end
        else:
            outlist.loc[i,'Chr'] = chrom
            outlist.loc[i,'Start'] = end
            outlist.loc[i,'End'] = start

    else:
        skipped.append(features.loc[i,'GRCh38'])

skipped = list(set(skipped))
keeps = list(set(keeps))


print(f"total amount of regions liftover: {len(keeps)}")
print(f"total amount of regions that failed liftover: {len(skipped)}")

outlist = pd.DataFrame(outlist)
outlist['keep'] = outlist.index
outlist = outlist[['keep','Chr','Start','End']]
outlist.rename({'keep': 'GENE', 'Chr': 'CHR', 'Start': 'START', 'End': 'END'}, axis=1, inplace=True)
outlist.drop(skipped,axis=0,inplace=True)


outlist.to_csv('/home/bmanzato/ad_data/scatac/coord_files/COORD_file_ALL.tsv',sep='\t',header=True,index=False)