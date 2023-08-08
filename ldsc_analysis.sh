#!/bin/bash

#Set job requirements





# Setting up the conda environment
# First get access to conda (it is not transferred to the nodes)
source /gpfs/admin/_hpc/sw/arch/AMD-ZEN2/RHEL8/EB_production/2022/software/Mamba/4.14.0-0/etc/profile.d/conda.sh
# Deactivate current environments, otherwise locations to packages may not be not right; then load the snakemake environment
conda deactivate
conda activate ldsc


# Create SNP set:
for i in {1..22};
	do gunzip -c ../../LDSC_inputfiles/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.${i}.l2.ldscore.gz | tail -n +2 | awk '{ print $2 }' > ../../LDSC_inputfiles/Baseline_v2.2_snplists/snplist_chrom_${i}.snp; done

# Annotation for the feature Sets:

for diag in AD CONTROLS;
do for method in method1 method2;
do for CT in PER.END OPC ODC MG ASC EX INH;
do for chrom in {1..22};

        #scrna
        do python ../make_annot.py --gene-set-file ../../ad_data/scrna/genes_files_${method}_${diag}/${CT}.tsv \
        --gene-coord-file ../../ad_data/scrna/coord_files/ENSG_coord_${diag}.txt \
        --windowsize 100000 \
        --bimfile ../../LDSC_inputfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
        --annot-file ../../ldsc_rna/${diag}/${method}/${CT}.${chrom}.gz &

        #scatac
        python ../make_annot.py --gene-set-file ../../ad_data/scatac/peak_files_${method}_${diag}/${CT}.tsv \
        --gene-coord-file ../../ad_data/scatac/coord_files/COORD_file_${diag}.txt \
        --windowsize 0 \
        --bimfile ../../LDSC_inputfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
        --annot-file ../../ldsc_atac/${diag}/${method}/${CT}.${chrom}.annot.gz &

        #overlap set between scrna and scatac (annotation from bed file)
        python ../make_annot.py --bed-file ../../ad_data/overlap_analysis/overlap_${diag}_${method}_${CT}.txt \
        --bimfile ../../LDSC_inputfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
        --annot-file ../../ldsc_overlap/${diag}/${method}/${CT}.${chrom}.annot.gz &

        #merge
        python ../make_annot.py --bed-file ../../ad_data/overlap_analysis/merge_${diag}_${method}_${CT}.txt \
        --bimfile ../../LDSC_inputfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
        --annot-file ../../ldsc_merge/${diag}/${method}/${CT}.${chrom}.annot.gz &
done;
wait
done; done; done


# Annotation for all the features in the initial dataset

for diag in AD CONTROLS;
do for chrom in {1..22};

	#scrna
         do python ../make_annot.py --gene-set-file ../../ad_data/scrna/coord_files/ENSG_coord_${diag}.txt \
         --gene-coord-file ../../ad_data/scrna/coord_files/ENSG_coord_ALL.txt \
         --windowsize 100000 \
         --bimfile ../../LDSC_inputfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
         --annot-file ../../ldsc_rna/${diag}/ALL/all.${chrom}.annot.gz &
        
	#scatac
        python ../make_annot.py --gene-set-file ../../ad_data/scatac/coord_files/COORD_file_${diag}.txt \
         --gene-coord-file ../../ad_data/scatac/coord_files/COORD_file_ALL.txt \
         --windowsize 0 \
         --bimfile ../../LDSC_inputfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
         --annot-file ../../ldsc_atac/${diag}/ALL/all.${chrom}.annot.gz &
        
	#overlap and merge sets 
        python ../make_annot.py --bed-file ../../ad_data/overlap_analysis/merge_coord_${diag}.txt \
        --bimfile ../../LDSC_inputfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
        --annot-file ../../ldsc_overlap/${diag}/ALL/all.${chrom}.annot.gz;
done; done


# Compute LD scores

for sc in overlap merge scatac scrna;
do for CT in ASC EX INH MG OPC ODC PER.END;
do for method in method1 method2;
do for diag in AD CONTROLS;
do for chrom in {1..22};

        do python ../ldsc.py --l2 --bfile ../../LDSC_inputfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
		--ld-wind-cm 1 \
		--annot ../../ldsc_${sc}/${diag}/${method}/${CT}.${chrom}.annot.gz \
		--thin-annot \
		--out ../../ldsc_${sc}/${diag}/${method}/${CT}.${chrom} \
		--print-snps ../../LDSC_inputfiles/Baseline_v2.2_snplists/snplist_chrom_${chrom}.snp &
done; done;
wait
done; done; done


# Compute LD scores for annotation of all the features

for sc in overlap scrna scatac;
do for diag in AD CONTROLS;
do for chrom in {1..22};

	do python ../ldsc.py --l2 --bfile ../../LDSC_inputfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
		--ld-wind-cm 1 \
		--annot ../../ldsc_${sc}/${diag}/ALL/all.${chrom}.annot.gz \
		--thin-annot \
		--out ../../ldsc_${sc}/${disg}/ALL/all.${chrom} \
		--print-snps ../../LDSC_inputfiles/Baseline_v2.2_snplists/snplist_chrom_${chrom}.snp & 
done; done
wait
done


# Run LDSC Regression

for diag in AD CONTROLS;
do for sc in overlap merge scrna scatac;
do for method in method1 method2;
        
	do python ../ldsc.py --h2-cts ../../munged_sumstats/AD.sumstats.gz \
		--ref-ld-chr ../../LDSC_inputfiles/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD. \
		--out ../../LDSC_output/${sc}_${diag}_${method} \
		--ref-ld-chr-cts ../../LDSC_inputfiles/LDCTS_files/${sc}_${diag}_${method}.ldcts \
		--w-ld-chr ../../LDSC_inputfiles/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. &
done; done;
wait
done

