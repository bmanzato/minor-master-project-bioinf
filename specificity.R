# Script to compute cell-type specificity; for both single cell technologies (scatac)
# and both diagnosis groups (AD and controls)


library(EWCE)
library(tidyverse)
library(ggdendro)
library(dplyr)
library(purrr)



for (diag in c('AD','CONTROLS')) {
  for (sc in c('scrna','scatac')){
    
    # Take the count matrix avg within a cell type group
    file = paste('/home/bmanzato/ad_data/',sc,'/data/',diag,'_avg.csv',sep='')
	  exp <- read_csv(file)
	
	  feat <- exp[,1]
	  exp <- exp[,-1]
	  exp <- as.matrix(exp)
	  rownames(exp) <- feat[[1]]
	
	  annotlev <- list(l1 = colnames(exp))
	  specif <- generate_celltype_data(exp,
					  groupName = paste(diag,'_avg',sep=''),
					  annotLevels = annotlev,
					  savePath = '/home/bmanzato/ad_data/',sc,'/specificity_files')
	
	  load(file=paste('/home/bmanzato/ad_data/',sc,'/specificity_files/CellTypeData_',diag,'_avg.rda',sep=''))
	  c <- ctd[[1]]$specificity
	  c <- as.matrix(c)
	  write.csv(c,paste('/home/bmanzato/ad_data/',sc,'/specificity_files/specificity_avg_',diag,'.csv',sep=''))
  }
}


