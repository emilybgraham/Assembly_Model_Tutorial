#####written September 5, 2025 as a null model tutorial####
#loosely based on: 
#Knelman, Joseph E., et al. "Multiple, compounding disturbances in a forest ecosystem: fire increases susceptibility of soil edaphic properties, bacterial community structure, and function to change with extreme precipitation event." Soil Systems 3.2 (2019): 40.

#NOTE: this code uses rarefaction and OTUs, but any normalization and species-grouping will work. The tree should be trimmed to the OTUs or ASVs that are in the final table.

#####set up####
rm(list=ls());graphics.off()#clear all
#change line below to directory where your files are
setwd('/Users/grah930/Documents/WSU_fire_paper/assembly_model_code_and_tutorial')
#load libraries
library(permute);library(gee);library(vegan);library(ape);library(picante);library(ecodist);library(picante);
#set rarefaction depth and dataset name
rare.depth = 10730 #this parameter should be adjusted to balance keeping samples and them having enough reads; lowest recommended is 5-10k, but it could be much higher (~50)
data.set.name = 'Knelman_2019'

#reps are set to 100 for testing; this should be much higher, recommended 999
no.reps = 100#number of reps for null models

#load functions, note these can be altered if needed
source("Raup_Crick_Abundance_One_Comparison.R")
source("beta.nti.R")

#####Optional Pre-processing####
#for this tutorial, we are using an OTU table that is already rarefied and a tree that has been trimmed to match.
#if rarefaction and tree trimming has not been done, example code is below.

#otu = as.data.frame(read.table(paste(data.set.name,"_otu_table.txt",sep=""),sep="\t",header=T,row.names=1,skip = 1))
#phylo = read.tree(paste(data.set.name,"_final.tre",sep=""))

#examine OTU table:
# otu[1:5,1:5]
# dim(otu)
# 
## number of sequences per sample; may have to change 1 to 2
# sample.seq.counts = apply(otu,1,sum)
# range(sample.seq.counts)

## number of sequences per otu; may have to change 1 to 2
# otus.seq.counts = apply(otu,2,sum)
# length(otus.seq.counts)
# range(otus.seq.counts)

#drop singletons
# otus.to.drop = names(otus.seq.counts)[otus.seq.counts < 2]
# length(otus.to.drop)
# # 233
# 
# otu.abund = otu[-which(rownames(otu) %in% otus.to.drop),]
# dim(otu.abund)
# # 7599   81

#number of sequences per otu
# otus.seq.counts = apply(otu.abund,1,sum)
# length(otus.seq.counts)
# range(otus.seq.counts)
# # 2 271023
#rm('otu')

# #set rarefaction depth
# seqs.per.samp = apply(otu.abund,2,sum)
# hist(seqs.per.samp); range(seqs.per.samp)
# rare.depth = 13000 # 

#drop samples with less than rarefied depth
# samps.to.drop = names(seqs.per.samp)[seqs.per.samp < rare.depth]

# #rarefy and drop any OTUs that aren't present in rarefied table
# otu.rare = rrarefy(t(otu.samp.drop),rare.depth)
# min.abund.after.rare = 1
# otus.seq.counts = apply(otu.rare,2,sum)
# otus.to.drop = names(otus.seq.counts)[otus.seq.counts < min.abund.after.rare]
# otu.rare = otu.rare[,-which(colnames(otu.rare) %in% otus.to.drop)]
# seqs.per.samp = apply(otu.rare,1,sum)
# range(seqs.per.samp)#range should be rarefaction depth (both min and max)

# write.csv(otu.rare,paste(data.set.name,"_OTU_rarefied_",unique(seqs.per.samp),".csv",sep=""),quote=F)
# otu = otu.rare
# rm('otu.rare','otu.abund','otu.samp.drop','seqs.per.samp')

#trim tree
# match.phylo.otu = match.phylo.data(phylo, as.data.frame(t(otu))); # species as rows, samples as columns for otu table
# str(match.phylo.otu); 

####Read in OTU table, trimmed tree, set number of samples####
otu = t(read.csv(paste(data.set.name,"_OTU_rarefied_",rare.depth,".csv",sep=""),row.names=1))
colnames(otu) = gsub(pattern = "X",replacement = "",x = colnames(otu))

phylo = read.tree(paste(data.set.name,"_tree_matched_to_rarified_otu_table.tre",sep=""))
phylo$tip.label= gsub("'","",phylo$tip.label)# if needed, removing extra notion in tree$tip.label
match.phylo.otu = match.phylo.data(phylo, as.data.frame(t(otu))); # if needed, match tree to OTU table, species as rows, samples as columns for otu table

no.of.samples = nrow(otu)

####### Step 1: calculate observed betaMNTD #####
#this can be computationally intensive
date(); beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T)); date();
#dim(beta.mntd.weighted);beta.mntd.weighted[1:5,1:5];
write.csv(beta.mntd.weighted,paste(data.set.name,"_bMNTD_weighted.csv",sep=""),quote=F);

##### Step 2: calculate randomized betaMNTD #####
rep = 1:no.reps

for (i in rep){
  print(i)
  rand.weighted.beta.mntd = as.matrix(comdistnt(comm = t(match.phylo.otu$data), 
                                                dis = taxaShuffle(cophenetic(match.phylo.otu$phy)), 
                                                abundance.weighted = T)); ## randomized beta.mntd
  write.csv(rand.weighted.beta.mntd,paste(data.set.name ,"_bMNTD_weighted_rep_",rep[i],".csv",sep=""),quote=F);
  print(date());
  rm(rand.weighted.beta.mntd)
}

##### Step 3: calculate betaNTI (must have observed beta.mntd and all null model reps run) #####
beta.mntd.weighted = read.csv(paste(data.set.name,"_bMNTD_weighted.csv",sep=""),row.names=1);

#get all null model files
all.files = sort(as.numeric(gsub(".*_rep_|\\.csv", "", list.files("./", pattern = paste0(data.set.name, "_bMNTD_weighted_rep_\\d+\\.csv$")))))

## find beta.mntd reps that still need to be done; if there are any, need to do them before proceeding
reps.to.do = c(1:no.reps)[-which(1:no.reps %in% all.files)];reps.to.do#should be empty
#write.table(t(reps.to.do),"./reps.to.do.txt",sep=" ",quote=F,row.names=F,col.names=F)

beta.nti.weighted = beta.nti.calc.stegen(samp = otu, reps = all.files, 
                                         path.to.reps = paste('./',data.set.name,'_bMNTD_weighted_rep_',sep=""),
                                         beta.mntd.obs = beta.mntd.weighted);
write.csv(beta.nti.weighted,paste(data.set.name,'_bNTI_weighted.csv',sep=""),quote=F);

# optional plot to check distribution
# pdf(paste(data.set.name,'_bNTI_weighted.pdf',sep=""));
# hist(as.dist(beta.nti.weighted),xlim = c(min(c(-2,min(as.dist(beta.nti.weighted)))),max(c(2,max(as.dist(beta.nti.weighted))))),xlab="betaNTI - Weighted",cex.lab=1.3,main=""); abline(v=c(-2,2),col=2,lwd=2);
# dev.off();

##### Step 4: calculate Bray-Curtis ####
bray = as.matrix(distance(t(match.phylo.otu$data),method = 'bray-curtis'))
write.csv(bray,paste(data.set.name,"_bray_weighted.csv",sep=""),quote=F);

# optional quick check, should be a non-linearily increasing relationship
# identical(colnames(beta.mntd.weighted),colnames(bray.out))
# plot(as.dist(beta.mntd.weighted) ~ as.dist(bray.out))

#### Step 5 Raup-Crick Bray####
metric = 'RC'
abund = T #option to set to false and not have abundance-weighted RC (uncommon)
abund.for.names = 'weighted' # 'unweighted' (uncommon)

#generate null model output
for (i in 1:(no.of.samples - 1)){
  for (j in (i + 1):(no.of.samples)){
    raup.crick.out = raup_crick_abundance_one_comparison(null.one.use = i,null.two.use = j,otu, plot_names_in_col1=FALSE, classic_metric=FALSE, split_ties=TRUE, reps=999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE);
    print(raup.crick.out)
    print(c(i,j))
    write.csv(raup.crick.out,paste('RC_',abund.for.names,'_Diss_comm1_',i,'_comm2_',j,'.csv',sep=""),row.names=T,quote=F);
  }
}

all.RC <- grep("RC_weighted_Diss_comm1_", list.files(pattern = "RC_"), value = TRUE)

#create matrix of comparisons
raup.crick.out = matrix(NA, nrow = nrow(bray), ncol = ncol(bray), dimnames = dimnames(bray))

for (i in 1:length(all.RC)) {
  RC.temp = read.csv(paste(all.RC[i],sep=""),#reading in the file is more complicated because "T" and "F" were being converted to "TRUE" and "FALSE"; this fixes that
                     header = TRUE,check.names = FALSE,
                     stringsAsFactors = FALSE,colClasses = "character")
  row.names(RC.temp) <- RC.temp[[1]];RC.temp <- RC.temp[, -1, drop = FALSE]
  if (nrow(RC.temp) != 1 ) {print(c(i,"ERROR"))} 
  if (ncol(RC.temp) != 1 ) {print(c(i,"ERROR"))} 
  raup.crick.out[which(rownames(raup.crick.out) == rownames(RC.temp)),which(colnames(raup.crick.out) == colnames(RC.temp))] = RC.temp[1,1]
  print(RC.temp)
}

raup.crick.out = as.data.frame(as.matrix(as.dist(raup.crick.out)))
write.csv(raup.crick.out,paste(data.set.name,'_RC_weighted.csv', sep = ''),quote=F);