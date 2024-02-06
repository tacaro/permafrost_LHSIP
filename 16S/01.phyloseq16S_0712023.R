# Sierra Jech
# Barger Lab
# July 2023
# Using ASV Table from Dada2 on microbe on 07102023

# Purpose: Preliminary EDA for permafrost samples 

# Load libraries
library(tidyverse)
library(phyloseq) # For all bioinformatics analyses post Dada2
packageVersion('phyloseq')
library(ggplot2)
#library(vegan)
#library(car) # For Levene Test and anova
#library(PMCMRplus) # For Nemenyi posthoc test
#library("DESeq2")
#library(indicspecies) # For indicator species
#library(microbiome)

# Load data
### Load the ASV table (called OTU in phyloseq)
list.files()
otumat <- read.table(file = "data/seqtab_wTax_mctoolsr.txt", header=T)
dim(otumat) # 44 samples and 15,987 ASVs
otumat <- otumat %>%
  tibble::column_to_rownames("X.ASV_ID") # push the ASV_ID column into the rownames position

### Load the taxonomy table - this comes from two columns in the dada2 output
# Find the ASV names
asv <- rownames(otumat)
# Find the taxonomy info
tax <- otumat$taxonomy
# Combine them into a dataframe
taxmat <- as.data.frame(cbind(asv, tax)) 
# change the column names to be "ASV_ID" and taxonomy
colnames(taxmat) <- c("ASV_ID", "taxonomy")
# split the column into each taxonomy
taxmat <- separate(taxmat, 
                   col = taxonomy, 
                   into= c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

# Format tax_sep for phyloseq by making ASV_ID the rownames
taxmat <- taxmat %>%
  tibble::column_to_rownames("ASV_ID")
# Delete the species column because it is ASV_ID (it is in column 7)
taxmat <- taxmat[, 1:6]
# remove the taxonomy column from the OTU table (it is column 385)
otumat <- otumat[, -43]

#### Load the data mapping (metadata) file
#list.files()
mapping <- read.table("data/permafrost16S_mapping.txt", header = TRUE) # all metadata
colnames(mapping)
mapping <- remove_rownames(mapping)
sampledata <- column_to_rownames(mapping, var = "SampleID") # rownames should be the sample IDs

# Transform to phyloseq objects
# check that the names match
all.equal(sort(colnames(otumat)), sort(rownames(sampledata))) 

# they do not match because the negative sign for "-4C" was changed to a decimal for some of the sample names
row.names(sampledata)
colnames(otumat)
colnames(otumat)[which(names(otumat) == "p35m..4C.0d")] <- "p35m.-4C.0d"
colnames(otumat)[which(names(otumat) == "p35m..4C.180d")] <- "p35m.-4C.180d"
colnames(otumat)[which(names(otumat) == "p35m..4C.30d")] <- "p35m.-4C.30d"
colnames(otumat)[which(names(otumat) == "p35m..4C.7d")] <- "p35m.-4C.7d"
colnames(otumat)[which(names(otumat) == "p54m..4C.0d")] <- "p54m.-4C.0d"
colnames(otumat)[which(names(otumat) == "p54m..4C.180d")] <- "p54m.-4C.180d"
colnames(otumat)[which(names(otumat) == "p54m..4C.30d")] <- "p54m.-4C.30d"
colnames(otumat)[which(names(otumat) == "p54m..4C.7d")] <- "p54m.-4C.7d"
colnames(otumat)[which(names(otumat) == "p83m..4C.0d")] <- "p83m.-4C.0d"
colnames(otumat)[which(names(otumat) == "p83m..4C.180d")] <- "p83m.-4C.180d"
colnames(otumat)[which(names(otumat) == "p83m..4C.30d")] <- "p83m.-4C.30d"
colnames(otumat)[which(names(otumat) == "p83m..4C.7d")] <- "p83m.-4C.7d"
colnames(otumat)

# check that the names match
all.equal(sort(colnames(otumat)), sort(rownames(sampledata))) 

# Transform otu table and taxonomy tables into matrices
otumat <- as.matrix(otumat)
taxmat <- as.matrix(taxmat)
#View(otu_mat)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
samples = sample_data(sampledata)
permafrost_16S_raw.ps <- phyloseq(OTU, TAX, samples)
permafrost_16S_raw.ps
# checks out

#### Check on the naming
colnames(otu_table(permafrost_16S_raw.ps)) # names associated with the otu table
sample_names(permafrost_16S_raw.ps) #names associated with the mapping table
sample_variables(permafrost_16S_raw.ps)
sample_sums(permafrost_16S_raw.ps) #returns the ASV count for each sample
# Blanks have 20 and 36 reads
# samples have min = 94,743, max = 248,921

#### Get Phylum counts ####
rank_names(permafrost_16S_raw.ps)
table(tax_table(permafrost_16S_raw.ps)[, "Phylum"], exclude = NULL)
# probably want to remove the NA phyla and the phyla with only 1 feature
permafrost_16S_removeNA.ps <- subset_taxa(permafrost_16S_raw.ps, !is.na(Phylum)) # I am not sure that this works


#### Prevalence ####
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(permafrost_16S_removeNA.ps),
               MARGIN = ifelse(taxa_are_rows(permafrost_16S_removeNA.ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(permafrost_16S_removeNA.ps),
                    tax_table(permafrost_16S_removeNA.ps))
head(arrange(prevdf,-TotalAbundance))
head(arrange(prevdf,-Prevalence)) # only a few ASVs can be found in every sample!
#prevdf
# So the prevalence column gives me a count of the number of samples which have that ASV. We do not really have a reason to limit the prevalence of ASVs in the analysis...? 
# note, the Total Abundance column is "the sum of all reads observed for each ASV". When people say that they removed rare species which had < 1% abundance...they are talking about an abundance cutoff. The phyloseq manual uses a prevalence cutoff of 5%. 
# Note that there are many ASVs present in the list which have a prevalence of 1

# CONTINUE PREVAENCE CODE HERE IF DESIRED
prevalenceThreshold = 0.05 * nsamples(permafrost_16S_removeNA.ps)
prevalenceThreshold # so they have to be in at least 2.1 samples to be included in the analysis

#Are there phyla that are comprised of mostly low-prevalence features? Compute the total and average prevalence of the features in each phylum.
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# now prune those with a prevalence < 5%
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
length(keepTaxa) # 2052 ASVs to keep
permafrost_16S_removeNA_prevalent.ps = prune_taxa(keepTaxa, permafrost_16S_removeNA.ps)

# check your work
table(tax_table(permafrost_16S_removeNA_prevalent.ps)[, "Phylum"], exclude = NULL)
# this table is smaller, I think it tells us the number of ASVs in each phylum


##################################################################
# FILTER OUT MITOCHONDRIA, CHLOROPLASTS, and KINGDOM EUKARYOTA
# chloroplast is an order and mitochondria is a family in this version of SILVA
##################################################################
# Get taxonomy table out of phyloseq:
permafrost16S_taxTable <- psmelt(permafrost_16S_removeNA_prevalent.ps) # this takes some time, be patient

# First check what will be removed
tax_filt <- permafrost16S_taxTable %>%
  filter(Order != "Chloroplast") %>%
  filter(Family != "Mitochondria") %>%
  filter(Kingdom != "Eukaryota") %>%
  filter(Kingdom != "NA") %>% #Remove ASVs where kingdom is unknown ("NA" in first column)
  filter(Phylum != "NA") #Remove ASVs where phylum is unknown ("NA" in second column) # many will be removed!

# How many were removed?
dim(permafrost16S_taxTable)[1] - dim(tax_filt)[1]
# We will filter out a total of 490,992 instances of ASVs
1302/86184 # We are removing 1.5% of ASV instances
84882/86184 # We are keeping 98.5% of ASV instances


# Now actually remove them from the pyloseq object. First find and save each set as their own object. We will then remove all the ones we do not want below. 
# Chloroplasts (order)
chloros <- permafrost16S_taxTable %>%
  filter(Order == "Chloroplast")
dim(chloros) #168 chloroplast ASVs
chloro_names <- chloros$OTU
chloro_names <- unique(chloro_names)
length(chloro_names) #4 unique chloroplast ASVs are being removed 

# Mitochondria (family)
mitos <- permafrost16S_taxTable %>%
  filter(Family == "Mitochondria")
dim(mitos) #378 mitochondria ASVs
mito_names <- mitos$OTU
length(mito_names)
mito_names <- unique(mito_names)
length(mito_names) #9 unique chloroplast ASVs are being removed 

# Eukaryota (kingdom)
kingdomEuks <- permafrost16S_taxTable %>%
  filter(Kingdom == "Eukaryota")
dim(kingdomEuks) #0
euks_names <- rownames(kingdomEuks)
euks_names

# NA's (kingdom) - these were removed above as NA's
kingdomNAs <- permafrost16S_taxTable %>%
   filter(Kingdom == "NA")
dim(kingdomNAs) #0

# NA's (phyla)
PhylumNAs <- permafrost16S_taxTable %>%
  filter(Phylum == "NA")
dim(PhylumNAs) #756
pNAs_names <- PhylumNAs$OTU
pNAs_names <- unique(pNAs_names)
length(pNAs_names) #18

# join all ASV IDs that should be removed into one list
remove_ASVs <- c(chloro_names, mito_names, euks_names, pNAs_names)
length(remove_ASVs) # removing thousands instances of ASVs, SDJ updated this to be just the unique ones which is 31 ASVs removed
#check it with real math
4+9+18 # 31 instances of ASV's
# are there duplicates?
temp <- remove_ASVs[!duplicated(remove_ASVs)]
length(temp) # there are only 31 ASVs that we are actually removing

# Removing ASVs in phyloseq comes from: Joey711's code: https://github.com/joey711/phyloseq/issues/652
# Remove the ASVs that identified:
all_Taxa <- taxa_names(permafrost_16S_removeNA_prevalent.ps) #get all tax names in original, uncleaned dataset
ASVstoKeep <- all_Taxa[!(all_Taxa %in% temp)]
length(ASVstoKeep) #2,021
permafrost_16S_ps <- prune_taxa(ASVstoKeep, permafrost_16S_removeNA_prevalent.ps) # new phyloseq object with just the taxa we want!

#### Get Phylum counts again ####
rank_names(permafrost_16S_ps)
table(tax_table(permafrost_16S_ps)[, "Phylum"], exclude = NULL)
#### Prevalence ####
# Compute prevalence of each feature, store as data.frame
prevdf_filt = apply(X = otu_table(permafrost_16S_ps),
                    MARGIN = ifelse(taxa_are_rows(permafrost_16S_ps), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_filt = data.frame(Prevalence = prevdf_filt,
                         TotalAbundance = taxa_sums(permafrost_16S_ps),
                         tax_table(permafrost_16S_ps))
head(arrange(prevdf_filt,-TotalAbundance))
# there should be no chloroplast, mitochondria reads in the dataset anymore. We should also not have any ASVs without reads or with NA kingdom or phyla assignments


##################################################################
# Checking the data
# Code here comes from Cliff Tutorial
##################################################################
# again calculate the number of ASV per sample
sort(phyloseq::sample_sums(permafrost_16S_ps)) 
# minimum reads in blanks = 20
# maximum reads in blanks = 28
# minimum reads in a sample = 68455
# maximum reads in a sample = 235349
mean(phyloseq::sample_sums(permafrost_16S_ps))
# mean reads per sample = 125409.9
median(phyloseq::sample_sums(permafrost_16S_ps))
# median reads per sample = 126256.5

# mean without the blanks
mean(sort(phyloseq::sample_sums(permafrost_16S_ps))[3:42])
# 131679.1
# median without the blanks
median(sort(phyloseq::sample_sums(permafrost_16S_ps))[3:42])
#128529.5

# How are number of reads distributed across the samples?
seqcounts <- as.data.frame(sort(colSums(otu_table(permafrost_16S_ps)))) %>%
  rename("seqs" = "sort(colSums(otu_table(permafrost_16S_ps)))") %>%
  rownames_to_column(var = "sampleID")
head(seqcounts, 10)
tail(seqcounts,10)

# Now we have a dataframe with two columns, seqs and sampleID which we can plot
ggplot(seqcounts, aes(reorder(sampleID, seqs, mean), seqs)) + # Dataframe and variables
  geom_bar(stat = "identity") + # Type of graph
  labs(y = "# Reads", x = "Sample") + # Axes labels
  coord_flip() + # Flip axes
  geom_hline(yintercept = 68450, color = "blue") + # this seems backwards because of the coord_flip() command above
  theme_classic() + 
  theme(axis.text.y = element_text(size = 2)) # this might help you decide to rarify
colSums(otu_table(permafrost_16S_ps)) # all read counts are different
sum(colSums(otu_table(permafrost_16S_ps))) # total number of reads in the dataset is 5,267,214

# visualization to make sure things look right...these are too large and need to be subset before running or R will crash
#plot_bar(mayberrycyano_16S_ps, fill = "Phylum")
# Alternatively
#ps.phylum = tax_glom(pits_16S_filt_ps, taxrank="Phylum", NArm=FALSE)
#plot_bar(ps.phylum, fill="Phylum")

####################################################################
# Save these Pre-Processing Steps so that you do not have to do them again
save(permafrost_16S_ps, file = "permafrost16S_processed_filt_ps.RData")
saveRDS(permafrost_16S_ps, file = "permafrost16S_processed_filt_ps.rds")
####################################################################


#### Check which ASVs that are in the blanks ###
# subset for the blank samples
samplesToKeep <- c("Blank1","Blank2")
blanks <- prune_samples(sample_names(permafrost_16S_ps) %in% samplesToKeep, permafrost_16S_ps)
colSums(otu_table(blanks)) # number of reads in these samples
table(tax_table(blanks)[, "Phylum"], exclude = NULL) # 35 phyla shown???
melt_blanks <- psmelt(blanks)
# plot it 
blanks.phylum = tax_glom(blanks, taxrank="Phylum", NArm=FALSE)
plot_bar(blanks.phylum, fill="Phylum")
blanks.phylum@otu_table
# I see 20 reads for ASV_2416 in Blank 1, 26 reads for ASV_2416 in Blank2, and 2 reads for ASV_221 in Blank 2
blanks.phylum@tax_table['ASV_2416', ]
# ASV_2416 = Deinococcota
blanks.phylum@tax_table['ASV_221', ]
# ASV_221 = Chloroflexi

# so there is some consistency in the blanks, we might consider removing 20-25 reads from each sample for ASV_2416 if it becomes important.

