# Alpha and Beta Diversity Analysis
# Permafrost Project

# Follows after running code in '01.phyloseq16S_07122023.R' and obtaining the processed phyloseq object

# The purpose of this code is to:

### Libraries ###
library(phyloseq) # For all bioinformatics analyses post Dada2
packageVersion('phyloseq')
library(ggplot2)
library(tidyverse)
library(vegan)
library(car) # For Levene Test and anova
#library(PMCMRplus) # For Nemenyi posthoc test
library(FSA) # dunnTest for post-hoc
#library("DESeq2")
#library(indicspecies) # For indicator species MUTIPATT
#library(microbiome)
#library(ggpubr) #for ggarrange
#library(stringr)
library(RVAideMemoire) # for pairwise.perm.manova
#library(scales)


### Functions ###
se <- function(x) {sd(x,na.rm=TRUE)/sqrt(length(x))}

### Load Data ###
permafrost_16S_ps <- readRDS("permafrost16S_processed_filt_ps.rds")

#### Get to know the data ####
permafrost_16S_ps@sam_data@names #colnames for the sample_data information in this phyloseq object
permafrost_16S_ps@sam_data[['names']] <- row.names(permafrost_16S_ps@sam_data) # view all samples

# Subsetting by sample name:
samplesToRemove <- c("Blank1","Blank2")

# now actually remove them in the ps object
permafrost_16S_ps_noblanks <- subset_samples(permafrost_16S_ps, !names %in% samplesToRemove)

# factor the label levels for plotting
permafrost_16S_ps_noblanks@sam_data$depth <- factor(permafrost_16S_ps_noblanks@sam_data$depth, levels = c("US", "35", "54", "83"))
permafrost_16S_ps_noblanks@sam_data$incubation_temp_C <- factor(permafrost_16S_ps_noblanks@sam_data$incubation_temp_C, levels = c("-4", "4", "12"))
#
# permafrost_16S_ps_noblanks@sam_data$incubation_time_d # this one should be integer


### Alpha Diversity ###

# Is richness different?
pAlpha = plot_richness(permafrost_16S_ps_noblanks,
                       x = "incubation_time_d", 
                       measures = c("Observed"), 
                       title = "Alpha Diversity") 
#color = "comboLabel")
# can add others here like ("Shannon,"Chao1", "InvSimpson")

pAlpha + 
  geom_point(aes(shape = depth), size = 2) + 
  theme_classic() + 
  labs(x = "Incubation Time (days)", 
       y = "Alpha Diversity") + 
  #scale_color_manual(values = init_biocrust_colors2) + 
  theme(axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12), 
        axis.title.y = element_text(size=14), 
        axis.title.x = element_text(size=14))+
  facet_wrap(~incubation_temp_C, scales = "free_x")



# estimate richness
rich <- estimate_richness(permafrost_16S_ps_noblanks)




# Ordination
bray_dist <- phyloseq::distance(permafrost_16S_ps_noblanks, method="bray", weighted= F)
ordination <- ordinate(permafrost_16S_ps_noblanks, method="PCoA", distance=bray_dist)
#ordination.nmds <- ordinate(year2_cp, method="NMDS", distance="bray") # insufficient data here

permafrost_16S_ps_noblanks@sam_data$incubation_time_d <- factor(permafrost_16S_ps_noblanks@sam_data$incubation_time_d, levels = c("0", "7", "30", "180"))

#Visualization
plot_ordination(permafrost_16S_ps_noblanks, ordination, shape = "depth", color = "incubation_time_d") +
  theme(aspect.ratio=1) + 
  geom_point(size = 3.5) +
  theme_classic()+
  scale_color_manual("Incubation Time (days)", values = c("grey80", "grey59", "grey23","black"))+
  scale_shape_manual("Depth (m)", values = c(3,19,17, 15), labels = c("US","35", "54", "83"))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face="bold"),
        axis.title.y = element_text(size = 14, face="bold"),
        strip.text.x = element_text(size = 12))+
  facet_wrap(~incubation_temp_C, scales = "free")

#ggsave("output/PCoA.pdf", width = 10, height = 3)



# PERMANOVA - test for differences in centroid of different pCoA hull groups
permafrost_16S_ps_noblanks@sam_data$incub_temptime <- paste(permafrost_16S_ps_noblanks@sam_data$incubation_temp_C, permafrost_16S_ps_noblanks@sam_data$incubation_time_d, sep = ".")
# are there differences in centroids based on each of the parameters?

# Depth
adon_depth <- adonis2(bray_dist ~ sample_data(permafrost_16S_ps_noblanks)$depth, method = "bray") # centroid locations significantly different for the groups
adon_depth # significant 
pairwise.perm.manova(bray_dist, sample_data(permafrost_16S_ps_noblanks)$depth, test = "Wilks", p.method = "holm", F = TRUE, R2 = TRUE)
# yes all of these are significant

# Incubation Temperature
adon_temp <- adonis2(bray_dist ~ sample_data(permafrost_16S_ps_noblanks)$incubation_temp_C, method = "bray") # centroid locations significantly different for the groups
adon_temp # not significant 

# Incubation Time
adon_time <- adonis2(bray_dist ~ sample_data(permafrost_16S_ps_noblanks)$incubation_time_d, method = "bray") # centroid locations significantly different for the groups
adon_time # not significant 
 


# PERMDISP - multivariate version of Levenne Test. Difference in dispersion for each factor level
m1 <- betadisper(bray_dist, sample_data(permafrost_16S_ps_noblanks)$depth)
anova(m1) # Dispersion is not different by depth
# if permdisp is significant, then you would want to use a different test than PERMANOVA to compare the groups
#TukeyHSD(m1) # post hoc test
#scores(m1) # access the scores


#### Plotting Phyla Relative Abundance ####
# dataframe of the OTUs for later analyses
OTU <- as.data.frame(permafrost_16S_ps_noblanks@otu_table)

#Get count of reads by phyla
table(phyloseq::tax_table(permafrost_16S_ps_noblanks)[, "Phylum"]) #35 Phyla included
# get the total read count per sample
reads_per_sample <- as.data.frame (sample_sums(permafrost_16S_ps_noblanks)) # they are all different because we did not rarefy in this dataset
#copy the rownames to a column
reads_per_sample <- tibble::rownames_to_column(reads_per_sample, "Sample")
# change the name of the second column
colnames(reads_per_sample) <- c("Sample", "totalReads")
nrow(reads_per_sample) # 40 samples are included

# Think about the phyla 
#get_taxa_unique(initial_biocrust, "Phylum") # 34 unique Phyla
# I only want to plot the top 9 and the rest go into an Other category

# first melt the data out of phyloseq format and into long dataframe
mb_psmelt <- psmelt(permafrost_16S_ps_noblanks)

# Get a list of the most abundant Phyla by depth
phyla_summary <- mb_psmelt %>%
  group_by(depth, Phylum) %>%
  summarise(totalAbund = sum(Abundance)) %>%
  arrange(-totalAbund)

# Get a list of the most abundant Phyla
temp <- mb_psmelt %>%
  group_by(Phylum) %>%
  summarise(totalAbund = sum(Abundance)) %>%
  arrange(-totalAbund)

temp <- temp[1:9,]

top9Phyla <- temp$Phylum # list that holds the Top9Phyla names
top9Phyla

#subset the psmelt dataframe for only the Top9Phyla
phylaSubset <- mb_psmelt %>%
  filter(Phylum %in% top9Phyla)

# check that we got the right ones
unique(phylaSubset$Phylum) # yes

# calculate the reads per Phylum, grouped by Sample
phylaSubset_abund <- phylaSubset %>%
  group_by(Sample, Phylum) %>%
  summarise(abund = sum(Abundance))

#make a new column were I divide each of the values by the correct one in the reads_per_sample list
# first merge the reads_per_sample data onto the dataframe
phylaSubset_abund2 <- inner_join(phylaSubset_abund, reads_per_sample, by = "Sample")
# divide the abund column by the total reads column
phylaSubset_abund2$relAbund <- phylaSubset_abund2$abund / phylaSubset_abund2$totalReads
# check that it does not add to 1 because there should be some in the Other category
total_relAbund <- phylaSubset_abund2 %>%
  group_by(Sample) %>%
  summarise(sum = sum(relAbund))
# add a new column which contains the relative abundance value for the Other category
total_relAbund$other <- (1 - total_relAbund$sum)
# delete the sum column
total_relAbund <- total_relAbund %>%
  select(Sample, other)
#add column which is "Other" repeated
total_relAbund$Phylum <- "Other"
#change column header "other" to relAbund
colnames(total_relAbund)[which(names(total_relAbund) == "other")] <- "relAbund"
# select columns to keep in the dataframe we want
phylaSubset_abund2 <- phylaSubset_abund2 %>%
  select(Sample, relAbund, Phylum)
# rbind the other values to the phylaSubset_abund2 dataframe
phylaSubset_abund3 <- rbind(phylaSubset_abund2, total_relAbund)

# Now check that they sum to 1
total_relAbund2 <- phylaSubset_abund3 %>%
  group_by(Sample) %>%
  summarise(sum = sum(relAbund)) # YES!!!!
head(total_relAbund2)

# plot it!
# Plot Relative Abundance
temp2 <- unique(mb_psmelt %>% select(Sample, depth, incubation_temp_C, incubation_time_d))
phylaSubset_abund3 <- left_join(phylaSubset_abund3, temp2, by = "Sample")

phylaSubset_abund3$Sample <- factor(phylaSubset_abund3$Sample, levels = c("pUS.12C.0d", "pUS.12C.7d", "pUS.12C.30d", "pUS.12C.180d", "p35m.-4C.0d", "p35m.-4C.7d", "p35m.-4C.30d", "p35m.-4C.180d", "p54m.-4C.0d", "p54m.-4C.7d", "p54m.-4C.30d", "p54m.-4C.180d", "p83m.-4C.0d", "p83m.-4C.7d", "p83m.-4C.30d", "p83m.-4C.180d", "p35m.4C.0d", "p35m.4C.7d", "p35m.4C.30d", "p35m.4C.180d", "p54m.4C.0d", "p54m.4C.7d", "p54m.4C.30d", "p54m.4C.180d", "p83m.4C.0d", "p83m.4C.7d", "p83m.4C.30d", "p83m.4C.180d", "p35m.12C.0d", "p35m.12C.7d", "p35m.12C.30d", "p35m.12C.180d", "p54m.12C.0d", "p54m.12C.7d",  "p54m.12C.30d", "p54m.12C.180d", "p83m.12C.0d", "p83m.12C.7d", "p83m.12C.30d", "p83m.12C.180d"))

  

# Plot Relative Abundance - basic
phylaSubset_abund3 %>%
  #filter(depth == "US") %>%
  #filter(depth == "35") %>%
  #filter(depth == "54") %>%
  filter(depth == "83") %>%
  ggplot(aes(y = relAbund, x = incubation_time_d, fill = Phylum)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.25) +
  labs(x = "Incubation Time (days)", y = "Relative Abundance", fill = "Phyla") +
  ggtitle(label = "Depth: US")+
  theme_classic()+
  coord_cartesian(ylim = c(0,1), expand = FALSE)+
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(face = "italic"),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 12))+
  facet_grid(cols = vars(incubation_temp_C), scales = "free")

#ggsave("output/phyla_stackedbar_US.pdf", width = 4, height = 4)





# Archaea
# Crenarchaeota, Euryarchaeota, Korarchaeota, Nanoarchaeota, and Thaumarchaeota
table(phyloseq::tax_table(permafrost_16S_ps_noblanks)[, "Phylum"])
# Crenarchaeota = 10
# Euryarchaeota = 2
# Nanoarchaeota = 2
# Thaumarchaeota = 0

taxa_table <- as.data.frame(tax_table(permafrost_16S_ps_noblanks))

# filter the OTU table for these archaea and see how abundant they were
archaea <- taxa_table %>% filter(Phylum %in% c("Crenarchaeota", "Euryarchaeota","Nanoarchaeota", "Thaumarchaeota"))

# create a list of ASVs which are archaea
archaeaASVs <- row.names(archaea)

# subset the mb_psmelt object for these ASVs
archaeaOTU <- mb_psmelt %>% filter(OTU %in% archaeaASVs)

# Plot all ASVs 
# make a OTU label
archaeaOTU$taxa <- paste(archaeaOTU$Phylum,";", archaeaOTU$Class,";",archaeaOTU$Order, ";",archaeaOTU$Family, ";", archaeaOTU$Genus, ";", archaeaOTU$OTU)

archaeaOTU %>%
  ggplot()+
  geom_line(mapping = aes(x = incubation_time_d, y = relAbund, group = OTU, color = taxa)) +
  geom_point(mapping = aes(x = incubation_time_d, y = relAbund, group = OTU, color = taxa)) +
  theme_classic()+
  labs(x = "Incubation Time (days)", y = "Relative Abundance")+
  ggtitle("All Archaea")+
  facet_grid(rows = vars(depth), cols = vars(incubation_temp_C))






# INDICATORS OF THAW
# Get list of ASVs that are abundant at higher temperatures and at longer time scales
# Plot mean abundance of individual ASVs over time, color the lines 
# Add the reads_per_sample to the dataframe for calculating the ASV-level relative abundances
mb_psmelt <- left_join(mb_psmelt, reads_per_sample, by = "Sample")
# Calculate the relative abundance
mb_psmelt$relAbund <- mb_psmelt$Abundance / mb_psmelt$totalReads
# Pick the top 10 ASVs
top10ASV <- mb_psmelt %>% group_by(OTU, depth) %>% summarise(sumRelAbund = sum(relAbund))
top10ASV_US <- top10ASV %>% filter(depth == "US") %>% arrange(desc(sumRelAbund))
top10ASV_35 <- top10ASV %>% filter(depth == 35) %>% arrange(desc(sumRelAbund))
top10ASV_54 <- top10ASV %>% filter(depth == 54) %>% arrange(desc(sumRelAbund))
top10ASV_83 <- top10ASV %>% filter(depth == 83) %>% arrange(desc(sumRelAbund))

# if I take the top 10 in each category, what proportion of the reads am I capturing?
sum(top10ASV_US$sumRelAbund[1:10])
sum(top10ASV_35$sumRelAbund[1:10])
sum(top10ASV_54$sumRelAbund[1:10])
sum(top10ASV_83$sumRelAbund[1:10])

# get lists of the OTU names for each depth
US_OTUs <- top10ASV_US$OTU[1:10]
D35_OTUs <- top10ASV_35$OTU[1:10]
D54_OTUs <- top10ASV_54$OTU[1:10]
D83_OTUs <- top10ASV_83$OTU[1:10]

# subset the main dataframe for these OTUs
US_top10ASVs <- mb_psmelt %>% filter(depth == "US") %>% filter(OTU %in% US_OTUs) 
# there are not replicates, simply plot the values 

# make a OTU label
US_top10ASVs$taxa <- paste(US_top10ASVs$Phylum,";", US_top10ASVs$Class,";",US_top10ASVs$Order, ";",US_top10ASVs$Family, ";", US_top10ASVs$Genus, ";", US_top10ASVs$OTU)

US_top10ASVs %>%
  ggplot()+
  geom_line(mapping = aes(x = incubation_time_d, y = relAbund, group = OTU, color = taxa)) +
  geom_point(mapping = aes(x = incubation_time_d, y = relAbund, group = OTU, color = taxa)) +
  theme_classic()+
  labs(x = "Incubation Time (days)", y = "Relative Abundance")+
  ggtitle("Surface Soils")



# subset the main dataframe for these OTUs
D35_top10ASVs <- mb_psmelt %>% filter(depth == "35") %>% filter(OTU %in% US_OTUs) 
# there are not replicates, simply plot the values 

# make a OTU label
D35_top10ASVs$taxa <- paste(D35_top10ASVs$Phylum,";", D35_top10ASVs$Class,";",D35_top10ASVs$Order, ";",D35_top10ASVs$Family, ";", D35_top10ASVs$Genus, ";", D35_top10ASVs$OTU)

D35_top10ASVs %>%
  ggplot()+
  geom_line(mapping = aes(x = incubation_time_d, y = relAbund, group = OTU, color = taxa)) +
  geom_point(mapping = aes(x = incubation_time_d, y = relAbund, group = OTU, color = taxa)) +
  theme_classic()+
  labs(x = "Incubation Time (days)", y = "Relative Abundance")+
  ggtitle("Depth: 35 m")+
  facet_grid(cols = vars(incubation_temp_C))




# subset the main dataframe for these OTUs
D54_top10ASVs <- mb_psmelt %>% filter(depth == "54") %>% filter(OTU %in% US_OTUs) 
# there are not replicates, simply plot the values 

# make a OTU label
D54_top10ASVs$taxa <- paste(D54_top10ASVs$Phylum,";", D54_top10ASVs$Class,";",D54_top10ASVs$Order, ";",D54_top10ASVs$Family, ";", D54_top10ASVs$Genus, ";", D54_top10ASVs$OTU)

D54_top10ASVs %>%
  ggplot()+
  geom_line(mapping = aes(x = incubation_time_d, y = relAbund, group = OTU, color = taxa)) +
  geom_point(mapping = aes(x = incubation_time_d, y = relAbund, group = OTU, color = taxa)) +
  theme_classic()+
  labs(x = "Incubation Time (days)", y = "Relative Abundance")+
  ggtitle("Depth: 54 m")+
  facet_grid(cols = vars(incubation_temp_C))



# subset the main dataframe for these OTUs
D83_top10ASVs <- mb_psmelt %>% filter(depth == "83") %>% filter(OTU %in% US_OTUs) 
# there are not replicates, simply plot the values 

# make a OTU label
D83_top10ASVs$taxa <- paste(D83_top10ASVs$Phylum,";", D83_top10ASVs$Class,";",D83_top10ASVs$Order, ";",D83_top10ASVs$Family, ";", D83_top10ASVs$Genus, ";", D83_top10ASVs$OTU)

D83_top10ASVs %>%
  ggplot()+
  geom_line(mapping = aes(x = incubation_time_d, y = relAbund, group = OTU, color = taxa)) +
  geom_point(mapping = aes(x = incubation_time_d, y = relAbund, group = OTU, color = taxa)) +
  theme_classic()+
  labs(x = "Incubation Time (days)", y = "Relative Abundance")+
  ggtitle("Depth: 83 m")+
  facet_grid(cols = vars(incubation_temp_C))
