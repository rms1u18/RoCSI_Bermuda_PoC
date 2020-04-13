### 18S Bermuda pen ocean and coastal heat tree comparisons
# Omiting all blank samples and 2823 (open ocean FF sample that plotted as a coastal sample in nMDS)

# first download packages
# Install phyloseq from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

# Install the rest of the packages from CRAN
install.packages(c("vegan", "metacoder", "taxa", "ggplot2", "dplyr", "readr", "stringr", "agricolae", "ape"),
                 repos = "http://cran.rstudio.com",
                 dependencies = TRUE)


# convert QIIME taxonomy.qza and feature-table.biom to tsv
setwd('OneDrive - University of Southampton/PhD/Bermuda/Coding_folders/Bermuda_18S_QIIME/exported')

# import taxa and OTU data
library(readr)
otu_data <- read_tsv("feature-table.txt")
taxa_data <- read_tsv("taxonomy.tsv")

# Combine taxa and OTU data
library(dplyr) # Loads the dplyr package so we can use `left_join`
taxa_data$`Feature ID` <- as.character(taxa_data$`Feature ID`) # Must be same type for join to work
otu_data$"OTU_ID" <- as.character(otu_data$"OTU_ID") # Must be same type for join to work
otu_data <- left_join(otu_data, taxa_data,
                      by = c("OTU_ID" = "Feature ID")) # identifies cols with shared IDs
print(otu_data)

# Load the sample metadata
setwd("..")
sample_data <- read_tsv("metadata.txt",
                        col_types = "ccccccccc")
print(sample_data)

# Now lets remove samples 2823, 2677a, 2945a, 2626a, 2888a and blanks from the dataset 
# (2823 was an anomoly in the nMDS and the rest are replicate samples with low total reads and were therofore removed so as not to limit the rarfied reads as drastically)
# first from otu_data
otu_data <- select(otu_data, -c("2823-18S","2677a-18S","2945a-18S","2626a-18S","2888b-18S","2500-18S","2570-18S","2629-1-18S","2800-18S","2833-18S","2892-18S","Blank-18S"))# removes 2823,2677a,2945a & blanks
# then from sample_data
sample_data <- sample_data %>% filter(SampleID != "2823-18S") 
sample_data <- sample_data %>% filter(SampleID != "2677a-18S")
sample_data <- sample_data %>% filter(SampleID != "2945a-18S")
sample_data <- sample_data %>% filter(SampleID != "2626a-18S")
sample_data <- sample_data %>% filter(SampleID != "2888b-18S")

#Now lets remove all blanks from sample_data
not_blank <- sample_data$Blank == "No"
sample_data <- sample_data[not_blank,]

# Converting otu_data to taxmap format
library(taxa)
obj <- parse_tax_data(otu_data,
                      class_cols = "Taxon", # The column in the input table
                      class_sep = ";") # What each taxon is seperated by
print(obj)
print(obj$data$tax_data)

# now split out rank information while parsing
obj <- parse_tax_data(otu_data,
                      class_cols = "Taxon",
                      class_sep = ";",
                      class_regex = "^(D_[0-9]{0,1}[0-9]{0,1})_{0,2}(.*)$", # see https://grunwaldlab.github.io/metacoder_documentation/workshop--03--parsing.html for how to use regex
                      class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
print(obj)
head(taxon_ranks(obj))

#So we don’t really need the “class_data” table, so lets get rid of it
obj$data$class_data <- NULL

#Lets also rename the “tax_data” table to something more informative
names(obj$data) <- "otu_counts"
print(obj)

# remove taxon without names
obj <- filter_taxa(obj, taxon_names != "")
print(obj)

head(taxon_names(obj))
head(all_names(obj),20)
length(all_names(obj))

# Let filter any non eukaryota taxa
obj <- filter_taxa(obj, taxon_names == "Eukaryota", subtaxa = TRUE)
print(obj)
filter_taxa(obj, taxon_names == "Eukaryota")

# remove taxa with no reads
has_no_reads <- rowSums(obj$data$otu_counts[, sample_data$SampleID]) == 0
sum(has_no_reads)
filter_obs(obj, "otu_counts", ! has_no_reads) # note the ! negation operator
obj <- filter_obs(obj, "otu_counts", ! has_no_reads, drop_taxa = TRUE)
print(obj)

# now let's plot some heat trees
library(metacoder)
obj %>% 
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "D_3", supertaxa = TRUE) %>% # subset to the order rank
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,
            node_color = n_obs,
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", initial_layout = "reingold-tilford")

### Nows lets try some additional data quality control
# removing low-abundance counts - this sets all read count <10 to 0
library(metacoder)
obj$data$otu_counts <- zero_low_counts(obj, "otu_counts", 
                                       min_count = 10,
                                       cols = c("2508-9-18S","2539-40-18S","2579-80-18S","2608-9-18S","2626b-18S","2677b-18S","2804-5-18S","2840-18S","2872-18S","2888a-18S","2945b-18S"), #all samples
                                       other_cols = TRUE) # keep OTU_ID column
# now filter out all OTUs and associated taxa with a read count of 0
no_reads <- rowSums(obj$data$otu_counts[, sample_data$SampleID]) == 0
sum(no_reads) # when `sum` is used on a TRUE/FALSE vector it counts TRUEs
obj <- filter_obs(obj, "otu_counts", ! no_reads, drop_taxa = TRUE)
print(obj)

# Rarefaction - simulate even number of reads per sample
# Lets look at the distribution of read depths
hist(colSums(obj$data$otu_counts[, sample_data$SampleID]))

obj$data$otu_rarefied <- rarefy_obs(obj,
                                    "otu_counts",
                                    cols = c("2508-9-18S","2539-40-18S","2579-80-18S","2608-9-18S","2626b-18S","2677b-18S","2804-5-18S","2840-18S","2872-18S","2888a-18S","2945b-18S"), #all samples
                                    other_cols = TRUE)
print(obj)

# Some OTUs probab;y now have no read in the rarefied dataset, so we can remove them like before.  However, this time we don't want to drop the taxa
no_reads <- rowSums(obj$data$otu_rarefied[, sample_data$SampleID]) == 0
obj <- filter_obs(obj, "otu_rarefied", ! no_reads)
print(obj)

# Alternatively can calculate the proportion of reads per-sample
obj$data$otu_props <- calc_obs_props(obj,
                                     "otu_counts",
                                     cols = c("2508-9-18S","2539-40-18S","2579-80-18S","2608-9-18S","2626b-18S","2677b-18S","2804-5-18S","2840-18S","2872-18S","2888a-18S","2945b-18S"), #all samples
                                     other_cols = TRUE)
print(obj)

# Comparing taxon abundance in two groups
# first calculate the abundance of each taxon for set of samples from our OTU abundance data
obj$data$tax_abund <- calc_taxon_abund(obj, 
                                      "otu_props",
                                      cols = c("2508-9-18S","2539-40-18S","2579-80-18S","2608-9-18S","2626b-18S","2677b-18S","2804-5-18S","2840-18S","2872-18S","2888a-18S","2945b-18S") #all samples
                                      )

# Then find mean abundance per group using sample characteristics
obj$data$type_abund <- calc_group_mean(obj, "tax_abund",
                                       cols = sample_data$SampleID,
                                       groups = sample_data$Region)
print(obj$data$type_abund)
# now we can use the per-taxon abundances to make a heat tree of the primary taxa present in oligo and coastal
set.seed(2)
obj %>%
  taxa::filter_taxa(Coastal> 0.001) %>% # taxa:: needed because of phyloseq::filter_taxa
  heat_tree(node_label = taxon_names,
            node_size = Coastal, 
            node_color = Coastal, 
            layout = "da", initial_layout = "re", 
            title = "Taxa in Coastal",
            output_file = "Coastal_TREE.pdf")

obj %>%
  taxa::filter_taxa(Oligo> 0.001) %>% # taxa:: needed because of phyloseq::filter_taxa
  heat_tree(node_label = taxon_names,
            node_size = Oligo, 
            node_color = Oligo, 
            layout = "da", initial_layout = "re", 
            title = "Taxa in Open Ocean",
            output_file = "OpenOcean_TREE.pdf")

# Now lets create a differntial heat tree
obj$data$diff_table <- compare_groups(obj, data = "tax_abund",
                                      cols = sample_data$SampleID,
                                      groups = sample_data$Region)
                                      
print(obj$data$diff_table)

# now correct for multiple comparisions and set non-significant difference to zero (false discovery rate detection)
obj <- mutate_obs(obj, "diff_table",
                  wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))



# create a heat tree with a neutral colour in the middle - like grey
set.seed(1)
heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs, # number of OTUs
          node_color = log2_median_ratio, # difference between groups
          node_color_interval = c(-10, 10), # symmetric interval
          node_color_range = c("cyan", "gray", "magenta"), # diverging colors
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Log 2 ratio of median counts")

# plot is difficult to read so let's filter out some more taxa
set.seed(1)
obj %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$")) %>%
  heat_tree(node_label = cleaned_names,
            node_size = n_obs, # number of OTUs
            node_color = log2_median_ratio, # difference between groups
            node_color_interval = c(-10, 10), # symmetric interval
            node_color_range = c("green", "gray", "orange"), # diverging colors
            node_size_axis_label = "OTU count",
            node_color_axis_label = "Log 2 ratio of median counts",
            layout = "da", initial_layout = "re", # good layout for large trees
            title = "Open Ocean vs Coastal samples",
            output_file = "OOvC_HeatTREE.pdf")

# find the log ratio of mean abundances in the two groups - we'll set any non significant differences to 0 and replot to only show significant differences
obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0  # use when comparing sites

set.seed(8)
obj %>%
mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$")) %>%
heat_tree(node_label = cleaned_names,
node_size = n_obs, # number of OTUs
node_color = log2_median_ratio, # difference between groups
node_color_interval = c(-10, 10), # symmetric interval
node_color_range = c("orange", "gray", "lightseagreen"), # diverging colors color blind safe orange coastal, blue open ocean
node_size_axis_label = "OTU count",
node_color_axis_label = "Log 2 ratio of median counts",
layout = "da", initial_layout = "re", # good layout for large trees
title = "18S",
output_file = "18S_OvC_HeatTREEsig.pdf")

# Playing with node labels
set.seed(8)
obj %>%
  mutate_obs("cleaned_names", gsub(taxon_names, pattern = "\\[|\\]", replacement = "")) %>%
  taxa::filter_taxa(grepl(cleaned_names, pattern = "^[a-zA-Z]+$")) %>%
  heat_tree(node_label = cleaned_names,
            node_size = n_obs, # number of OTUs
            node_color = log2_median_ratio, # difference between groups
            node_label_max = 12,
            node_label_size_range = c(0.02, 0.045),
            node_color_interval = c(-10, 10), # symmetric interval
            node_color_range = c("orange", "gray", "lightseagreen"), # diverging colors color blind safe orange coastal, blue open ocean
            node_size_axis_label = "OTU count",
            node_color_axis_label = "Log 2 ratio of median counts",
            layout = "da", initial_layout = "re", # good layout for large trees
            output_file = "18S_OvC_HeatTREEsig_pub.pdf")
