#extracting blank sequences using microDecon
#before starting open feature-table.tsv in excell delete first line and save as .xlsx

# set working directory
setwd("/Users/rms1u18/Desktop/Bermuda_16S_QIIME/exported/")

install.packages("readxl")
library ("readxl")
install.packages("devtools")
devtools::install_github("donaldtmcknight/microDecon")
library(microDecon)
install.packages("tidyverse")
library(tidyverse)

# read in table
OTUtable <- read_excel("feature-table.xlsx")

#split table by blank grouping into data frames$
# starting with JR method extraction blank
JRBlank <- as.data.frame(OTUtable[c("OTU_ID", "Blank-16S", "2626a-16S", "2626b-16S", "2677a-16S", "2677b-16S", "2888a-16S", "2888b-16S", "2945a-16S", "2945b-16S")])

# use decon function to decontaminate the data
JRDecon <- decon(data = JRBlank, numb.blanks = 1, numb.ind = 8, taxa = F)

Blank2800 <- as.data.frame(OTUtable[c("OTU_ID", "2800-16S", "2804-5-16S", "2823-16S")])
Decon2800 <- decon(data = Blank2800, numb.blanks = 1, numb.ind = 2, taxa = F)

Blank2500 <- as.data.frame(OTUtable[c("OTU_ID", "2500-16S", "2508-9-16S", "2539-40-16S")])
Decon2500 <- decon(data = Blank2500, numb.blanks = 1, numb.ind = 2, taxa = F)

Blank2833 <- as.data.frame(OTUtable[c("OTU_ID", "2833-16S", "2840-16S", "2872-16S")])
Decon2833 <- decon(data = Blank2833, numb.blanks = 1, numb.ind = 2, taxa = F)

Blank2570 <- as.data.frame(OTUtable[c("OTU_ID", "2570-16S", "2579-80-16S", "2608-9-16S")])
Decon2570 <- decon(data = Blank2570, numb.blanks = 1, numb.ind = 2, taxa = F)

# reassemble the decontaminated OTU in to one table
require(tidyverse)
lst <- lst(JRDecon$decon.table, Decon2500$decon.table, Decon2570$decon.table, Decon2800$decon.table, Decon2833$decon.table)
DeconOTUs <- reduce(lst, full_join, by = "OTU_ID") %>% replace(., is.na(.), 0);

# export table to tab delimited .txt
write.table(DeconOTUs, "/Users/rms1u18/Desktop/Bermuda_16S_QIIME/R/DeconOTUs.txt", quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)

# convert back to .biom in terminal use Export_import_feature_table pipe
