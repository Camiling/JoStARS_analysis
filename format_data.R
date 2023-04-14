rm(list=ls())



# Load TCGA BRCA, UCEC and OVAC data from USC XENABROWSER 



# Load BRCA data --------------------------------------------------------------------------

# We use RPPA (replicate-base normalization) data
# Data from https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FRPPA_RBN&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
options(timeout=200)
utils::download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FRPPA_RBN.gz", "data/TCGA.BRCA.sampleMap%2FRPPA_RBN.gz")
BRCA_download <- utils::read.table(gzfile("data/TCGA.BRCA.sampleMap%2FRPPA_RBN.gz"), header = T)
rownames(BRCA_download) = BRCA_download[,1]
BRCA_download = t(BRCA_download[,-1])
dim(BRCA_download) # 747 x 131
# Controlled that all Ids end with '01', i.e., no controls. 

# Get PAM50 phenotypes and select only luminal A and B patients
utils::download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix", "data/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix")
BRCA_Pam50_download <- utils::read.csv("data/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix", sep='\t')
rownames(BRCA_Pam50_download) = BRCA_Pam50_download[,1]
BRCA_Pam50_download = BRCA_Pam50_download[,-1]
IDs.lumA = rownames(BRCA_Pam50_download)[which(BRCA_Pam50_download$PAM50_mRNA_nature2012=='Luminal A')]
IDs.lumB = rownames(BRCA_Pam50_download)[which(BRCA_Pam50_download$PAM50_mRNA_nature2012=='Luminal B')]
IDs.lumA = gsub("-", ".",IDs.lumA)
IDs.lumB = gsub("-", ".",IDs.lumB)

IDs.lumA.rppa = rownames(BRCA_download)[which(rownames(BRCA_download) %in% IDs.lumA)]
IDs.lumB.rppa = rownames(BRCA_download)[which(rownames(BRCA_download) %in% IDs.lumB)]
length(IDs.lumA.rppa) # 173
length(IDs.lumB.rppa) # 100

brca_dat = BRCA_download[which(rownames(BRCA_download) %in% c(IDs.lumA.rppa, IDs.lumB.rppa)),]
dim(brca_dat) # 273 x 131

# Load UCEC data --------------------------------------------------------------------------

# We use RPPA (replicate-base normalization) data
# Data from https://xenabrowser.net/datapages/?dataset=TCGA.UCEC.sampleMap%2FRPPA_RBN&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
options(timeout=200)
utils::download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.UCEC.sampleMap%2FRPPA_RBN.gz", "data/TCGA.UCEC.sampleMap%2FRPPA_RBN.gz")
UCEC_download <- utils::read.table(gzfile("data/TCGA.UCEC.sampleMap%2FRPPA_RBN.gz"), header = T)
rownames(UCEC_download) = UCEC_download[,1]
UCEC_download = t(UCEC_download[,-1])
dim(UCEC_download) # 404 x 131
# Controlled that all Ids end with '01', i.e., no controls. 
all(colnames(UCEC_download) == colnames(BRCA_download)) # Same set of proteins, in same order.
ucec_dat = UCEC_download

# Load OVAC data --------------------------------------------------------------------------

# We use RPPA (replicate-base normalization) data
# Data from https://xenabrowser.net/datapages/?dataset=TCGA.OV.sampleMap%2FRPPA_RBN&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
utils::download.file("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.OV.sampleMap%2FRPPA_RBN.gz", "data/TCGA.OV.sampleMap%2FRPPA_RBN.gz")
OVAC_download <- utils::read.table(gzfile("data/TCGA.OV.sampleMap%2FRPPA_RBN.gz"), header = T)
rownames(OVAC_download) = OVAC_download[,1]
OVAC_download = t(OVAC_download[,-1])
dim(OVAC_download) # 412 x 131
# Controlled that all Ids end with '01', i.e., no controls. 
all(colnames(OVAC_download) == colnames(BRCA_download)) # Same set of proteins, in same order.
ovac_dat = OVAC_download

# Mapping proteins to genes ---------------------------------------------------------------

# Map the names of the antibodies used to identify the proteins to the names of the genes that encode them
# Use file showing which genes the proteins detected by the antibodies are encoded from (can be downloaded from TCGA website)
protein.names = colnames(OVAC_download)
rppa.to.gene = read.table("data/RPPA_to_gene.txt", sep = "\t", stringsAsFactors = F)
mapping.frame = data.frame(protein = protein.names,gene=rppa.to.gene[match(protein.names,rppa.to.gene[,1]),2])
length(unique(mapping.frame$gene)) # 102 unique proteins


# Save data ----------------------------------------

save(ovac_dat, brca_dat, ucec_dat, mapping.frame, file='data/PanCancer_data.RData')





