#Installing packages
install.packages("tidyverse")
install.packages("ape")
ninstall.packages("pegas")
install.packages("phangorn")
install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("ggtree")     

#Resetting and calling packages
rm(list = ls())
library(tidyverse)
library(ape)
library(pegas)
library(phangorn)
library(Biostrings)
library(ggtree)

#Get Working directory
getwd()

#Set Working directory
setwd("C:/Users/adero/OneDrive/Documents/Exercise 3B")

#Location = NIGERIA
#look at the input fasta file
file.show("nigeria.fasta")

#read in the file
flyDNA_nigeria<- read.dna("nigeria.fasta" , format = "fasta")

#Get summary details on sequence
myflyalign_nigeria <- as.alignment(flyDNA_nigeria)
myflyalign_nigeria$seq

#Simplified view of data
alview(flyDNA_nigeria)


#Read the FASTA file to calculate sequence legth
fasta_data <- readDNAStringSet("nigeria.fasta")

#length of the first sequence
first_sequence <- fasta_data[[1]]

# Calculate the length of the first sequence
sequence_length <- nchar(as.character(first_sequence))
sequence_length


#Calculate segregating sites
nuc.div(flyDNA_nigeria)
S <- length(seg.sites(flyDNA_nigeria))
L <- 1488
s <- S/L
s

#Calculate segregating sites
tajima.test(flyDNA_nigeria)

#Count number of segregating sites
S



#Location = BURKINA FASO
#look at the input fasta file
file.show("Burkina_Faso.fasta")

#read in the file
flyDNA_B_Faso<- read.dna("Burkina_Faso.fasta" , format = "fasta")

#Get summary details on sequence
myflyalign_B_Faso <- as.alignment(flyDNA_B_Faso)
myflyalign_B_Faso$seq

#Simplified view of data
alview(flyDNA_B_Faso)

#Calculate segregating sites
nuc.div(flyDNA_B_Faso)
S <- length(seg.sites(flyDNA_B_Faso))
L <- 1488
s <- S/L
s

#Calculate segregating sites
tajima.test(flyDNA_B_Faso)

#Count number of segregating sites
S


#Location = Cote d'Ivoire
#look at the input fasta file
file.show("Cote_d_ivoire.fasta")

#read in the file
flyDNA_cdi<- read.dna("Cote_d_ivoire.fasta" , format = "fasta")

#Get summary details on sequence
myflyalign_cdi <- as.alignment(flyDNA_cdi)
myflyalign_cdi$seq

#Simplified view of data
alview(flyDNA_cdi)

#Calculate segregating sites
nuc.div(flyDNA_cdi)
S <- length(seg.sites(flyDNA_cdi))
L <- 1488
s <- S/L
s

#Calculate segregating sites
tajima.test(flyDNA_cdi)

#Count number of segregating sites
S



#Location = Nigeria, Burkina Faso, Cote d'Ivoire
#look at the input fasta file
file.show("Combined_locations.fasta")

#read in the file
flyDNA_combined<- read.dna("Combined_locations.fasta" , format = "fasta")

#Get summary details on sequence
myflyalign_combined <- as.alignment(flyDNA_combined)
myflyalign_combined$seq

#Simplified view of data
alview(flyDNA_combined)

#Calculate segregating sites
nuc.div(flyDNA_combined)
S <- length(seg.sites(flyDNA_combined))
L <- 1488
s <- S/L
s

#Calculate segregating sites
tajima.test(flyDNA_combined)

#Count number of segregating sites
S



# Generate a phylogenetic tree with an outgroup
fasta_file <- "Combined_with_outgroup.fasta"
aligned_sequences <- read.dna(fasta_file, format = "fasta")

# Calculate the distance matrix
distance_matrix <- dist.dna(aligned_sequences)

# Construct the phylogenetic tree using the Neighbor Joining (NJ) method
nj_tree <- nj(distance_matrix)

# Increase the branch length
nj_tree$edge.length <- nj_tree$edge.length * 1

# Plot the phylogenetic tree without the title
# The only thing wrong with this tree is that the text is too big. 
# We just need to make one small change.
# plot(nj_tree, cex = 0.6) # Adjust text size

plot(nj_tree, cex = 0.2) # Made the text smaller. 

## For the title, we have to get some distance parameters.
# This is the x-axis range of the plot:
x_range <- range(nj_tree$xx)

# Calculate the midpoint of the x-axis range.
midpoint <- mean(x_range)

# Now we can add the title at the bottom of the plot.
# This works on my machine, but you may need to adjust the value of the "line" argument to greater than 0 on yours.
# It determines how far the text is from the bottom of the screen.
mtext(expression(paste("Phylogenetic tree of ", italic("Bactrocera dorsalis")," from the combined locations")),
      side = 1, line = 0, cex = 0.5, at = midpoint)


# Save the plot as a PNG file
png("phylogenetic_tree.png", width = 800, height = 1000) # Adjust the width and heigh parameters so that it isn't too compacted.
plot(nj_tree, cex = 0.6) 
dev.off() 

##### End Simple Tree #####

##### Begin ggtree  ######
rm(list = ls())
library(tidyverse)
library(ape)
#library(pegas)
#library(phangorn)
library(Biostrings)
library(ggtree)
library(ggplot2)

fly_seqs <- "Combined_with_outgroup.fasta"
aligned_sequences <- read.dna(fly_seqs, format = "fasta")

# Calculate the distance matrix
distance_matrix <- dist.dna(aligned_sequences)

# Construct the phylogenetic tree using the Neighbor Joining (NJ) method
nj_tree <- nj(distance_matrix)

# You can plot a basic tree with ggtree this way.
q <- ggtree(nj_tree) +
  geom_tiplab(size=1.2, align=FALSE, linesize=.5)  
print(q)

# Add a title:
q <- ggtree(nj_tree) +
  geom_tiplab(size=1.2, align=FALSE, linesize=.5)  +
  labs(caption = expression(paste("Phylogenetic tree of ", italic("Bactrocera dorsalis"), "from the combined locations"))) +
  theme(plot.caption = element_text(hjust = 0.5, size = rel(.75)))  # Adjust position and size
print(q)


# You can also change the layout of the plot in many ways.
q <- ggtree(nj_tree, layout = "roundrect") +
  geom_tiplab(size=1.2, align=FALSE, linesize=.5)  +
  labs(caption = "Put the plot title here") +
  theme(plot.caption = element_text(hjust = 0.5, size = rel(.75)))  # Adjust position and size
print(q)

# The possible layouts  are 'rectangular', 'dendrogram', 'slanted', 'ellipse', 'roundrect',
# 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle', 'daylight' or 'ape'
# Try changing the layout value of the plot below to some of these options.
q <- ggtree(nj_tree, layout = "roundrect", branch.length = 'none') +
  geom_tiplab(size=1.2, align=FALSE, linesize=.5)  +
  labs(caption = "Put the plot title here") +
  theme(plot.caption = element_text(hjust = 0.5, size = rel(.75)))  
print(q)


# Example: fan shaped cladogram
q <- ggtree(nj_tree, layout = "fan", branch.length = 'none') +
  geom_tiplab(size=1.5, align=FALSE, linesize=.5)  +
  labs(caption =  expression(paste("Phylogenetic tree of ", italic("Bactrocera dorsalis"), " from the combined locations"))) +
  theme(plot.caption = element_text(hjust = 0.5, size = rel(.75)))  
print(q)


## You can also add colors, and do a lot of other things with ggtree!

