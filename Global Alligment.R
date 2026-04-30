
#This script is made to analyze 2 samples. 
#Forward and Reverse by each sample, in total 4 files.

####- 1-- Packages -------------------------------####
#-##################################################-#

pacman::p_load(
  sangerseqR, Biostrings, dplyr)

#-##################################################-#
####- 2-- Set up your working directory ----------####
#-##################################################-#

setwd("COPY_YOUR_PATH_HERE")

#-##################################################-#
####- 3-- Reading AB1 files-----------------------####
#-##################################################-#

# Reading AB1 files (TOP band F/R)
s_f <- read.abif("1-5F_ID3.ab1")
s_r <- read.abif("1-178R_ID4.ab1")
sseq_f <- sangerseq(s_f)
sseq_r <- sangerseq(s_r)


# Reading AB1 files (BOT band F/R)
s_f <- read.abif("2-5F_ID5.ab1")
s_r <- read.abif("2-178R_ID6.ab1")
sseq_f <- sangerseq(s_f)
sseq_r <- sangerseq(s_r)



#-##################################################-#
####- 4-- F - Selecting Trimming positions -------####
#-##################################################-#

## Manual visualization in SnapGene and star= and end= possition annotation 
# Trimming AB1 files (TOP band F)

start_f <- 18
end_f   <- 139


# Reading AB1 files (BOT band F)

#start_f <- 19
#end_f   <- 142


#-##################################################-#
####- 5-- F - Base calling and trimming ----------####
#-##################################################-#

# Base-calling depends on ratio (30%)
seq_f_called <- makeBaseCalls(sseq_f, ratio = 0.33)
raw_f <- as.character(primarySeq(seq_f_called))
trimmed_f <- substring(raw_f, first = start_f, last = end_f)

PolyPeakParser()
#-##################################################-#
####- 6-- R - Selecting Trimming positions -------####
#-##################################################-#


## Manual visualization in SnapGene and star= and end= possition annotation 
#Keep in mind that selected trimming possitions in SnapGene are the same as original for reverse product

# Trimming AB1 files (TOP band R)
start_r <- 28
end_r   <- 150



# Trimming AB1 files (BOT band R)
#start_r <- 30
#end_r   <- 145

#-##################################################-#
####- 7-- R - Base calling and trimming ----------####
#-##################################################-#

# Base-calling depends on ratio (30%)
seq_r_called <- makeBaseCalls(sseq_r, ratio = 0.33)
raw_r <- as.character(primarySeq(seq_r_called))

trimmed_r <- substring(raw_r, first = start_r, last = end_r)


#-##################################################-#
####- 8-- R - Reverse Complement ----------------####
#-##################################################-#


# Reverse complement de la R both Top and Bottom reverse trimmed
trimmed_r_rc <- as.character(reverseComplement(DNAString(trimmed_r)))


#-##################################################-#
####- 9-- F/R - Alignment -----------------------####
#-##################################################-#

# Alinear 
aln <- pairwiseAlignment(DNAString(trimmed_f), DNAString(trimmed_r_rc), type="global")


# How was the alignment executed
writePairwiseAlignments(aln)

#-##################################################-#
####- 10-- F/R - Consensus -----------------------####
#-##################################################-#

#  Preliminary Aligment result
aligned(aln)

# Aligment in a matrix 
consensus <- consensusMatrix(aln) # inspeccionar #
consensus

# Aligment consensus (stronger) based on the consensusMatrix
x <- consensusString(aln)
x
# Real bases numer of consensus sequence
nchar(x)

# Aligment Score 
score(aln)



#-##################################################-#
####- 11-- F/R - Saving as Fasta file-------------####
#-##################################################-#

# Converting consensus sequence into a DNAStringSet object

consensus_fasta <- DNAStringSet(x)

names(consensus_fasta) <- "TopBand"  # Header FASTA name 

# Saving fasta file
writeXStringSet(consensus_fasta, filepath = "TopBand.fasta")
















