#' Position Frequency Matrix for CTCF Motif
#' 
#' A Position Frequency Matrix (PFM) for the CTCF Motif
#' obatined from the Jaspar motif database. 
#' 
#' @format A 4x19 matrix giving the nucleotide frequency 
#' at each of 19 positions. 
#' 
#' @source \url{http://jaspar.genereg.net/}
"PFM"

#' BaMM transition probabilities for CTCF Motif
#' 
#' A 2nd order Markov chain representation of the 
#' CTCF motif. The parameters are obtained using the
#' BaMM software and using ENCODE ChIP-seq peaks for the
#' CTCF binding site in MCF-7 cells.  
#' 
#' See also the vignette BaMM Analysis.
#' 
#' @format A list of transition matrices.
#' 
#' @source \url{https://github.com/soedinglab/BaMMmotif} \url{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/}
"BammTransProb"

#' BaMM transition probabilities background
#' 
"BammTransProbBg"

