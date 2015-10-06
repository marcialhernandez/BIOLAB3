#install.packages("ggplot2")
#install.packages("psych")
#source ("http://www.bioconductor.org/biocLite.R")
#biocLite()
#biocLite("ALL")
#biocLite("Biostrings")
#biocLite("multtest")
library(multtest)
library(Biostrings)

claseDNA<-function(entrada){
  dnaSeq = readDNAStringSet(entrada)
  return (list(dnaSeq=dnaSeq,
               largo=length(dnaSeq), #Establece el largo de la secuencia
               complemento=reverseComplement(dnaSeq), #Genera una cadena complementaria reversa
               frecuenciaAlfabetica=alphabetFrequency(dnaSeq), #Indica la cardinalidad de cada caracter
               dinucleotideFrequency(dnaSeq)))
}

rattus = claseDNA("FASTA/rattusnorvegicuschromosomeY.fasta")
print (rattus)