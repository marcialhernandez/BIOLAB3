library(Biostrings)

waterSmithSencillo<-function(s1,s2,puntajeMatch,descuentoMiss){
  mat <- nucleotideSubstitutionMatrix(match = puntajeMatch, mismatch = descuentoMiss, baseOnly = TRUE)
  return (pairwiseAlignment(s1, s2, type = "local", substitutionMatrix = mat))
}

waterSmithDetallado<-function(s1,s2,puntajeMatch,descuentoMiss){
  mat <- nucleotideSubstitutionMatrix(match = puntajeMatch, mismatch = descuentoMiss, baseOnly = TRUE)
  resultado<-pairwiseAlignment(s1, s2, type = "local", substitutionMatrix = mat)
  return (list(sumario=summary(resultado),alineado=aligned(resultado)))
}

print ("WaterSmith from BioConductor")
print (waterSmithSencillo(DNAString("GAATTCCTACTACGGAATTCCCCTCCCATAATTCCTACTACGA"),DNAString("GAATTCCTACTACGAAGAATTCCTACTACGAAACTACGAAAATTCCTACTACGA"),1,-1))