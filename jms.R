#Carga BioStrings...
library(Biostrings);

#Parametros de Entrada: Secuencias, bonus y castigo.
Seq1 <- "GAATTCCTACTACGAAGAATTCCTACTACGAAACTACGAAAATTCCTACTACGA";
Seq2 <- "GAATTCCTACTACGGAATTCCCCTCCCATAATTCCTACTACGA";
#Seq1 <- "PELICAN";
#Seq2 <- "COELACANTH";
valorMatch <- 1;
valorMismatch <- -1;
valorGap <- 0;

#Parametro de control
escribirMatriz <- 1;

#Algoritmo.

#Paso 1: Funciones de apoyo, busca el mayor de los tres anteriores.
maximo <- function(i,j) {
  n <- valores[i, j];
  if (n < valores[i+1, j]) {
    n <- valores[i+1, j] + valorGap;
  }
  if (n < valores[i, j+1]) {
    n <- valores[i, j+1] + valorGap;
  }
  return(n); 
}

#Paso 2: Crear la matriz inicial y sus datos
columnas <- nchar(Seq1);
filas <- nchar(Seq2);
valores <- matrix(nrow = (filas+1), ncol = (columnas+1));
               
#Valores de control para maximo local.
maxval <- 0;
maxi <- 0;
maxj <- 0;

#Paso 3: Inicializar primera fila y primera columna
for (i in 1:(filas+1)) {
  valores[i,1] <- 0;
}
for (j in 1:(columnas+1)) {
  valores[1,j] <- 0;
}

#Paso 4: Calcular matchs / mismatchs
# y conocer el punto de maximo encuentro
for (i in 1:filas) {
  for (j in 1:columnas) {
    n <- maximo(i,j);
    if (substr(Seq2, i, i) == substr(Seq1, j, j)) {
      n <- n + valorMatch;
    } else {
      n <- n + valorMismatch;
      if (n < 0) {
        n <- 0;
      }
    }
    valores[i+1, j+1] <- n;
    
    #verifica si alcanzo maximo local.
    if (n > maxval) {
      maxval <- n;
      maxi <- i;
      maxj <- j;
    }
  }
}

#Paso 5: Escribe datos de matriz.
if (escribirMatriz == 1) {
  print("Matriz de datos:");
  colnames(valores) <- strsplit(paste("-", Seq1, sep=""), "")[[1]];
  rownames(valores) <- strsplit(paste("-", Seq2, sep=""), "")[[1]];
  print(valores);
}

#Paso6: Obtener las salidas
sseq1 <- ""; #Asociada a Seq1, va por columnas.
sseq2 <- ""; #Asociada a Seq2, va por filas.

direccion <- "";
j <- maxj;
i <- maxi;
n <- maxval;
lastdir <- "*";
lastone <- 2;
while(i > 1 & j > 1 & (n + lastone) > 0) {
  if (n > 0) {
    n <- maximo(i,j);
  }
  if (n == valores[i,j]) {
    sseq1 = paste(substr(Seq1, j, j), sseq1, sep="");
    sseq2 = paste(substr(Seq2, i, i), sseq2, sep="");
    i <- i - 1;
    j <- j - 1;
    lastdir <- "D";
  } else if (n == valores[i, j+1]) {
    sseq1 = paste("-", sseq1, sep="");
    sseq2 = paste(substr(Seq2, i, i), sseq2, sep="");
    lastdir <- "A";
    i <- i - 1;
  } else {
    sseq1 = paste(substr(Seq1, j, j), sseq1, sep="");
    sseq2 = paste("-", sseq2, sep="");
    lastdir <- "I";
    j <- j - 1;
  }
  direccion <- paste(direccion, lastdir);
  if (n == 0) {
    lastone <- lastone - 1;
  }
}
if (n > 0) {
  sseq1 = paste(substr(Seq1, j, j), sseq1, sep="");
  sseq2 = paste(substr(Seq2, i, i), sseq2, sep="");
}

#Paso 7: Imprimir resultados.
print("Algorithm from Jose Santibanez")
print(paste("MaxScore:",maxval))
print(paste("Secuencia 1:", sseq1))
print(paste("Secuencia 2:", sseq2))