#Maria Teresa Escobar
#30/09/2015
#Alineamientos de Secuencias


#Secuencias ingresadas, donde la secuencia 2 es la mas larga

string1<-"GAATTCCTACTACGGAATTCCCCTCCCATAATTCCTACTACGA"
string2<-"GAATTCCTACTACGAAGAATTCCTACTACGAAACTACGAAAATTCCTACTACGA"
penalty<-1 #Penalizacion para la matriz.

#Convierte los strings en vectores
seq1<-unlist(strsplit(string1,""))
seq2<-unlist(strsplit(string2,""))

#strings para nombrar filas y columnas en la matriz
fila<-"-GAATTCCTACTACGGAATTCCCCTCCCATAATTCCTACTACGA"
columna<-"-GAATTCCTACTACGAAGAATTCCTACTACGAAACTACGAAAATTCCTACTACGA"
fil<-unlist(strsplit(fila,""))
col<-unlist(strsplit(columna,""))

#Estos vectores son para la salida de la comparacion de secuencia
A<-c()
B<-c()

#Creacion de la matriz de panalizacion
H=matrix(data = NA, nrow = length(seq1), ncol = length(seq2));

H<-matrix(c(0),nrow = length(seq1)+1, ncol = length(seq2)+1)
n=(length(seq1)+1)
m=(length(seq2)+1)

#Se ha optado para este ejemplo penalizacion "2"

for(i in 2:n)
  for(j in 2:m)
  {
    if(seq1[i-1]==seq2[j-1])
      H[i,j] <- (penalty + H[i-1,j-1]) #penalty: Opcion de Cambiar la penalizacion. 
    else 
    { 
      aux<-max(H[i-1,j-1],H[i,j-1],H[i-1,j])
      aux<-aux-1
      if(aux < 0){
        aux<-0
      }
      H[i,j] <- (aux)
    }
  }

#Algortimo para la secuencias de salida con gaps

i<-n-1
j<-m-1
k<-1

while(k < m)  
{  
  if(seq1[i] == seq2[j]) #Los nucleotidos son iguales
  {
    A[k]<-c(seq1[i])
    B[k]<-c(seq2[j])
    i<-i-1
    j<-j-1
    k<-k+1
  }else   #Los nucleotidos son distintos y en base a la matriz sigue hacia "el lado" 
  {
    if(H[i+1,j]>=H[i,j+1])
    {
      if(H[i+1,j]==0) #Si el alineamiento terminas antes de H[1,1] quiebra el ciclo.
      {
        A[k]<-c(seq1[i])
        B[k]<-c("-")
        break
      }
      A[k]<-c('-')
      B[k]<-c(seq2[j])
      j<-j-1
      k<-k+1
    }else  #Los nucleotidos son distintos y en base a la matriz sigue hacia "arriba"
    {
      A[k]<-c(seq1[i])
      B[k]<-c('-')
      i<-i-1
      k<-k+1
      m<-m+1
    }
  }
}

#Asignar nombre a las fila de la Matriz
rownames(H)<- c(fil)
#Asignar nombre a las columnas de la Matriz
colnames(H)<- c(col)
#Muestra la Matriz
(H)

#Creacion de la salida de la comparacion con forma de matriz
k=length(A)
S=matrix(data = NA, nrow = 2, ncol = length(A));
j<-1

for(i in k:1)
{
  S[1,j]=A[i]
  S[2,j]=B[i]
  j<-j+1
}

colnames(S)<-c(1:length(A))
rownames(S)<- c("A","B")
score<-H[length(fil),length(col)]

#Salida final de Matriz de penalizacion, alineamientos y score
(H)
(S)

print("Algorithm from Maria Escobar")
print(paste0("El score de esta alineacion (con penalizacion 1) fue:  ",score))
print(A)
print(B)

