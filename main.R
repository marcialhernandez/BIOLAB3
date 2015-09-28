#install.packages("ggplot2")
#install.packages("psych")
#source ("http://www.bioconductor.org/biocLite.R")
#biocLite()
#biocLite("Biostrings")
#biocLite("multtest")
library (psych)
library(multtest)


inicializaMatriz<-function(largo,ancho){
  return (matrix(0, nrow = largo, ncol = ancho))
}

#Funcion que tiene como entrada 2 parametros de tipo string 
#Y retorna una lista con los valores 'parametro1' y 'parametro2'
#Donde siempre el 'parametro1' tendra largo mayor que 'parametro2'
#tanto parametro1 como parametro2 tienen los atributos 'valor' y 'largo'
#En caso que una entrada no sea String, entonces retorna el bool FALSE
reubicaStrings<-function(string1,string2){
  #Validacion de tipo
  #En caso que si sea
  if (typeof(string1)=="character" & typeof(string2)=="character"){
    
    #Se asocia su largo como atributo
    string1.largo=nchar(string1)
    string2.largo=nchar(string2)
    
    if (string1.largo>=string2.largo){
      #salida.string1=string1
      #salida.string2=string2
      return (list(parametro1=list(valor=unlist(strsplit(string1,"")),largo=string1.largo),parametro2=list(valor=unlist(strsplit(string2,"")),largo=string2.largo)))
    } #fin if (string1.largo>string2.largo)
    
    else{
      return (list(parametro1=list(valor=unlist(strsplit(string2,"")),largo=string2.largo),parametro2=list(valor=unlist(strsplit(string1,"")),largo=string1.largo)))
    } #fin else
    
  } #fin if (typeof(string1)=="character" & typeof(string2)=="character") 
  else{
    return(FALSE)
  }
}

generaMatrizMatch<-function(salidaDeReubicaStrings,bonoMatch,penalizacionMiss,penalizacionGap){
  #Genero matriz
  matriz=inicializaMatriz(salidaDeReubicaStrings$parametro2$largo+1,salidaDeReubicaStrings$parametro1$largo+1)
  colnames(matriz)<-c("-",salidaDeReubicaStrings$parametro1$valor)
  rownames(matriz)<-c("-",salidaDeReubicaStrings$parametro2$valor)
  
  for (ejeX in 1:salidaDeReubicaStrings$parametro2$largo+1){
    
    for (ejeY in 1:salidaDeReubicaStrings$parametro1$largo+1){
      
      #Valida correcto cruce de informacion para la comparacion
      #print (paste("Pos (",toString(ejeX),",",toString(ejeY),")",": Cruce ->",salidaDeReubicaStrings$parametro1$valor[ejeX-1],"-",salidaDeReubicaStrings$parametro2$valor[ejeY-1]))
      
      #CASO MATCH
      if(salidaDeReubicaStrings$parametro1$valor[ejeY-1]==salidaDeReubicaStrings$parametro2$valor[ejeX-1]){ #casoMatch
        #print (paste("Pos (",toString(ejeX),",",toString(ejeY),")",": Cruce ->",salidaDeReubicaStrings$parametro1$valor[ejeY-1],"-",salidaDeReubicaStrings$parametro2$valor[ejeX-1]))
        
        
#         # IZQ MAYOR
#         if (matriz[ejeX-1,ejeY]>matriz[ejeX-1,ejeY-1]&matriz[ejeX-1,ejeY]>matriz[ejeX,ejeY-1]){
#           matriz[ejeX,ejeY]<-matriz[ejeX-1,ejeY]-penalizacionGap
#         }
        
        # CRUZADO MAYOR
#         if (matriz[ejeX-1,ejeY-1]>matriz[ejeX-1,ejeY]&matriz[ejeX-1,ejeY-1]>matriz[ejeX,ejeY-1]){
          #matriz[ejeX,ejeY]<-matriz[ejeX-1,ejeY-1]+bonoMatch
        matriz[ejeX,ejeY]<-max(matriz[ejeX-1,ejeY-1],matriz[ejeX-1,ejeY],matriz[ejeX,ejeY-1])+bonoMatch
        # }
        
#         # ARRIBA MAYOR
#         if (matriz[ejeX,ejeY-1]>matriz[ejeX-1,ejeY]&matriz[ejeX,ejeY-1]>matriz[ejeX-1,ejeY-1]){
#           matriz[ejeX,ejeY]<-matriz[ejeX,ejeY-1]-penalizacionGap
#         }

      }#FIN CASO MATCH
      
      else{# casoMissMatch
        
        mayorVecino=max(matriz[ejeX,ejeY-1]-penalizacionGap,matriz[ejeX-1,ejeY-1],matriz[ejeX-1,ejeY]-penalizacionGap)
        if (mayorVecino>penalizacionMiss){
          matriz[ejeX,ejeY]<-mayorVecino-penalizacionMiss
        }
        
        else{
          matriz[ejeX,ejeY]<-0
        }
      }
      
      
     
    } #Fin for (posLetra1 in 1:salidaDeReubicaStrings$parametro1$largo)
  }#Fin for (posLetra2 in 1:salidaDeReubicaStrings$parametro2$largo)
  
  return (matriz)
}

validaStrings<-function(salidaDeReubicaStrings){
  if (typeof(salidaDeReubicaStrings)=="logical"){
    return (FALSE)
  }
  else{
    return (TRUE)
  }
}

preprocesoPorLaDerecha<-function(stringEntrada,largo,posXActual){
  stringSalida<-stringEntrada
  if (posXActual>=largo){
    return (stringSalida)
  }
  
  else{
    #print (posXActual+1:largo)
    while(posXActual<largo){
      posXActual<-posXActual+1
      stringSalida[posXActual]<-"_"
    }
    return (stringSalida)
  }
}

preprocesoPorLaIzquierda<-function(stringEntrada,posXActual){
  if (posXActual==1){
    return (stringEntrada)
  }
  else{
    stringEntrada[posXActual]<-"_"
    return (preprocesoPorLaIzquierda(stringEntrada,posXActual-1))
  }
}

obtieneSimilitud<-function(stringPreprocesadoPorLaDerecha,posX,posY,matrizMatches){
  if (matrizMatches[posX,posY]==0){
    return (list(stringSalida=stringPreprocesadoPorLaDerecha, posX=posX))
  }
  else{
    #CasoMatch
    if (matrizMatches[posX-1,posY-1]>=matrizMatches[posX-1,posY]){
      return (obtieneSimilitud(stringPreprocesadoPorLaDerecha,posX-1,posY-1,matrizMatches))
    }
    
    #Caso Gap
    else{
      stringPreprocesadoPorLaDerecha[posX]<-"_"
      return (obtieneSimilitud(stringPreprocesadoPorLaDerecha,posX-1,posY,matrizMatches))
    }
  }
}

#Parametros de entrada
#Se asigna string entrada 1
stringEntrada1<-"ACACACTA" #"ATTAAG"
#Se asigna string entrada 2
stringEntrada2<-"AGCACACA" #"GATTACAAG"
#Bonos y Penalizaciones
bonoMatch<-2
penalizacionMiss<-1
penalizacionGap<-1


#----------------------

entrada<-reubicaStrings(stringEntrada1,stringEntrada2);
if (validaStrings(entrada)==TRUE){
  matrizM=generaMatrizMatch(entrada,bonoMatch,penalizacionMiss,penalizacionGap)
  maxPosMatriz<-which(matrizM == max(matrizM), arr.ind = TRUE)
  cantidadMax<-length(maxPosMatriz[,1]) #1 es y, 2 es x
  print(matrizM)
  
  #Para cada valor Maximo encontrado
  while(cantidadMax>0){
    
    #Validacion - imprime cada columna
    #print (maxPosMatriz[cantidadMax,])
    #Rellena con gaps por la derecha cuando lo requiera
    stringPreprocesado=preprocesoPorLaDerecha(entrada$parametro1$valor,entrada$parametro1$largo,maxPosMatriz[cantidadMax,][2])
    
    #Termina el proceso y se imprime
    print("--------------------------------------------------------------------------------------")
    print(entrada$parametro1$valor)
    print(entrada$parametro2$valor)
    stringPreprocesado<-obtieneSimilitud(stringPreprocesado,maxPosMatriz[cantidadMax,][1],maxPosMatriz[cantidadMax,][2],matrizM)
    stringPreprocesado<-preprocesoPorLaIzquierda(stringPreprocesado$stringSalida,stringPreprocesado$posX)
    print(stringPreprocesado)
    cantidadMax<-cantidadMax-1
  }
  
} else{
  print ("Uno de los parametros no es de tipo String y no se puede ejecutar")
}

#warnings()
#options(error=recover)
