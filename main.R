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
      if(salidaDeReubicaStrings$parametro1$valor[ejeX-1]==salidaDeReubicaStrings$parametro2$valor[ejeY-1]){ #casoMatch
        
#         # IZQ MAYOR
#         if (matriz[ejeX-1,ejeY]>matriz[ejeX-1,ejeY-1]&matriz[ejeX-1,ejeY]>matriz[ejeX,ejeY-1]){
#           matriz[ejeX,ejeY]<-matriz[ejeX-1,ejeY]-penalizacionGap
#         }
        
        # CRUZADO MAYOR
#         if (matriz[ejeX-1,ejeY-1]>matriz[ejeX-1,ejeY]&matriz[ejeX-1,ejeY-1]>matriz[ejeX,ejeY-1]){
          matriz[ejeX,ejeY]<-matriz[ejeX-1,ejeY-1]+bonoMatch
        # }
        
#         # ARRIBA MAYOR
#         if (matriz[ejeX,ejeY-1]>matriz[ejeX-1,ejeY]&matriz[ejeX,ejeY-1]>matriz[ejeX-1,ejeY-1]){
#           matriz[ejeX,ejeY]<-matriz[ejeX,ejeY-1]-penalizacionGap
#         }

      }#FIN CASO MATCH
      
      else{# casoMissMatch
        
        mayorVecino=max(matriz[ejeX,ejeY-1],matriz[ejeX-1,ejeY-1],matriz[ejeX-1,ejeY])
        if (mayorVecino>penalizacionMiss){
          matriz[ejeX,ejeY]<-mayorVecino-penalizacionGap
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

preprocesoPorLaIzquierda<-function(salidaObtieneSimilitud){
  stringSalida<-list(stringSalida=salidaObtieneSimilitud$stringSalida,puntaje=salidaObtieneSimilitud$puntaje)
  if (salidaObtieneSimilitud$posX==0){
    return (stringSalida)
  }
  else{
    for (posLetra in salidaObtieneSimilitud$posX:0){
      stringSalida[posLetra]<-"_"
    }
    return (stringSalida)
  }
}

obtieneSimilitud<-function(stringPreprocesadoPorLaDerecha,puntaje,posX,posY,matrizMatches){
  puntaje=matrizMatches[posY,posX]
  while(posY>1&posX>1){
    if (matrizMatches[posY,posX]==0){
      return (list(stringSalida=stringPreprocesadoPorLaDerecha,puntaje=puntaje, posX=posX,posY=posY))
    }
    else{
      if (matrizMatches[posX-1,posY]>matrizMatches[posX-1,posY-1]){ #Corresponde a Gap
          stringPreprocesadoPorLaDerecha[posX]<-"_"
          posX<-posX-1
          }
      else{
        posX<-posX-1
        posY<-posY-1
        puntaje=matrizMatches[posY,posX]
      }
      
    }
    
  }#fin while (posY>1&posX>1)
  
  return (list(stringSalida=stringPreprocesadoPorLaDerecha,puntaje=puntaje, posX=posX,posY=posY))
#   if (posX==1| posY==1){
#     return (list(stringSalida=stringPreprocesadoPorLaDerecha,puntaje=puntaje, posX=posX,posY=posY))
#   }
#   
#   else if (matrizMatches[posY-1,posX-1]==0 & matrizMatches[posX-1,posY-1]==0){
#     return (list(stringSalida=stringPreprocesadoPorLaDerecha,puntaje=puntaje, posX=posX,posY=posY))
#   }
#   
#   else{ #posX!=1 && posY!=1
#       if (matrizMatches[posX-1,posY]>matrizMatches[posX-1,posY-1]){ #Corresponde a Gap
#         stringPreprocesadoPorLaDerecha[posY]<-"_"
#         return (obtieneSimilitud(stringPreprocesadoPorLaDerecha,puntaje,posX-1,posY,matrizMatches))
#       }
#     else{
#       return (obtieneSimilitud(stringPreprocesadoPorLaDerecha,puntaje+1,posX-1,posY-1,matrizMatches))
#       
#     }
#  }
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
    print(stringPreprocesado)
    cantidadMax<-cantidadMax-1
  }
  #stringPreprocesado=preprocesoPorLaDerecha(entrada$parametro1$valor,entrada$parametro1$largo,maxPosMatriz[2])
  #stringPreprocesado=obtieneSimilitud(stringPreprocesado,0,5,8,matrizM)
  #stringPreprocesado=preprocesoPorLaIzquierda(obtieneSimilitud(stringPreprocesado,0,9,6,bonoHorizontal,bonoCruzado))
  #print(stringPreprocesado)
  
} else{
  print ("Uno de los parametros no es de tipo String y no se puede ejecutar")
}

#warnings()
#options(error=recover)
