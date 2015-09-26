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

generaMatrizMatch2<-function(salidaDeReubicaStrings,bonoMatch,penalizacionMiss,penalizacionGap){
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

generaMatrizMatch<-function(salidaDeReubicaStrings,bonoMatch,penalizacionMissMatch,penalizacionGap){
  #Genero matriz
  matriz=inicializaMatriz(salidaDeReubicaStrings$parametro2$largo,salidaDeReubicaStrings$parametro1$largo)
  #print (matriz)
  ejeX<-0
  ejeY<-0
  
  for (posLetra2 in 1:salidaDeReubicaStrings$parametro2$largo){
    ejeX<-ejeX+1
    ejeY<-0
    for (posLetra1 in 1:salidaDeReubicaStrings$parametro1$largo){
      ejeY<-ejeY+1
      #Valida correcto cruce de informacion para la comparacion
      #print (paste("Pos (",toString(ejeX),",",toString(ejeY),")",": Cruce ->",salidaDeReubicaStrings$parametro1$valor[posLetra1],"-",salidaDeReubicaStrings$parametro2$valor[posLetra2]))
      
      #-Formacion Matriz----------------------------------------------------------------# 
      if (ejeX==1){ 
        
        if (ejeY==1){ #ejeX==1 && ejeY==1
          
          if(salidaDeReubicaStrings$parametro1$valor[posLetra1]==salidaDeReubicaStrings$parametro2$valor[posLetra2]){ #casoMatch
            print
            matriz[ejeX,ejeY]<-matriz[ejeX,ejeY]+bonoMatch
          }
          
          else{# casoMissMatch
            matriz[ejeX,ejeY]<-0
          }
          
        }#Fin caso #ejeX==1 && ejeY==1
        
        else{ #ejeX==1 && ejeY!=1
          
          if(salidaDeReubicaStrings$parametro1$valor[posLetra1]==salidaDeReubicaStrings$parametro2$valor[posLetra2]){ #casoMatch
            matriz[ejeX,ejeY]<-matriz[ejeX,ejeY-1]+bonoMatch
          }
          
          else{# casoMissMatch
            if (matriz[ejeX,ejeY-1]-penalizacionMissMatch>0){
              matriz[ejeX,ejeY]<-matriz[ejeX,ejeY-1]-penalizacionMissMatch
            }
            else{
              matriz[ejeX,ejeY]<-0
            }
          }
          
        }# Fin caso #ejeX==1 && ejeY!=1
        
      }
      
      else if(ejeY==1){ #ejeX!=1 && ejeY==1
        
        if(salidaDeReubicaStrings$parametro1$valor[posLetra1]==salidaDeReubicaStrings$parametro2$valor[posLetra2]){ #casoMatch
          matriz[ejeX,ejeY]<-matriz[ejeX-1,ejeY]+bonoMatch
        }
        
        else{# casoMissMatch
          if (matriz[ejeX-1,ejeY]-penalizacionMissMatch>0){
            matriz[ejeX,ejeY]<-matriz[ejeX-1,ejeY]-penalizacionMissMatch
          }
          else{
            matriz[ejeX,ejeY]<-0
          }
        }
        
      }#Fin #ejeX!=1 && ejeY==1
      
      else{ #ejeX!=1 && ejeY!=1
        
        if(matriz[ejeX-1,ejeY-1]>matriz[ejeX-1,ejeY] & matriz[ejeX-1,ejeY-1]>matriz[ejeX,ejeY-1]){ #Caso Cruzado
          
          if(salidaDeReubicaStrings$parametro1$valor[posLetra1]==salidaDeReubicaStrings$parametro2$valor[posLetra2]){ #casoMatch
            matriz[ejeX,ejeY]<-matriz[ejeX-1,ejeY-1]+bonoMatch
          }
          
          else{# casoMissMatch
            if (matriz[ejeX-1,ejeY-1]>penalizacionMissMatch){
              matriz[ejeX,ejeY]<-matriz[ejeX-1,ejeY-1]-penalizacionMissMatch
            }
            else{
              matriz[ejeX,ejeY]<-0
            }
          }
          
        } #Fin Caso Cruzado
        
        else if (matriz[ejeX-1,ejeY]>matriz[ejeX-1,ejeY-1] & matriz[ejeX-1,ejeY]>matriz[ejeX,ejeY-1]){ #Caso Izq (Gap)
          
          if(salidaDeReubicaStrings$parametro1$valor[posLetra1]==salidaDeReubicaStrings$parametro2$valor[posLetra2]){ #casoMatch
            matriz[ejeX,ejeY]<-matriz[ejeX-1,ejeY]-penalizacionGap
          }
          
          else{# casoMissMatch
            if (matriz[ejeX-1,ejeY]>penalizacionMissMatch){
              matriz[ejeX,ejeY]<-matriz[ejeX-1,ejeY]-penalizacionMissMatch
            }
            else{
              matriz[ejeX,ejeY]<-0
            }
          }
          
        } #Fin Caso Izq (Gap)
        
        else{ #Caso Arriba (Gap)
          
          if(salidaDeReubicaStrings$parametro1$valor[posLetra1]==salidaDeReubicaStrings$parametro2$valor[posLetra2]){ #casoMatch
            matriz[ejeX,ejeY]<-matriz[ejeX,ejeY-1]-penalizacionGap
          }
          
          else{# casoMissMatch
            if (matriz[ejeX,ejeY-1]>penalizacionMissMatch){
              matriz[ejeX,ejeY]<-matriz[ejeX,ejeY-1]-penalizacionMissMatch
            }
            else{
              matriz[ejeX,ejeY]<-0
            }
          }
          
        }# Fin Caso Arriba (Gap)
        
      }#Fin ejeX!=1 && ejeY!=1
      
      #---------------------------------------------------------------------------------#
      
    } #Fin for (posLetra1 in 1:salidaDeReubicaStrings$parametro1$largo)
  }#Fin for (posLetra2 in 1:salidaDeReubicaStrings$parametro2$largo)
  
  return (matriz)
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
  matrizM=generaMatrizMatch2(entrada,bonoMatch,penalizacionMiss,penalizacionGap)
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
