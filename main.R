# Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)
# Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use. 
# NonCommercial — You may not use the material for commercial purposes.
# NoDerivatives — If you remix, transform, or build upon the material, you may not distribute the modified material. 
# No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits. 
# https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode#languages
# Copyright (C)
# author: Marcial Hernandez Sanchez
# date: 30/9/2015
# University of Santiago, Chile (Usach)

#Historial de cambios: https://github.com/marcialhernandez/BIOLAB3/commits/master

#install.packages("ggplot2")
#install.packages("psych")
#source ("http://www.bioconductor.org/biocLite.R")
#biocLite()
#biocLite("Biostrings")
#biocLite("multtest")
library (psych)
library(multtest)

#Funciones###################################################################

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
      
      if(salidaDeReubicaStrings$parametro1$valor[ejeY-1]==salidaDeReubicaStrings$parametro2$valor[ejeX-1]){ #casoMatch
        #print (paste("Pos (",toString(ejeX),",",toString(ejeY),")",": Cruce ->",salidaDeReubicaStrings$parametro1$valor[ejeY-1],"-",salidaDeReubicaStrings$parametro2$valor[ejeX-1]))
        
        matriz[ejeX,ejeY]<-max(matriz[ejeX-1,ejeY-1],matriz[ejeX-1,ejeY],matriz[ejeX,ejeY-1])+bonoMatch

      }#FIN CASO MATCH
      
      else{# casoMissMatch
        
        mayorVecino=max(matriz[ejeX,ejeY-1],matriz[ejeX-1,ejeY-1],matriz[ejeX-1,ejeY])
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

obtieneSimilitud<-function(entrada1,entrada2,salida1,salida2,posX,posY,matrizMatches){
  if (matrizMatches[posX,posY]==0){
  #if (posX==1 | posY==1){
    return (list(salida1=salida1, salida2=salida2))
  }
  else{
    #Caso Cruzado Mayor #MATCH
    maximoVecino=max(matrizMatches[posX-1,posY-1],matrizMatches[posX-1,posY],matrizMatches[posX,posY-1])
    if (matrizMatches[posX-1,posY-1]==maximoVecino){
      salida1<-append(entrada1[posY-1],salida1)
      salida2<-append(entrada2[posX-1],salida2)
      return (obtieneSimilitud(entrada1,entrada2,salida1,salida2,posX-1,posY-1,matrizMatches))
    }
    #Caso Izquierda Mayor #GAP
    else if (matrizMatches[posX-1,posY]==maximoVecino){
      salida1<-append("_",salida1)
      return (obtieneSimilitud(entrada1,entrada2,salida1,salida2,posX-1,posY,matrizMatches))
  }
    #Caso Arriba Mayor #indel
    else{
      salida2<-append("_",salida2)
      return (obtieneSimilitud(entrada1,entrada2,salida1,salida2,posX,posY-1,matrizMatches))
    }
  }
}

smithWaterman<-function(string1,string2, pMiss,pMatch){
  #Parametros de entrada
  #Se asigna string entrada 1
  stringEntrada1<-string1 #"ATTAAG"
  #Se asigna string entrada 2
  stringEntrada2<-string2 #"GATTACAAG"
  #Bonos y Penalizaciones
  bonoMatch<-pMatch
  penalizacionMiss<-pMiss
  penalizacionGap<-1
  entrada<-reubicaStrings(stringEntrada1,stringEntrada2);
  
  if (validaStrings(entrada)==TRUE){
    matrizM=generaMatrizMatch(entrada,bonoMatch,penalizacionMiss,penalizacionGap)
    maxPosMatriz<-which(matrizM == max(matrizM), arr.ind = TRUE)
    cantidadMax<-length(maxPosMatriz[,1]) #1 es y, 2 es x
    print("--------------------------------------------------------------------------------------")
    print("Matriz Miss-Match")
    print("--------------------------------------------------------------------------------------")
    print(matrizM)
    print("--------------------------------------------------------------------------------------")
    print(paste("Puntaje Max:",matrizM[maxPosMatriz][1]))
    #Para cada valor Maximo encontrado
    while(cantidadMax>0){
      
      print("--------------------------------------------------------------------------------------")
      stringPreprocesado<-obtieneSimilitud(entrada$parametro1$valor,entrada$parametro2$valor,c(),c(),maxPosMatriz[cantidadMax,][1],maxPosMatriz[cantidadMax,][2],matrizM)
      print (stringPreprocesado)
      cantidadMax<-cantidadMax-1
    }
    
  } else{
    print ("Uno de los parametros no es de tipo String y no se puede ejecutar")
  }
  #warnings()
  #options(error=recover)
}

#Para ejecutar ingrese los parametros como se muestra
##################################################
#MAIN()
#string1,string2,puntajeMiss,puntajeMatch
smithWaterman("AGCACACA","ACACACTA",1,2)

##################################################