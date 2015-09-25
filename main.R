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

validaStrings<-function(salidaDeReubicaStrings){
  if (typeof(salidaDeReubicaStrings)=="logical"){
    return (FALSE)
  }
  else{
    return (TRUE)
  }
}

#Parametros de entrada
#Se asigna string entrada 1
stringEntrada1<-"ATTAAG"
#Se asigna string entrada 2
stringEntrada2<-"GATTACAAG"
#Bonos y Penalizaciones
penalizacionIzquierda<-1
penalizacionArriba<-1
penalizacionCruzada<-1
penalizacionMissMatch<-1
bonoMatch<-1

#----------------------

entrada<-reubicaStrings(stringEntrada1,stringEntrada2);
if (validaStrings(entrada)==TRUE){
  #Genero matriz
  matriz=inicializaMatriz(entrada$parametro2$largo,entrada$parametro1$largo)
  #print (matriz)
  ejeX<-0
  ejeY<-0
  
  for (posLetra2 in 1:entrada$parametro2$largo){
    ejeX<-ejeX+1
    ejeY<-0
    for (posLetra1 in 1:entrada$parametro1$largo){
      ejeY<-ejeY+1
      #Valida correcto cruce de informacion para la comparacion
      #print (paste("Pos (",toString(ejeX),",",toString(ejeY),")",": Cruce ->",entrada$parametro1$valor[posLetra1],"-",entrada$parametro2$valor[posLetra2]))
    
    #-Formacion Matriz----------------------------------------------------------------# 
      if (ejeX==1){ 
        
        if (ejeY==1){ #ejeX==1 && ejeY==1
          
          if(entrada$parametro1$valor[posLetra1]==entrada$parametro2$valor[posLetra2]){ #casoMatch
            matriz[ejeX,ejeY]<-matriz[ejeX,ejeY]+bonoMatch
          }
          
          else{# casoMissMatch
            matriz[ejeX,ejeY]<-0
          }
          
        }#Fin caso #ejeX==1 && ejeY==1
        
        else{ #ejeX==1 && ejeY!=1

          if(entrada$parametro1$valor[posLetra1]==entrada$parametro2$valor[posLetra2]){ #casoMatch
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
        
        if(entrada$parametro1$valor[posLetra1]==entrada$parametro2$valor[posLetra2]){ #casoMatch
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
        
        mejorValorVecinos=max(matriz[ejeX,ejeY],matriz[ejeX-1,ejeY],matriz[ejeX,ejeY-1])
        
        if(entrada$parametro1$valor[posLetra1]==entrada$parametro2$valor[posLetra2]){ #casoMatch
          matriz[ejeX,ejeY]<-mejorValorVecinos+bonoMatch
        }
        
        else{# casoMissMatch
          if (mejorValorVecinos>penalizacionMissMatch){
            matriz[ejeX,ejeY]<-mejorValorVecinos-penalizacionMissMatch
          }
          else{
            matriz[ejeX,ejeY]<-0
          }
        }
        
      }#Fin ejeX!=1 && ejeY!=1
    
    #---------------------------------------------------------------------------------#
      
    } #Fin for (posLetra1 in 1:entrada$parametro1$largo)
  }#Fin for (posLetra2 in 1:entrada$parametro2$largo)
  
  print (matriz)
} else{
  print ("Uno de los parametros no es de tipo String y no se puede ejecutar")
}

warnings()
options(error=recover)