
################################################################################
# Fecha de edición: 4 de junio del 2025
# Autor: Jorge Emmanuel Rosales Silva
# 
# En este código se presenta un camino para simular de la densidad que se propone 
# en el modelo de BYM. Se usan datos de casos registrados de COVID-19 en la CDMX. 
# Este códifo hace con instrucciones "for" un muestreo de Gibbs que, además, emplea 
# aceptación rechazo para generar realizaciones de cada condicional del algoritmo 
# de Gibbs.
################################################################################



################################################################################
#                       1.    Lectura de datos 
################################################################################
{

# Dirección donde están los datos 
# Acá se establece la ruta "X" donde tiene los registros de contagio
setwd("X")
#install.packages("dplyr")
library(dplyr)
library(tidyr)

data <- read.csv("X/Casos_Diarios_Municipio_Confirmados_20230625.csv") # Se leen los datos de contagios de covid
# Estos datos fueron tomados de:
# https://datos.covid-19.conacyt.mx/#DownZCSV
data <- as.data.frame(data)                                             # Se les da formato de data frame

clave_del <- c(9002,9003,9004,9005,9006,9007,9008,9009,9010,9011,9012,9013,9014,9015,9016,9017) # Las claves de las 16 alcaldías de la CDMX

# Nos quedamos sólo con los renglones de esas alcaldías
data <- data %>%
  filter(cve_ent %in% clave_del ) %>%
  arrange(nombre)                          # Los ordenamos en orden alfabético, pues el grafo está enumerado de esa manera.

poblacion<- data$poblacion                 # almacenamos los datos de población, por si se llegan a ocupar. Aunque usaremos los del CENSO 2020

data <- data[,-c(1,2)]                     # eliminamos la columna de la clave y la de población 

######################################
# Dar formato a las fechas 
######################################
# Es necesario dar formato a las fechas, para que puedan ser intepretadas adecuadamente.
# Nos quedamos con los nombres de las columnas 

column_names <- names(data)

# Sólo nos quedamos con las columnas de fechas: quitamos la primer columna, pues a esa NO se le hará un tratamiento especial de formato fecha.

date_columns <- column_names[column_names != "nombre"]

# Le quitamos la "X" a estas columnas de fechas 

clean_date_names <- gsub("X", "", date_columns)

# Ahora, damos formato de fecha a estas 

formatted_date_names <- as.Date(clean_date_names, format = "%d.%m.%Y")


# Ahora que tienen formato de fecha, las volvemos String para que puedan ser leídas como encabezado de columna
formatted_date_names <- format(formatted_date_names, "%d.%m.%Y")

# Asignamos estas "fechas-string" a los nombres de las columnas,con excepción de la primer columna. Estra columna la renombramos como "delegacion"
# Se trabajó con la etiqueta delegación, por una vieja costumbre. 
# Pero eso es irrelevante, nos referimos a las actuales alcaldías.
names(data) <- c("delegacion", formatted_date_names)


#######################################################################
# Sumar los datos por MES y por RENGLÓN (delegación)
######################################################################
# En este código se toman los registros agrupados/ sumados por mes. 
# En lugar de los registros diarios.

# Paso 1: Para todas las columnas, MENOS la de "delegacion", nos quedamos sólo con el año y con el mes.
data_long <- data %>%
  pivot_longer(cols = -delegacion, names_to = "date", values_to = "value") %>%
  mutate(month = format(as.Date(date, format = "%d.%m.%Y"), "%Y-%m"))


# Paso 2: Agrupamos por delegacion y por mes, para obtener la suma de cada delegación a lo largo de un mes.
data_monthly <- data_long %>%
  group_by(delegacion, month) %>%
  summarise(monthly_sum = sum(value, na.rm = TRUE)) %>%
  ungroup()

# Paso 3: Volvemos al formato de las delegaciones como renglón, pero ahora en las columnas está el acumulado por mes.
data_wide <- data_monthly %>%
  pivot_wider(names_from = month, values_from = monthly_sum)

################################################################################
# Población por delegación con base al CENSO 2020
################################################################################
# Como se comentó, se tomará como población al mes "t" en a delegación "i"
# la población del CENSO 2020.

# Para este vector, se toma que la entrada "i-ésima" se refiere al nodo "i-ésimo" como se definió en el GMRF
poblacion_CENSO_2020 <- c(759137, 432205, 434153, 614447, 217686, 545884, 1173351, 404695, 1835486, 247622, 414470, 152685, 392313, 699928, 443704, 442178  )
# Estas cifras fueron tomadas de:
# https://www.sedeco.cdmx.gob.mx/storage/app/media/uploaded-files/resultados-del-censo-pob-y-viv-2020-1.pdf
# 

poblacion_Total_CENSO_2020 <- sum(poblacion_CENSO_2020)



}


################################################################################
#                    2.   Definición de parámetros fijos
################################################################################
# Para hacer simulación del modelo BYM, es necesario tener codificadas las relaciones 
# de qué delegación es vecina de cuál delegación, para así poder calcular ciertos 
# términos. 


{
  # Se define qué vecinos tiene cada "i".
  u1vecinos<-c(3,11,5,10,4,14)
  u2vecinos<-c(11,6,7)
  u3vecinos<-c(11,6,8,9,4,1)
  u4vecinos<-c(1,3,9,16,14)
  u5vecinos<-c(11,1)
  u6vecinos<-c(2,7,15,8,3,11)
  u7vecinos<-c(2,6,15)
  u8vecinos<-c(15,6,3,9)
  u9vecinos<-c(8,3,4,16,13)
  u10vecinos<-c(1,14)
  u11vecinos<-c(2,6,3,1,5)
  u12vecinos<-c(14,16,13)
  u13vecinos<-c(9,16,12)
  u14vecinos<-c(10,1,4,16,12)
  u15vecinos<-c(7,6,8)
  u16vecinos<-c(4,9,14,13,12)

  
  #Para poder recorrer en un cilo la lista de vecinos de cada i-ésima región, estos se almacenan en una lista.
  vecinos<-list(u1vecinos,u2vecinos,u3vecinos,u4vecinos,u5vecinos,u6vecinos,u7vecinos,u8vecinos,u9vecinos,u10vecinos,u11vecinos,u12vecinos,u13vecinos,u14vecinos,u15vecinos,u16vecinos)
  
  # cardinalidad del conjunto de vecinos para cada región i
  ni<-c(length(u1vecinos),length(u2vecinos),length(u3vecinos),length(u4vecinos),length(u5vecinos),length(u6vecinos),length(u7vecinos),length(u8vecinos),length(u9vecinos),length(u10vecinos),length(u11vecinos),length(u12vecinos),length(u13vecinos),length(u14vecinos),length(u15vecinos),length(u16vecinos))
  
  # Se da un valor inicial a épsilon.
  eps<-0.01
  
  
  s1<-c(1,1,1,1,1,1,2,2,2,3,3,3,3,3,4,4,4,5,6,6,6,6,7,8,8,9,9,10,12,12,12,13,14)  #Cada entrada de este vector representa la "i-ésima zona": Las primeras 6 entradas hacen referencia a la "zona 1"
  s2<-c(3,4,5,10,11,14,6,7,11,4,6,8,9,11,9,14,16,11,7,8,11,15,15,9,15,13,16,14,13,14,16,16,16) # La j-ésima entrada de este vector, nos dice que eso es vecino de lo que está en la j-ésima entrada del vector de arriba (s1)

  # Es decir, los pares ordenados de la forma (s1[i],s2[i]) nos da los 33 posibilidades de "vecinos", sin repetición. Lo que sería igual a tener la diagonal superior, o inferior, de la matriz de precisión.
  
}



################################################################################
#
#                     3.      Función que hace Gibbs
#
################################################################################
#Se define una función que hace Gibbs con AR, como se menciona en el paper BYM.
# Esta función ocupa TRES parámetros: el vector de casos observados "y"
# que en nuestro caso estará dado por las entradas dada una fecha y una delegación del data frame "data_wide".
# Y, además, ocupa el tamaño de las muestras que se quiere calcular.
# Finalmente, requiere de un vector "theta" como input inicial a ser actualizado.

BYM_AR<- function(y,tamaño, theta){
  
  #Defino una matriz en la que guardaré los resultados del muestreo de Gibbs
  tamaño<- tamaño
  muestras<-matrix(0,tamaño,34) #Será una matriz que hace "tamaño" simulaciones del vector theta
  muestras<-as.data.frame(muestras) # Es un data frame de "tamaño" renglones y 34 columnas. El data frame comienza con ceros en cada entrada.
  # Se da nombre a las columnas del data frame que se generará
  colnames(muestras)<-c("u1","u2","u3","u4","u5","u6","u7","u8","u9","u10","u11","u12","u13","u14","u15","u16",
                        "v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15","v16",
                        "lambda_u","lambda_v")
                        

  # A continuación, haré un for para tener una muestra de tamaño "tamaño"
  
  # PARA CADA RENGLÓN del data frame se hará Gibbs sampling
  # Cada columna es una variable a actualizar
  for(k in 1:tamaño){ 

    c<-c() # Defino este vector para almacenar posteriormente el valor de "ci"
    
    for(i in 1:34){                           #Recorremos las 34 entrada del vector
      
      # Para cada "i" tendrémos una ci, la cual calculamos a continuación.
      # Todos los números involucrados son constantes dado un periodo de tiempo, PERO para cada "i"
      # cambia la suma de los casos (sum_ys)  y también cambia la población de la region i.
      # Es decir, se tomarán distintas entradas del vector "y" con base a la alcaldía en turno.
      
      
      #Primero, se actualizan las u1,...,u16
      if(1<=i & i<=16){        
        
        #Calculamos el promedio de las vecinas para esa ui (es un témrino que se ocupa en las distribuciones a simular)
        
        
        vecinosui<-vecinos[i]                         # Saber cuáles son los vecinos de esa ui
        vecinosui<-as.data.frame(vecinosui)           # Quedarme sólo con los elementos de la lista que necesito, lo convierto de lista a data frame
        vecinosui<-t(vecinosui)                       # Lo transpongo (igual por cuestiones de formato)
        
        uVecinasPromedio=0                            # Incializo la variable para calcular ese promedio
        

        for(j in vecinosui){
          uVecinasPromedio= uVecinasPromedio+theta[j] # voy sumando los valores de u_{j} para las vecinas de u_{i} 
        }
        uVecinasPromedio=uVecinasPromedio/ni[i]       # Finalmente, para tener el promedio, divido entre el número de vecinos.
        
        #Definimos eps
        eps=0.1
        
        #Calculamos lambda_u (ya está en el vector inicial/actual)
        lambda_u<-theta[33]
        
        #Calculamos yi
        yi<-y[i]                                    # Estos valores son los registros (suma mensual) de contagios de covid en la zona i-ésima.
      
        
        # Cálculo ci
        sum_ys <- sum(y)
        c[i]<- (sum_ys)*(poblacion_CENSO_2020[i]/poblacion_Total_CENSO_2020) 
        
        #####################################
        #Cálculo de a (media de la propuesta)
        #####################################
        # Con la función BFGS optim()
        # Para que la propuesta de aceptación rechazo se parezca a la objetivo
        # se calcula dónde llega a su máximo, para centrar la propuesta en ese valor.
        
        {
          
          hu <- function(u) {
            
              -c[i] * exp(u) * exp(theta[16 + i]) +
                u * yi -
                (ni[i] / (2 * lambda_u)) * u^2 +
                (ni[i] / lambda_u) * uVecinasPromedio * u
            
          }
          
          # Maximización con optim usando BFGS
          au_result <- optim(par = theta[i], fn = hu, method = "BFGS", control = list(fnscale = -1))
          au <- au_result$par
          
        }
        
        #Cálculo de b^2
        
        b2_u<- 2*lambda_u/ni[i]
        
        #Definición de g
        log_gu<-function(u){ return( -(1/(2*b2_u))*(u^2) + (au/(b2_u))*u ) }
        
        #####################################
        # Cálculo de M
        #####################################
        # Con la función BFGS optim(

        {
          
          log_ratio_pi_g_neg <- function(u) {
            hu(u) - log_gu(u)
          }
          
          optim_ratio_result <- optim(par = theta[i], fn = log_ratio_pi_g_neg, method = "BFGS", control = list(fnscale = -1))
          log_Mu <- optim_ratio_result$value
          
        }

        
        #####################################
        # Aceptación Rechazo para ui
        #####################################
        
        #Hacemos un "do while", pra hacerlo tantas veces se requiera para conseguir un ui
        
        repeat{
          x<-rnorm(1,mean=au,sd=sqrt(b2_u))     # 1.  Simulación de una variable de la densidad g 
          
          # Definimos pi tilde                  # 2. Cálculo de log_r
          vi<-theta[16+i]
          log_pix<- hu(x)

          #Defino g
          log_gx<- log_gu(x)
          
          # Defino el log de r
          log_r<- log_pix - log_gx - log_Mu

                                               # 3. Simular y asignar unif(0,1)
          u<-runif(1,0,1)  #Generamos y asignamos la variale aleatoria para comparar con la r

                                               # 4. Comparar u con r
          if(log(u)<=log_r){
            
            break
          }
          
        }                                     # FIN de AR para ui
        
        # Si ya encontramos una "x" tal que se cumple AR, nos salimos del Do while
        
        #Actualizar theta
        #Con esta x , actualizamos ui en el vector theta 
        
        theta[i]<-x
        
      }                                      # Aquí termina la simulación de las ui
      
    
      
      
      # Ahora, hacemos algo análogo para vi
      
      if(17<=i & i<=32){

        
        #Comenzamos por definir las variables que se necesitan para los cálculos de aceptación rechazo
        
        # Cálculo de lambda_v
        lambda_v<-theta[34]
        
        # Cálculo de yi
        yi<-y[i-16]     #pues a estas alturas la "i" ya pasó por las 16 primeras
        
        # Cálculo de epsilon
        eps=0.1
        

        
        # Cálculo de a
        # Con la función optim

        {
          # Se define el exponente 
          hv <- function(v) {
            
              -c[i - 16] * exp(v) * exp(theta[i - 16]) +
                v * yi -
                (1 / (2 * lambda_v)) * v^2
            
          }
          
          av_result <- optim(par = theta[i], fn = hv, method = "BFGS", control = list(fnscale = -1))
          av <- av_result$par

          
        }
        

        
        # Cálculo b2
        
        b2_v <- 2*lambda_v
        
        #Defino g
        log_gv<-function(v){ return( -(1/(2*b2_v))*(v^2) + (av/(b2_v))*v ) }
        
        # Cálculo de M
        # Con la función optim

        {
          
          log_ratio_pi_g_v_neg <- function(v) {
            hv(v) - log_gv(v)
          }
          
          optim_ratio_v_result <- optim(par = theta[i], fn = log_ratio_pi_g_v_neg, method = "BFGS", control = list(fnscale = -1))
          log_M_v <- optim_ratio_v_result$value
          
        }
        

        
        # Hacemos AR para la vi
        
        repeat{
          x<-rnorm(1,mean=av,sd=sqrt(b2_v))     # 1.  Simulación de una variable de la densidad g 
          
          
          
          # Definimos la función "pi tilde"     # 2. Cálculo de log_r
          
          log_pix_v<- hv(x)
          
          
          # Evalúo la g 
          log_gx_v<- log_gv(x)
          # 2. Calculo de r(x) 
          log_r_v<- log_pix_v - log_M_v - log_gx_v
          
                                               # 3. Simular y asignar unif(0,1)
          
          u<-runif(1,0,1)  #Generamos y asignamos la variale aleatoria para comparar con la r
          
                                               # 4. Comparar u con r
          if(log(u)<=log_r_v){
            break
          }
          
        }
        # Si ya encontramos una "x" tal que se cumple AR, nos salimos del Do while
        
        #Actualizar theta
        #Con esta x , actualizamos vi en el vector theta 
        
        theta[i]<-x

      }
      
      
      #Ahora, vamos a actualizar lambda_{u}
      
      
      if(i==33){
        
        #Cálculo de alpha
        alpha<- (16/2) - 1
        
        #Cálculo de la suma de los cuadrados de las diferencias entre TODAS las u's vecinas
        sumaT <- theta[s1]- theta[s2]
        sumaT <- sumaT^2
        sumaT <- sum(sumaT)
        
        # Definición de épsilon 
        eps<- 0.01
        
        # Cálculo de theta (scale)
        scale<- 2/(eps+sumaT)
        
        # Simular una IG con una Gamma de shape, scale
        q<-rgamma(1,shape = alpha, scale = scale)
        
        # Invertir el valor para tener una IG
        q<-1/q             
        
        # Actualizamos theta
        
        theta[33]<-q

        
      }
      
      
      # Finalmente, vamos a actualizar lambda_{v}
      
      if(i==34){
        
        #Cálculo de alpha
        alpha<- (16/2) - 1
        
        #Cálculo de la suma de los cuadrados de vi
        sumav<-theta[17:32]
        sumav<-sum(sumav^2)
        
        # Definición de épsilon 
        eps<- 0.01
        
        # Cálculo de theta(scale) 
        scale<- 2/(eps+sumav)
        
        # Simular una IG con una Gamma de shape, scale
        q<-rgamma(1,shape = alpha, scale = scale)
        
        # Invertir el valor para tener una IG
        q<-1/q             
        
        # Actualizamos theta
        
        theta[34]<-q

      }
      
      
    }
    
    # Con el vector theta totalmente actualizado/simulado, se guarda en la matriz definida al inicio.
    muestras[k,]<-theta

    print(paste("pasamos a la iteración", k))  # OPCIONAL: para tener una idea de cómo está avanzando el código.
    
    
    
    
  }
  
  
  
  return(muestras)
  
  
  
  
}


################################################################################
#
#                  4.   CASOS: Correr el código para 12 meses
#        
################################################################################

library(xlsx)

# A continuación, se define un bloque de código para cada mes para el que se simulará.
# Esta estructura por bloques consiste en lo siguiente:
# Hay un for para hacer las simulaciones para cada mes de marzo 2020 a febrero 2021.
# Y, además, dentro del bloque de cada mes hay un for que genera 100 000 simulaciones de 10 000 en 10 000.
# Esto se definió así por cuestiones de memoria, para que R no cargue con una matriz de 90 000 renglones, por ejemplo, 
# mientras calcula las últimas 10 000 simulaciones. 
# En su lugar, se va escribiendo el csv por bloques de 10 000, e inmediatamente después se eliminan de la memoria de R.


# Se fija una semilla.
set.seed(182435)

tiempo_ejecucion <- system.time({ # Para medir los tiempos necesarios.
  
  
  
  for(l in 1:12){ # Para recorrer los 12 meses.
    

    if(l ==1){ # Para 2020/03
      
      
      # Casos observados para este mes (Es la misma para las 10 iteraciones de este mes)
      y<- data_wide$`2020-03`
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000 # Los sub bloques antes mencionados
      
      
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){

        if (z == 1) { # Primera bloque de 10 k
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)

          
          M1<-BYM_AR(y,tamaño,theta)

          # Para la primera iteración, escribir con el encabezado
          write.table(M1, file = "X/Marzo_2020.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          # A este csv lo llamamos "Marzo_2020".
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M1[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilita el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                  ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                  ,theta$lambda_u,theta$lambda_v)
          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M1)
          
          #Volvemos a simular y a escribir los resultados
          M1<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M1, file = "X/Marzo_2020.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
      }
      
      rm(M1)
      
      
    }
    
    
    if(l ==2){ # Para 2020/04
      
      
      # Casos observados para este mes
      y<- data_wide$`2020-04`
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){
        
        if (z == 1) { # Primera iteración
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)
          
          M2<-BYM_AR(y,tamaño,theta)
        
          # Para la primera iteración, escribir con el encabezado
          write.table(M2, file = "X/Abril_2020.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M2[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilida el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                    ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                    ,theta$lambda_u,theta$lambda_v)
          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M2)
          
          #Volvemos a simular y a escribir los resultados
          M2<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M2, file = "X/Abril_2020.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
      }
      
      rm(M2)
      
      
    }
    
    
    if(l ==3){ # Para 2020/05
      
      
      # Casos observados para este mes
      y<- data_wide$`2020-05`
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){
        
        if (z == 1) { # Primera iteración
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)
          
          M3<-BYM_AR(y,tamaño,theta)
          
          # Para la primera iteración, escribir con el encabezado
          write.table(M3, file = "X/Mayo_2020.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M3[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilida el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                    ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                    ,theta$lambda_u,theta$lambda_v)
          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M3)
          
          #Volvemos a simular y a escribir los resultados
          M3<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M3, file = "X/Mayo_2020.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
      }
      
      rm(M3)
      
      
    }
    
    
    if(l ==4){ # Para 2020/06
      
      
      # Casos observados para este mes
      y<- data_wide$`2020-06`
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){)
        
        if (z == 1) { # Primera iteración
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)
          
          M4<-BYM_AR(y,tamaño,theta)
          
          # Para la primera iteración, escribir con el encabezado
          write.table(M4, file = "X/Junio_2020.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M4[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilida el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                    ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                    ,theta$lambda_u,theta$lambda_v)
          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M4)
          
          #Volvemos a simular y a escribir los resultados
          M4<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M4, file = "X/Junio_2020.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
      }
      
      rm(M4)
      
      
    }
    
    
    if(l ==5){ # Para 2020/07
      
      
      # Casos observados para este mes
      y<- data_wide$`2020-07`
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){
        
        if (z == 1) { # Primera iteración
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)
          
          M5<-BYM_AR(y,tamaño,theta)
          
          # Para la primera iteración, escribir con el encabezado
          write.table(M5, file = "X/Julio_2020.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M5[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilida el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                    ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                    ,theta$lambda_u,theta$lambda_v)
          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M5)
          
          #Volvemos a simular y a escribir los resultados
          M5<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M5, file = "X/Julio_2020.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
      }
      
      rm(M5)
      
      
    }
    
    
    if(l ==6){ # Para 2020/08
      
      # Casos observados para este mes
      y<- data_wide$`2020-08`
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){
        
        if (z == 1) { # Primera iteración
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)
          
          M6<-BYM_AR(y,tamaño,theta)
          
          # Para la primera iteración, escribir con el encabezado
          write.table(M6, file = "X/Agosto_2020.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M6[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilida el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                    ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                    ,theta$lambda_u,theta$lambda_v)
          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M6)
          
          #Volvemos a simular y a escribir los resultados
          M6<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M6, file = "X/Agosto_2020.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
      }
      
      rm(M6)
      
      
    }
    
    
    if(l ==7){ # Para 2020/09
      
      # Casos observados para este mes
      y<- data_wide$`2020-09`
      
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){
        
        if (z == 1) { # Primera iteración
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)
          
          M7<-BYM_AR(y,tamaño,theta)
          
          # Para la primera iteración, escribir con el encabezado
          write.table(M7, file = "X/Septiembre_2020.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M7[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilida el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                    ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                    ,theta$lambda_u,theta$lambda_v)
          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M7)
          
          #Volvemos a simular y a escribir los resultados
          M7<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M7, file = "X/Septiembre_2020.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
         
          
        }
      }
      
      rm(M7)
      
      
    }
    
    
    if(l ==8){  # Para 2020/10
      
      
      # Casos observados para este mes
      y<- data_wide$`2020-10`
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){
        
        if (z == 1) { # Primera iteración
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)
          
          M8<-BYM_AR(y,tamaño,theta)
          
          # Para la primera iteración, escribir con el encabezado
          write.table(M8, file = "X/Octubre_2020.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M8[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilida el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                    ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                    ,theta$lambda_u,theta$lambda_v)
          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M8)
          
          #Volvemos a simular y a escribir los resultados
          M8<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M8, file = "X/Octubre_2020.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
      }
      
      rm(M8)
      
      
    }
    
    
    if(l ==9){ # Para 2020/11
      
      # Casos observados para este mes
      y<- data_wide$`2020-11`
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){
        
        if (z == 1) { # Primera iteración
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)
          
          M9<-BYM_AR(y,tamaño,theta)
          
          # Para la primera iteración, escribir con el encabezado
          write.table(M9, file = "X/Noviembre_2020.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M9[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilida el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                    ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                    ,theta$lambda_u,theta$lambda_v)
          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M9)
          
          #Volvemos a simular y a escribir los resultados
          M9<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M9, file = "X/Noviembre_2020.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
        
      }
      
      rm(M9)
      
      
    }
    
    
    if(l ==10){ # Para 2020/12
      
      # Casos observados para este mes
      y<- data_wide$`2020-12`
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){
        
        if (z == 1) { # Primera iteración
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)
          
          M10<-BYM_AR(y,tamaño,theta)
          
          # Para la primera iteración, escribir con el encabezado
          write.table(M10, file = "X/Diciembre_2020.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M10[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilida el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                    ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                    ,theta$lambda_u,theta$lambda_v)
          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M10)
          
          #Volvemos a simular y a escribir los resultados
          M10<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M10, file = "X/Diciembre_2020.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
      }
      
      rm(M10)
      
      
    }
    
    
    if(l ==11){ # Para 2021/01
      
      
      # Casos observados para este mes
      y<- data_wide$`2021-01`
      
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){
        
        if (z == 1) { # Primera iteración
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)
          
          M11<-BYM_AR(y,tamaño,theta)
          
          # Para la primera iteración, escribir con el encabezado
          write.table(M11, file = "X/Enero_2021.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M11[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilida el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                    ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                    ,theta$lambda_u,theta$lambda_v)
          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M11)
          
          #Volvemos a simular y a escribir los resultados
          M11<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M11, file = "X/Enero_2021.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
      }
      
      rm(M11)
      
      
    }
    
    
    if(l ==12){ # Para 2021/02
      
      # Casos observados para este mes
      y<- data_wide$`2021-02`
      
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 100K en muestras de tamaño 10K
      for(z in 1:10){
        
        if (z == 1) { # Primera iteración
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
                  1,1)
          
          M12<-BYM_AR(y,tamaño,theta)
          
          # Para la primera iteración, escribir con el encabezado
          write.table(M12, file = "X/Febrero_2021.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,10}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M12[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilida el uso de este vector
          theta<- c(theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                    ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                    ,theta$lambda_u,theta$lambda_v)
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M12)
          
          #Volvemos a simular y a escribir los resultados
          M12<-BYM_AR(y,tamaño,theta)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M12, file = "X/Febrero_2021.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
      }
      
      rm(M12)
      
    }
        
    
  }
  
})
