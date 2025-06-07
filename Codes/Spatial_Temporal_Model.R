
################################################################################
# Fecha de edición: 5 de junio del 2025
# Autor: Jorge Emmanuel Rosales Silva
# 
# En este código se presenta un camino para simular de la densidad que se propone 
# en el modelo espacial y temporal. Se usan datos de casos registrados de COVID-19 en la CDMX. 
# Este código hace con instrucciones "for" un muestreo de Gibbs que, además, emplea 
# Metropolis-Hastings para simular de cada condicional del algoritmo de Gibbs.
# Luego, para los hiperparámetros se usa la Gamma para generar una inversa Gamma
# Es decir, este cófigo es para hacer el "Metropolis-within-Gibbs con inversas Gamma"
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
  

  # MATRIZ DE VECINOS
  
  # Esta matriz se define con la intención de optimizar el código.
  # Definimos W como una matriz de ceros 
  W <- matrix(data=0, nrow=16, ncol = 16)
  
  # La llenamos con "1" con base en los vectors "uivecinos"
  # Es decir, esta matriz tiene un 1 en la entrada (i,j), si i es vecina de j
  # 0 en otro caso
  
  for(i in 1:16){
    
    for(j in vecinos[i]){
      W[i,j] <- 1
    }
    
  }
  
  col_names_w <- c(1:16)
  colnames(W) <- col_names_w
  W <-as.data.frame(W)
  
  
  
}



################################################################################
#
#                     3.      Función que hace Gibbs
#
################################################################################
# A continuación, se define la función que hace Metropolis-within-Gibbs con inversas Gamma
# Antes, se defie una matriz "C" que contiene los cálculos de cit.
# Esta se define de manera global, pues es empleada a lo largo del código, y no depende de 
# alguna iteració para su cálculo.
# Como la matriz C se define de manera global, es un parámetro dado a la función.
  


#Definimos el vector inicial 
theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,1,1,1,1)

# Matriz C

C<-matrix(0,16,12)
C<- as.data.frame(C)
for(t in 1:12){
  
  suma_casos_i <- data_wide[ , 2+t]
  suma_casos_i <- as.data.frame(suma_casos_i)
  suma_casos_i <- sum(suma_casos_i)
  
  for(i in 1:16){
    
    C[i,t] <- suma_casos_i*(poblacion_CENSO_2020[i]/poblacion_Total_CENSO_2020)
    
  }
}


BYM_AR<- function(data_wide,tamaño, theta, C,W){
  
  #Defino una matriz en la que guardaré los resultados del muestreo gibbs
  tamaño<- tamaño
  muestras<-matrix(0,tamaño,60) #Será una matriz que hace "tamaño" simulaciones del vector theta
  muestras<-as.data.frame(muestras) # Es un data frame de "tamaño" renglones y 34 columnas. El data frame comienza con ceros en cada entrada.
  colnames(muestras)<-c("alpha_1", "alpha_2", "alpha_3", "alpha_4", "alpha_5", "alpha_6", "alpha_7", "alpha_8", "alpha_9", "alpha_10", "alpha_11", "alpha_12",
                        "gamma_1", "gamma_2", "gamma_3", "gamma_4", "gamma_5", "gamma_6", "gamma_7", "gamma_8", "gamma_9", "gamma_10", "gamma_11", "gamma_12",
                        "u1","u2","u3","u4","u5","u6","u7","u8","u9","u10","u11","u12","u13","u14","u15","u16",
                        "v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15","v16",
                        "lambda_alpha","lambda_gamma","lambda_u","lambda_v")
  
  
  # A continuación, haré un for para tener una muestra de tamaño "tamaño"
  

  for(k in 1:tamaño){ 
    
    
    for(i in 1:60){                           #Recorremos las 60 entrada del vector
      
      

      # Primero hacemos el caso alpha_1
      if(i ==1){
        
        # A continuación se almacenas o calculan ciertas constantes necesarias para esta densidad
        
        # Si i vale 1, estamos en alpha_{1}.
        # Para alphga_{1} se usa la primer columna de C
        # La tomamos 
        
        cit<-C[,i] 

        #Calculamos la pi
  
        pi <- c()
        
        #acá sólo almaceno los sumandos 
        for (j in 1:16) {
          pi[j]<-cit[j]*exp(theta[12+i]+ theta[24+j] + theta[40 + j])
        }
        #acá ya lo sumo 
        pi<-sum(pi)
        
        # Defino la qi
        
        qi <- sum(cit)
        
        
        # Definimos el exponente de pi_tilde
        
        h_alpha_t <- function(w){
          return(  qi*w - pi*exp(w) - (1/(2*theta[57]))*(theta[2] - w)^2   ) 
        }
        
        
      # ************************************************************************
      # 1. Simulamos "y" con la propuesta g
      # ************************************************************************  
      
      # Esta g tiene como media el valor inicial (anterior) de theta[i]
        
        b<- (0.8/6)^2
        y <- rnorm(1,mean=theta[i], sd=sqrt(b))
        
      # ************************************************************************
      # 2. Calculo de r
      # ************************************************************************  
        
        log_pi_y  <-h_alpha_t(y)
        log_q_y_x <-dnorm(theta[i],mean = y,sd=sqrt(b),log = TRUE)
        
        log_pi_x  <-h_alpha_t(theta[i])
        log_q_x_y <-dnorm(y,mean = theta[i],sd=sqrt(b),log = TRUE)
        
        log_r     <- log_pi_y + log_q_y_x - log_pi_x - log_q_x_y
        
      # ************************************************************************
      # 3. Simular y asignar u
      # ************************************************************************
        
        log_u <- log(runif(1,0,1))
        
      # ************************************************************************
      # 4. Comparar y actualizar 
      # ************************************************************************
        
        theta[i] <- ifelse( log_u <= log_r, y, theta[i]  )

        
      }
        

      ############################################
      
      
      # Ahora, haré los casos de i en {2, ..., 11} (sigo con las alpha_i)
      
      if(i %in% 2:11){
        
        
        # Si i vale t, estamos en alpha_{t}.
        cit<-C[,i] 
        # Es decir, "alpha_{t}" con el "mes t".
        
        # Calculamos la pi
        pi <- c()
        #acá sólo almaceno los sumandos 
        for (j in 1:16) {
          pi[j]<-cit[j]*exp(theta[12+i]+ theta[24+j] + theta[40 + j])
        }
        #acá ya lo sumo 
        pi<-sum(pi)
        
        # Defino la qi
        
        qi <- sum(cit)
        
        
        
        # Definimos el exponente de pi_tilde
        
        h_alpha_t <- function(w){
          return(  qi*w - pi*exp(w) -  (1/(2*theta[57]))*(  ( w - theta[i-1]  )^2  + ( theta[i+1] -w  )^2   )       ) 
        }
        
        # ************************************************************************
        # 1. Simulamos "y" con la propuesta g
        # ************************************************************************  
        
        # Esta g tiene como media el valor inicial (anterior) de theta[i]
        
        b<- (0.8/6)^2
        y <- rnorm(1,mean=theta[i], sd=sqrt(b))
        
        # ************************************************************************
        # 2. Calculo de r
        # ************************************************************************  
        
        log_pi_y  <-h_alpha_t(y)
        log_q_y_x <-dnorm(theta[i],mean = y,sd=sqrt(b),log = TRUE)
        
        log_pi_x  <-h_alpha_t(theta[i])
        log_q_x_y <-dnorm(y,mean = theta[i],sd=sqrt(b),log = TRUE)
        
        log_r     <- log_pi_y + log_q_y_x - log_pi_x - log_q_x_y
        
        # ************************************************************************
        # 3. Simular y asignar u
        # ************************************************************************
        
        log_u <- log(runif(1,0,1))
        
        # ************************************************************************
        # 4. Comparar y actualizar 
        # ************************************************************************
        
        theta[i] <- ifelse( log_u <= log_r, y, theta[i]  )
        
        
      } 
      
      
      if(i ==12){
        
        
        # Si i vale 12, estamos en alpha_{12}.
        # Para alphga_{12} se usa la 12-ésima columna de C
        # La tomamos 
        
        cit<-C[,i] 
        
        #Calculamos la pi
        
        pi <- c()
        #acá sólo almaceno los sumandos 
        for (j in 1:16) {
          pi[j]<-cit[j]*exp(theta[12+i]+ theta[24+j] + theta[40 + j])
        }
        #acá ya lo sumo 
        pi<-sum(pi)
        
        # Defino la qi
        
        qi <- sum(cit)
        
        
        # Definimos el exponente de pi_tilde
        
        h_alpha_t <- function(w){
          return(  qi*w - pi*exp(w) - (1/(2*theta[57]))*(w - theta[i-1] )^2   ) 
        }
        
        
        # ************************************************************************
        # 1. Simulamos "y" con la propuesta g
        # ************************************************************************  
        
        # Esta g tiene como media el valor inicial (anterior) de theta[i]
        
        b<- (0.8/6)^2
        y <- rnorm(1,mean=theta[i], sd=sqrt(b))
        
        # ************************************************************************
        # 2. Calculo de r
        # ************************************************************************  
        
        log_pi_y  <-h_alpha_t(y)
        log_q_y_x <-dnorm(theta[i],mean = y,sd=sqrt(b),log = TRUE)
        
        log_pi_x  <-h_alpha_t(theta[i])
        log_q_x_y <-dnorm(y,mean = theta[i],sd=sqrt(b),log = TRUE)
        
        log_r     <- log_pi_y + log_q_y_x - log_pi_x - log_q_x_y
        
        # ************************************************************************
        # 3. Simular y asignar u
        # ************************************************************************
        
        log_u <- log(runif(1,0,1))
        
        # ************************************************************************
        # 4. Comparar y actualizar 
        # ************************************************************************
        
        theta[i] <- ifelse( log_u <= log_r, y, theta[i]  )
        
      }
      
      
      ############################################
      # Ahora, sigue actualizar las gamma 
      ############################################
      
      if(i %in% 13:24){
        
        # Para gamma_{t} se usa la t-ésima columna de C
        # La tomamos 
        
        cit<-C[,i-12] 
        # Es decir, "gamma_{t}" con el "mes t".
        
        pi <- c()
        #acá sólo almaceno los sumandos 
        for (j in 1:16) {
          pi[j]<-cit[j]*exp(theta[i - 12]+ theta[24+j] + theta[40 + j])   # En la pt de las gamma, se usan las alphas, por eso ahora el subíndice es "i-12"
        }
        #acá ya lo sumo 
        pi<-sum(pi)
        
        # Defino la qi
        
        qi <- sum(cit)
        
        # Definimos el exponente de pi_tilde
        
        h_gamma_t <- function(w){
          return(  qi*w - pi*exp(w)-  (1/(2*theta[58]))*w^2   ) 
        }
        
        # ************************************************************************
        # 1. Simulamos "y" con la propuesta g
        # ************************************************************************  
        
        # Esta g tiene como media el valor inicial (anterior) de theta[i]
        
        b<- (0.8/6)^2
        y <- rnorm(1,mean=theta[i], sd=sqrt(b))
        
        # ************************************************************************
        # 2. Calculo de r
        # ************************************************************************  
        
        log_pi_y  <-h_gamma_t(y)
        log_q_y_x <-dnorm(theta[i],mean = y,sd=sqrt(b),log = TRUE)
        
        log_pi_x  <-h_gamma_t(theta[i])
        log_q_x_y <-dnorm(y,mean = theta[i],sd=sqrt(b),log = TRUE)
        
        log_r     <- log_pi_y + log_q_y_x - log_pi_x - log_q_x_y
        
        # ************************************************************************
        # 3. Simular y asignar u
        # ************************************************************************
        
        log_u <- log(runif(1,0,1))
        
        # ************************************************************************
        # 4. Comparar y actualizar 
        # ************************************************************************
        
        theta[i] <- ifelse( log_u <= log_r, y, theta[i]  )

        
      }
      

        ############################################
        # Ahora, sigue actualizar las u_i
        ############################################
      
        if(25<=i & i<=40){    
          
          # Estas variables son ESPACIALES
          # Entonces los datos de C no se leen por tiempo, si no por espacio. 
          # Dicho de otra forma, ahora nos traemos los renglones de C, 
          # pues cada renglón representa una delegación a lo largo del tiempo.

          # Calculamos pt
          pt<-c()
          for (t in  1:12) {
            
            pt[t]<-C[i-24 , t]*exp( theta[t] + theta[t+12] + theta[i+16]  ) # el primer sumando son las alpha // El segundo son las gamma (por eso se le suma 12, para desplazarnos a la posición de las gamma)
                                                                # el tercer sumando va desplazado 16, para tomarnos las "v" en lugar de las "u"
            }
          
          pt <- sum(pt) # Listo crack !!!!!
          
          # Defino la qt
          
          qt <- c()
          for(t in 1:12){
            
            qt[t]<-data_wide[i-24,2+t]
            
          }
          
          qt<-as.data.frame(qt)
          qt<-sum(qt)
          
          # Creamos un vector que almacenará la resta de los cuadrados para los vecinos de cada ui
          sumandos_vecinos_ui <-c()
          for (contador in 1:16) {
            sumandos_vecinos_ui[contador] <- (theta[i]- theta[contador + 24])*W[i-24,contador]
            # A cada ui, se le restan TODAS las otras uj, incluso ui. 
            # Pero por construcción de W, sólo son distinto a cero aquellas que son vecinas de ui. 
            # Se hace la resta de las ui con todas las uj
            # Pero sólo sobrevive el resultado de las ui menos sus uj vecinas 
            
          }
          
          sumandos_vecinos_ui <- sumandos_vecinos_ui^2
          sumandos_vecinos_ui <- sum(sumandos_vecinos_ui)
          
          # Con esto, ya tenemos los elementos necesarios para calcular el exponente de la función de densidad de P(ui|...)
        
          
          h_u_t <- function(w){
            
            return(     -pt*exp(w) + qt*w - (1/(2*theta[59]))*sumandos_vecinos_ui       )
           
          }
          
          
          # ************************************************************************
          # 1. Simulamos "y" con la propuesta g
          # ************************************************************************  
          
          # Esta g tiene como media el valor inicial (anterior) de theta[i]
          
          b<- (0.5/6)^2
          y <- rnorm(1,mean=theta[i], sd=sqrt(b))
          
          # ************************************************************************
          # 2. Calculo de r
          # ************************************************************************  
          
          log_pi_y  <-h_u_t(y)
          log_q_y_x <-dnorm(theta[i],mean = y,sd=sqrt(b),log = TRUE)
          
          log_pi_x  <-h_u_t(theta[i])
          log_q_x_y <-dnorm(y,mean = theta[i],sd=sqrt(b),log = TRUE)
          
          log_r     <- log_pi_y + log_q_y_x - log_pi_x - log_q_x_y
          
          # ************************************************************************
          # 3. Simular y asignar u
          # ************************************************************************
          
          log_u <- log(runif(1,0,1))
          
          # ************************************************************************
          # 4. Comparar y actualizar 
          # ************************************************************************
          
          theta[i] <- ifelse( log_u <= log_r, y, theta[i]  )
          
        }                                                                         
        

      
        ############################################
        #Ahora, sigue actualizar las v_i
        ############################################
        
        if(41 <= i & i<=56){
          
          
          # Calculamos pt
          pt<-c()
          for (t in  1:12) {
            
            pt[t]<-C[i-40 , t]*exp( theta[t] + theta[t+12] + theta[i-16]  ) # el primer sumando son las alpha // El segundo son las gamma (por eso se le suma 12, para desplazarnos a la posición de las gamma)
            # el tercer sumando va desplazado -16, para tomarnos las "u" en lugar de las "v"
          }
          
          pt <- sum(pt) 
          
          # Defino la qt
          
          qt <- c()
          for(t in 1:12){
            
            qt[t]<-data_wide[i-40,2+t]
            
          }
          
          qt<-as.data.frame(qt)
          qt<-sum(qt)
          
          
          # Calculamos el exponente de la exponenicla en pi tilde
          
          h_v <- function(w){
            
            return( -pt*exp(w) + qt*w - (1/(2*theta[60]))*w^2  )
          }
         
          
          # ************************************************************************
          # 1. Simulamos "y" con la propuesta g
          # ************************************************************************  
          
          # Esta g tiene como media el valor inicial (anterior) de theta[i]
          
          b<- (0.5/6)^2
          y <- rnorm(1,mean=theta[i], sd=sqrt(b))
          
          # ************************************************************************
          # 2. Calculo de r
          # ************************************************************************  
          
          log_pi_y  <-h_v(y)
          log_q_y_x <-dnorm(theta[i],mean = y,sd=sqrt(b),log = TRUE)
          
          log_pi_x  <-h_v(theta[i])
          log_q_x_y <-dnorm(y,mean = theta[i],sd=sqrt(b),log = TRUE)
          
          log_r     <- log_pi_y + log_q_y_x - log_pi_x - log_q_x_y
          
          # ************************************************************************
          # 3. Simular y asignar u
          # ************************************************************************
          
          log_u <- log(runif(1,0,1))
          
          # ************************************************************************
          # 4. Comparar y actualizar 
          # ************************************************************************
          
          theta[i] <- ifelse( log_u <= log_r, y, theta[i]  )

          
        } 
        
        
        
        ############################################
        # Ahora,vamos con los 4 hiperparámetros 
        ############################################

        # lambda_alpha
        
        if(i==57){

          # Cálculo de alpha  [parámetro shape]
          alpha_shape<- (12/2) - 1
          
          # Cálculo de theta [parámetro scale]
          sum_alpha_t <- c()
          for(t in 2:12){
            
            sum_alpha_t[t-1] <- theta[t] - theta[t-1]
            
          }
          
          sum_alpha_t <- sum_alpha_t^2
          sum_alpha_t <- sum(sum_alpha_t)
          
          eps<- 0.1
          
          theta_scale <- 2/(sum_alpha_t + eps) 
          
          
          # Simular una IG con una Gamma de shape, scale
          q<-rgamma(1,shape = alpha_shape, scale = theta_scale)
          
          # Invertir el valor para tener una IG
          q<-1/q             
          
          # Actualizamos theta
          theta[57]<-q

        }
        
        
        ############################################
        # lambda_gamma
        
        if(i==58){
          
          # Cálculo de alpha  [parámetro shape]
          alpha_shape<- (12/2) - 1
          
          # Cálculo de theta [parámetro scale]
          
          sum_gamma_t <- c()
          
          for(t in 1:12){
            
            sum_gamma_t[t] <- theta[t+12]^2
            
          }
          
          sum_gamma_t <- sum(sum_gamma_t)
          
          eps<- 0.1
          
          theta_scale <- 2/(sum_gamma_t + eps) 
          
          
          # Simular una IG con una Gamma de shape, scale
          q<-rgamma(1,shape = alpha_shape, scale = theta_scale)
          
          # Invertir el valor para tener una IG
          q<-1/q             
          
          # Actualizamos theta
          
          theta[58]<-q
          
        }
        
        

        ############################################
        # lambda_u
        
        if(i==59){
          
          # Cálculo de alpha  [parámetro shape]
          alpha_shape<- (16/2) - 1
          
          # Cálculo de theta [parámetro scale]
          
          suma_Todos_vs_Todos <- 0
          
          for(contador_renglon in 1:16){
            
            for(contador_columna in contador_renglon:16){
              
              auxiliar <- theta[contador_renglon+24] - theta[contador_columna + 24] # se toma la resta de ui contra las uj a la derecha de la diagonal
              auxiliar <- auxiliar^2 # eleva la resta al cuadrado
              auxiliar <- auxiliar*W[contador_renglon , contador_columna] # si son vecinos, lo deja; si no, lo manda a cero.
              
              suma_Todos_vs_Todos <- suma_Todos_vs_Todos + auxiliar # almacena el valor en el acumulado
              
              auxiliar <- 0
              
            }
          }
          
          
          # Definición de épsilon 
          eps<- 0.01
          
          # Cálculo de beta
          theta_scale<- 2/(eps+suma_Todos_vs_Todos)

          
          # Simular una IG con una Gamma de shape, scale
          q<-rgamma(1,shape = alpha_shape, scale = theta_scale)
          
          # Invertir el valor para tener una IG
          q<-1/q             
          ################  METER ALGUNA CONDICIÓN EN CASO DE QUE q=0
          
          # Actualizamos theta
          
          theta[59]<-q
          
        }
        
        

        ############################################
        # lambda_v
        
        if(i==60){

          
          # Cálculo de alpha  [parámetro shape]
          
          alpha_shape<- (16/2) - 1
          
          # Cálculo de theta [parámetro scale]
          
          sum_v_t <- c()
          
          for(t in 1:16){
            
            sum_v_t[t] <- theta[40 + t]^2
            
          }
          
          sum_v_t <- sum(sum_v_t)
          
          # Definición de épsilon 
          eps<- 0.01
          
          # Cálculo de beta
          theta_scale<- 2/(eps + sum_v_t)
          
          
          # Simular una IG con una Gamma de shape, scale
          q<-rgamma(1,shape = alpha_shape, scale = theta_scale)
          
          # Invertir el valor para tener una IG
          q<-1/q             
          
          # Actualizamos theta
          
          theta[60]<-q
          
        }
        
   
        
  
          
        }
              # Acá cierra el for de tamaño 60 , es decir, acá dentro ocurre una iteración completa de Gibbs
    
  
    # Con el vector simulado por completo, se almacena en la matriz definida al inicio.
    muestras[k,]<-theta

    

      } # Acá cierra el for de tamaño "tamaño". Es decir, acá ya se llenó la matriz. 
        # Ya se hizo "tamaño" veces el for de tamaño 60.
      
  
  return(muestras)  # Ya se tiene como output un data frame con "tamaño" simulaciones del vector de tamaño 60 
  
  
  
  
  
}     





################################################################################
#
#          4. CASOS: Correr el código y generar 200 000 muestras
#      
################################################################################
rm(data)
rm(data_long)
rm(data_monthly)
options(java.parameters = "-Xmx8g")
library(xlsx)
# A continuación, haré 200K simulaciones para cada mes
# Se harán a partir de bloques de tamaño 10 000


# Se fija una semilla
set.seed(182435)

tiempo_ejecucion <- system.time({
  
  
  { 
     
      
      #Definimos el número de muestras que queremos simular
      tamaño <- 10000
      
      
      # Dentro de este if, que hará la simulación para el primer mes, hacemos un for para hacer 200K en muestras de tamaño 10K
      for(z in 1:20){
    
        
        print(paste("Vamos en la iteración ", z, " de 20"))

        if (z == 1) { # Primera bloque de 10k
          
          
          
          #theta inicial 
          theta=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,1,1,1,1)
          
          M1<-BYM_AR(data_wide,tamaño, theta,C,W)
          # Para la primera iteración, escribir con el encabezado
          write.table(M1, file = "X/Espacial_Temporal_200.csv", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
          
          
        } else { # Acá entra para z en {2,3,...,19,20}
          
          #Actualizamos el valor de "theta" al último en la iteración anterior
          theta <- M1[10000,]
          # Es necesario hacer este cambio para tener un objeto en el formato de "c()" y no de un data frame. Esto facilita el uso de este vector
          theta<- c(theta$alpha_1,theta$alpha_2,theta$alpha_3,theta$alpha_4,theta$alpha_5,theta$alpha_6,theta$alpha_7,theta$alpha_8,theta$alpha_9,theta$alpha_10,theta$alpha_11,theta$alpha_12
                    ,theta$gamma_1,theta$gamma_2,theta$gamma_3,theta$gamma_4,theta$gamma_5,theta$gamma_6,theta$gamma_7,theta$gamma_8,theta$gamma_9,theta$gamma_10,theta$gamma_11,theta$gamma_12
                          ,theta$u1,theta$u2,theta$u3,theta$u4,theta$u5,theta$u6,theta$u7,theta$u8,theta$u9,theta$u10,theta$u11,theta$u12,theta$u13,theta$u14,theta$u15,theta$u16
                          ,theta$v1,theta$v2,theta$v3,theta$v4,theta$v5,theta$v6,theta$v7,theta$v8,theta$v9,theta$v10,theta$v11,theta$v12,theta$v13,theta$v14,theta$v15,theta$v16
                          ,theta$lambda_alpha,theta$lambda_gamma,theta$lambda_u,theta$lambda_v)

          
          #Liberamos esa memoria, pues ya lo escribimos y ya tomamos el último valor
          rm(M1)
          
          #Volvemos a simular y a escribir los resultados
          M1<-BYM_AR(data_wide,tamaño, theta,C,W)
          # Ya no ponemos encabezado, pues no hace falta, ya que lo vamos a escribir en el archivo que ya tiene encabezado.
          write.table(M1, file = "X/Espacial_Temporal_200.csv", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
          
          
        }
        
        
        
      }
      
      rm(M1)
      
      
  }
  
  
})


