
################################################################################

# Autor: Emmanuel Rosales
# Fecha: 7 de junio de 2025
# En este código se generan 16 bloques de 500K simulaciones.
# Antes de correr cada bloque se inicializará la semilla en 182435
# El primer bloque se inicializa con un vector theta inicial, de valores dados por mí.
# Después de esto, los bloques posteriores se inicializarán con la última simulación del ciclo anterior 
# es decir, el renglón 500K de la matriz que almacena las simulaciones.

# A DIFERENCIA del código para el modelo espacial y temporal, en este caso 
# el Metropolis-within-gibbs es multivariado. 
# Es decir, las simulaciones se hacen por bloques y no variable a variable.
# Por ejemplo, para el bloque alpha la propuesta se da con una normal de dimensión 12.
# De modo que, se actualiza con el valor propuesto, o se queda igual, para las 12 alphas en un solo paso.

# Además, en este código se define la función objetivo casi completa (sin la parte de los hiperparámetros)
# Esto se hace para optimizar el código.
# Dada una variable, el cociente de las condicionales es lo mismo que evaluar la densidad completa con todas las otras variables sin cambiar. 
# Eso es lo que se hace, por lo tanto no es necxesario definir todas las condicionales. 
# Como los hiperparámetros no se simuln con MH, y ese bloque se simplifica a "1" en el cociente de MH, se decidió no agregar ese bloque. 
# Además, como esta simulaciones son computacionalmente intensivas, se trató de disminuir las operaciones tanto como fuera posible.

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
# Para hacer simulación del modelo BYM extendido, es necesario tener codificadas las relaciones 
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
  
  
  # Cálculo de la Matriz C
  
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
  
  
  
}




################################################################################
#
#             3.    CODIFICACIÓN DE LA POSTERIOR COMPLETA
#
################################################################################

y <- data_wide[3:14]
y <- as.matrix(y)

# Quitamos los encabezados para sólo tener los valores 

rownames(y) <- NULL
colnames(y) <- NULL




# *****************************************
# DEFINICIÓN DE LA POSTERIOR
# *****************************************

# Como se mencionó al inicio, este MH-within-Gibbs se hará por bloques. 
# Por lo tanto, lo óptimo es aprovechar la estructura matricial y definir los exponentes 
# de la densidad de manera matricial (bloques).

{
  
  
  # ******************************
  # Función de Poisson (verosimilitud)
  # ******************************
  poisson_prior <- function(theta, y, C) {
    # Bloques espaciales, temporales y espacio temporales
    
    alpha_t <- theta[1:12]       # Efecto temporal estructurado
    gamma_t <- theta[13:24]      # Efecto temporal no estructurado
    u_i <- theta[25:40]          # Efecto espacial estructurado
    v_i <- theta[41:56]          # Efecto espacial no estructurado 
    delta_it <- theta[57:(57 + 16 * 12 - 1)]  # Efectos espacio-temporales de Tipo I
    
    #Definición de las matrices (bloques de efectos)

    # Se van tomando entradas de 16 en 16, las cuales se ponen por columnas.
    # De modo que esta matriz sigue la estructura delta_(i,t), donde i es la región y t el tiempo.
    # Dicho de otro modo, cada columna de esta matriz hace referencia a un tiempo, y cada renglón a una alcaldía.
    delta_matrix <- matrix(delta_it, nrow = 16, ncol = 12, byrow = FALSE)
    

    alpha_matrix <- matrix(rep(alpha_t, each = 16), nrow = 16, ncol = 12, byrow = FALSE)
    gamma_matrix <- matrix(rep(gamma_t, each = 16), nrow = 16, ncol = 12, byrow = FALSE)
    u_matrix <- matrix(rep(u_i, times = 12), nrow = 16, ncol = 12, byrow = FALSE)
    v_matrix <- matrix(rep(v_i, times = 12), nrow = 16, ncol = 12, byrow = FALSE)
    
     # Se define el exponente de la "e" en la verosimilitud
     eta <- alpha_matrix + gamma_matrix + u_matrix + v_matrix + delta_matrix
     
     # se define el logaritmo de la verosimilitud
     log_poisson <- sum(y * eta - C * exp(eta))
     return(log_poisson)
    
    
  }
  
  
  # ******************************
  # Priori de alpha
  # ******************************
  alpha_prior <- function(theta) {
    
    resta_alpha <- diff(theta[1:12])^2  # Resta de consecutivos al cuadrado 
    
    log_alpha <- (-12 / 2) * log(theta[249]) - ( 1 / ( 2 *theta[249]  ) ) * sum(resta_alpha) #exponente del bloque alpha
    
    return(log_alpha)
  }
  
  # ******************************
  # Priori de gamma
  # ******************************
  gamma_prior <- function(theta) {
    gamma_bloque <- sum(theta[13:24]^2)  #Suma de cuadrados
    

    
    log_gamma <- (-12 / 2) * log(theta[250]) - (1 / (2*theta[250]) ) * gamma_bloque
    return(log_gamma)
  
    
    }
  
  # ******************************
  # Priori del bloque u
  # ******************************
  u_prior <- function(theta, W) {
    u_values <- theta[25:40]  # Se toman los valores del vector inicial
    
    # suma de los cuadrados de la resta de las u vecinas
    u_suma <- sum((outer(u_values, u_values, "-")^2) * W) / 2  # Se divide entre dos porque en W está repetida la relación de vecindad
    
    u_aux <- (-16 / 2) * log(theta[251]) - (1 / (2 *theta[251]) ) * u_suma
    return(u_aux)
    
    
  }
  
  # ******************************
  # Priori del bloque v
  # ******************************
  v_prior <- function(theta) {
    v_suma <- sum(theta[41:56]^2)  # suma de cuadrados
    
    v_aux <- (-16 / 2) * log(theta[252]) - (1 / (2 *theta[252]) ) * v_suma
    return(v_aux)
  }
  
  # ******************************
  # Priori del bloque de las deltas
  # ******************************
  delta_prior <- function(theta) {
    delta_values <- theta[57:248]  # Se toman los valores del vector inicial
    
    delta_suma <- sum(delta_values^2)  # Suma de cuadrados
    
    aux_delta <- (-16 * 12 / 2) *log(theta[253]) - (1 / (2 *theta[253]) ) * delta_suma
    return(aux_delta)
  }
  

  
  # ******************************
  # Posterior "truncada"
  # ******************************
  # Acá se define el logarimto de la posterior objetivo, pero sin lo que corresponde a los hiperparámetros.
  
  pi_posterior <- function(theta, y, C, W) {
    poisson_prior_bloque <- poisson_prior(theta, y, C)
    alpha_prior_bloque <- alpha_prior(theta)
    gamma_prior_bloque <- gamma_prior(theta)
    u_prior_bloque <- u_prior(theta, W)
    v_prior_bloque <- v_prior(theta)
    delta_prior_bloque <- delta_prior(theta)
    
    posterior_value <- poisson_prior_bloque + alpha_prior_bloque + gamma_prior_bloque + 
      u_prior_bloque + v_prior_bloque + delta_prior_bloque 
    
    return(posterior_value)
  }
  
}



################################################################################
#                 4.  Función que hace MH-within-Gibbs
################################################################################


run_MwG <- function(n_iter, y, C, W, initial_values, proposal_sds) {
  
  # Se inicializa el vector theta y unas variables auxiliares para llevar un registro de cómo va correindo el código
  theta <- initial_values$theta  
  theta <- as.numeric(theta)
  trace_theta <- matrix(0, nrow = n_iter, ncol = length(theta))  
  accepted_counts <- rep(0, length(update_blocks))  # Servirá para imprimir el radio de aceptación de las propuestas MH
  
  for (iter in 1:n_iter) {
    
    # Para cada bloque, menos el de hiperparámetros, se tomará una desviación estándar para la propuesta
    for (b in seq_along(update_blocks)) {
      
      # Para todos los bloques, menos para el de hiperparámetros se hace MwG
      if(b %in% 1:5){ 
                      
        
        
        indices <- update_blocks[[b]]  # Se establece en qué blque estamos
        sd_block <- proposal_sds[b]  # Setoma la desviación estándar para ese bloque
        
        # Se proponen valores para ese bloque
        proposed_theta <- theta
        proposed_theta[indices] <- theta[indices] + rnorm(length(indices), mean = 0, sd = sd_block)
        
        # Se calculan los elementos a comparar en MH
        current_posterior <- pi_posterior(theta, y, C, W)
        new_posterior <- pi_posterior(proposed_theta, y, C, W)
        

        # Cociente de mh
        acceptance_ratio <- exp(new_posterior - current_posterior)
        
        # En caso de aceptar la propuesta
        if (runif(1) < acceptance_ratio) {                # Sólo se acepta si la uniforme es menor. Si no, no se cambia y se queda igual.
          theta[indices] <- proposed_theta[indices]       # Así, theta trae lo más actual que se tiene como base.
          accepted_counts[b] <- accepted_counts[b] + 1    # Contador de cuántas propuestas se han aceptado para cada bloque .
        }
        
        
      }  # Cierre de los primeros 5 bloques
      
      # Ahora, acá tenemos lo que se hace para el sexto bloque (el de hiperparámetros)

      if(b==6){ 
        
        # Hasta acá, el estado actual está en theta. Sólo hay que atualizar theta[249:253]
        
        for (i in 249:253) {
          
          # lambda_alpha
          if(i==249){
            n <- 12
            # Cálculo de alpha  [parámetro shape]
            alpha_shape<- (n/2) - 1
            # Cálculo de theta [parámetro scale]
            sum_alpha_t <- diff(theta[1:12])^2
            sum_alpha_t <- sum(sum_alpha_t)
            
            eps<- 0.01
            theta_scale <- 2/(sum_alpha_t + eps) 
            
             # Simular una IG con una Gamma de shape, scale
            q<-rgamma(1,shape = alpha_shape, scale = theta_scale)
            
            # Invertir el valor para tener una IG
            q<-1/q             

            # Actualizamos theta
            theta[i]<-q

          }
          
          
          # lambda_gamma
          if(i==250){
            n <- 12
            # Cálculo de alpha  [parámetro shape]
            alpha_shape<- (n/2) - 1
            # Cálculo de theta [parámetro scale]
            
            sum_gamma_t <- sum(theta[13:24]^2)
            
            eps<- 0.01
            theta_scale <- 2/(sum_gamma_t + eps) 
            
            # Simular una IG con una Gamma de shape, scale
            q<-rgamma(1,shape = alpha_shape, scale = theta_scale)
            
            # Invertir el valor para tener una IG
            q<-1/q             
            
            # Actualizamos theta
            
            theta[i]<-q
            
          }
          
          
          # lambda_u
           if(i==251){
            
            n <- 16
            
            # Cálculo de alpha  [parámetro shape]
            alpha_shape<- (n/2) - 1
            
            # Cálculo de theta [parámetro scale]
            scale_u_aux <- theta[25:40]  # Extract spatial structured effects
            
            # Suma de los cuadrados de la resta entre us vecinas
            scale_u <- sum((outer(scale_u_aux, scale_u_aux, "-")^2) * W) / 2  # Se divide entre dos, porque la relación de vecindad se repite en w
            
            # Definición de épsilon 
            eps<- 0.01
            
            # Cálculo de beta
            theta_scale<- 2/( eps + scale_u )
            
            
            # Simular una IG con una Gamma de shape, scale
            q<-rgamma(1,shape = alpha_shape, scale = theta_scale)
            
            # Invertir el valor para tener una IG
            q<-1/q             
            
            # Actualizamos theta
            
            theta[i]<-q

           }
          
          
          # lambda_v
          if(i==252){
            
            n <- 16
            
            # Cálculo de alpha  [parámetro shape]
            alpha_shape<- (n/2) - 1
            
            # Cálculo de theta [parámetro scale]
            
            v_aux_hiper <- sum(theta[41:56]^2)  # Squared sum
            
            # Definición de épsilon 
            eps<- 0.01
            
            # Cálculo de beta
            theta_scale<- 2/(eps + v_aux_hiper)
            
            
            # Simular una IG con una Gamma de shape, scale
            q<-rgamma(1,shape = alpha_shape, scale = theta_scale)
            
            # Invertir el valor para tener una IG
            q<-1/q             
            
            # Actualizamos theta
            
            theta[i]<-q
            
          }
          
          
          # lambda_delta
          if(i==253){
            
            n <- 16*12
            
            # Cálculo de alpha  [parámetro shape]
            alpha_shape<- (n/2) - 1
            
            # Cálculo de theta [parámetro scale]
            
            delta_values_aux <- theta[57:248]  # Se toman los deltas necesarios para el cálculo
            
            delta_suma_hiper <- sum(delta_values_aux^2)  # Suma de cuadrados
            
            # Definición de épsilon 
            eps<- 0.01
            
            # Cálculo de beta
            theta_scale<- 2/(eps + delta_suma_hiper)
            
            
            # Simular una IG con una Gamma de shape, scale
            q<-rgamma(1,shape = alpha_shape, scale = theta_scale)
            
            # Invertir el valor para tener una IG
            q<-1/q             
            
            # Actualizamos theta
            
            theta[i]<-q
            
          }
          
          
        } # Se termina el for que actualiza el último bloque 
        
        
        
        
      }  # Cierre del sexto (último) bloque
      
      

    }
    
    # Se actualiza la matriz de las simulaciones
    trace_theta[iter, ] <- theta                       # Este theta ya va actualizado para el b-ésimo bloque, luego se pasa al siguiente. Es decir, en cada renglón se actualiza 6 veces.
    
    # Para monitorear cómo van las simulaciones
    if (iter %% 1000 == 0) { # Se imprime info cada 1000 simulaciones completas
      acceptance_rates <- accepted_counts / iter
      print(paste("Bloque", z,"Iteration", iter, "of", n_iter, "- Acceptance rates:", paste(round(acceptance_rates, 3), collapse = " ")))
      

      
    }
  }
  
  return(list(theta = trace_theta, acceptance_rates = accepted_counts / n_iter))
}



################################################################################
#               5.  Iterar para tener 8 M simulaciones
################################################################################

# A continuación, se define un "for" que corre 16 veces la función anterior para 
# así llegar a 8 M simulaciones. Esto se hace simulando 16 bloques de tamaño 500K

# Si se desean más o menos simulaciones, basta con modificar esos parámetos en el "for".
# Además, es posible modificar el tamaño de los bloques.



# Los nombres de las columnas es GLOBAL
# Acá se defienne los nombres de las variables en el vector theta

columnas_resultados <- c("alpha_1", "alpha_2", "alpha_3", "alpha_4", "alpha_5", "alpha_6", "alpha_7", "alpha_8", "alpha_9", "alpha_10", "alpha_11", "alpha_12",
                         "gamma_1", "gamma_2", "gamma_3", "gamma_4", "gamma_5", "gamma_6", "gamma_7", "gamma_8", "gamma_9", "gamma_10", "gamma_11", "gamma_12",
                         "u1","u2","u3","u4","u5","u6","u7","u8","u9","u10","u11","u12","u13","u14","u15","u16",
                         "v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15","v16",
                         "delta_1_1","delta_2_1","delta_3_1","delta_4_1","delta_5_1","delta_6_1","delta_7_1","delta_8_1","delta_9_1","delta_10_1","delta_11_1","delta_12_1","delta_13_1","delta_14_1","delta_15_1","delta_16_1",
                         "delta_1_2","delta_2_2","delta_3_2","delta_4_2","delta_5_2","delta_6_2","delta_7_2","delta_8_2","delta_9_2","delta_10_2","delta_11_2","delta_12_2","delta_13_2","delta_14_2","delta_15_2","delta_16_2",
                         "delta_1_3","delta_2_3","delta_3_3","delta_4_3","delta_5_3","delta_6_3","delta_7_3","delta_8_3","delta_9_3","delta_10_3","delta_11_3","delta_12_3","delta_13_3","delta_14_3","delta_15_3","delta_16_3",
                         "delta_1_4","delta_2_4","delta_3_4","delta_4_4","delta_5_4","delta_6_4","delta_7_4","delta_8_4","delta_9_4","delta_10_4","delta_11_4","delta_12_4","delta_13_4","delta_14_4","delta_15_4","delta_16_4",
                         "delta_1_5","delta_2_5","delta_3_5","delta_4_5","delta_5_5","delta_6_5","delta_7_5","delta_8_5","delta_9_5","delta_10_5","delta_11_5","delta_12_5","delta_13_5","delta_14_5","delta_15_5","delta_16_5",
                         "delta_1_6","delta_2_6","delta_3_6","delta_4_6","delta_5_6","delta_6_6","delta_7_6","delta_8_6","delta_9_6","delta_10_6","delta_11_6","delta_12_6","delta_13_6","delta_14_6","delta_15_6","delta_16_6",
                         "delta_1_7","delta_2_7","delta_3_7","delta_4_7","delta_5_7","delta_6_7","delta_7_7","delta_8_7","delta_9_7","delta_10_7","delta_11_7","delta_12_7","delta_13_7","delta_14_7","delta_15_7","delta_16_7",
                         "delta_1_8","delta_2_8","delta_3_8","delta_4_8","delta_5_8","delta_6_8","delta_7_8","delta_8_8","delta_9_8","delta_10_8","delta_11_8","delta_12_8","delta_13_8","delta_14_8","delta_15_8","delta_16_8",
                         "delta_1_9","delta_2_9","delta_3_9","delta_4_9","delta_5_9","delta_6_9","delta_7_9","delta_8_9","delta_9_9","delta_10_9","delta_11_9","delta_12_9","delta_13_9","delta_14_9","delta_15_9","delta_16_9",
                         "delta_1_10","delta_2_10","delta_3_10","delta_4_10","delta_5_10","delta_6_10","delta_7_10","delta_8_10","delta_9_10","delta_10_10","delta_11_10","delta_12_10","delta_13_10","delta_14_10","delta_15_10","delta_16_10",
                         "delta_1_11","delta_2_11","delta_3_11","delta_4_11","delta_5_11","delta_6_11","delta_7_11","delta_8_11","delta_9_11","delta_10_11","delta_11_11","delta_12_11","delta_13_11","delta_14_11","delta_15_11","delta_16_11",
                         "delta_1_12","delta_2_12","delta_3_12","delta_4_12","delta_5_12","delta_6_12","delta_7_12","delta_8_12","delta_9_12","delta_10_12","delta_11_12","delta_12_12","delta_13_12","delta_14_12","delta_15_12","delta_16_12",
                         "lambda_alpha","lambda_gamma","lambda_u","lambda_v","lambda_delta")








# Correr por bloques 

tiempo_ejecucion <- system.time({
  
  
  { 
    
    # for para hacer 16 bloques de tamaño 500K
    for(z in 1:16){
      
      
      print(paste("Vamos en el bloque ", z, " de 16"))
      
      # Primera bloque de 500k
      if (z == 1) { 
        
        #Se definen los valores iniciales (los que van cambiando entre bloques y que no se definen en la sección "2.   Definición de parámetros fijos")
        
        initial_values<-list(theta = rep(0.01,12+12+16+16+12*16+5))
        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_1 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_1 <- mcmc_results_1$theta # se toman las simulaciones
        
        resultados_1 <- as.data.frame(resultados_1) #se les da formato de data frame
        
        colnames(resultados_1) <- columnas_resultados # se le asigna el nombre a las columnas
        
        
        # Guardar como csv
        write.csv(resultados_1, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_1.csv", row.names = FALSE)
        

        
      } 
      
      # Segundo bloque de 500k
      if (z == 2) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values <- list(theta = as.numeric(resultados_1[500000, ]))
        rm(resultados_1) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_2 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_2 <- mcmc_results_2$theta
        
        resultados_2 <- as.data.frame(resultados_2)
        
        colnames(resultados_2) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_2, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_2.csv", row.names = FALSE)
        
        
      } 
      
      # Tercer bloque de 500k
      if (z == 3) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values <- list(theta = as.numeric(resultados_2[500000, ]))
        rm(resultados_2) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_3 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_3 <- mcmc_results_3$theta
        
        resultados_3 <- as.data.frame(resultados_3)
        
        colnames(resultados_3) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_3, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_3.csv", row.names = FALSE)
        
        
      } 
      
      # Cuarto bloque de 500k
      if (z == 4) { 
        
        #Se dan valores para iniciar el siguiente bloque 

        initial_values <- list(theta = as.numeric(resultados_3[500000, ]))
        rm(resultados_3) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_4 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_4 <- mcmc_results_4$theta
        
        resultados_4 <- as.data.frame(resultados_4)
        
        
        colnames(resultados_4) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_4, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_4.csv", row.names = FALSE)
        
        
      } 
      
      # Quinto bloque de 500k
      if (z == 5) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_4[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_4) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_5 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_5 <- mcmc_results_5$theta
        
        resultados_5 <- as.data.frame(resultados_5)
        
        
        colnames(resultados_5) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_5, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_5.csv", row.names = FALSE)
        
        
      } 
      
      # Sexto bloque de 500k
      if (z == 6) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_5[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_5) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_6 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_6 <- mcmc_results_6$theta
        
        resultados_6 <- as.data.frame(resultados_6)
        
        
        colnames(resultados_6) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_6, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_6.csv", row.names = FALSE)
        
        
      } 
      
      # Séptimo bloque de 500k
      if (z == 7) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_6[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_6) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_7 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_7 <- mcmc_results_7$theta
        
        resultados_7 <- as.data.frame(resultados_7)
        
        
        colnames(resultados_7) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_7, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_7.csv", row.names = FALSE)
        
        
      } 
      
      # Octavo bloque de 500k
      if (z == 8) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_7[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_7) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_8 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_8 <- mcmc_results_8$theta
        
        resultados_8 <- as.data.frame(resultados_8)
        
        
        colnames(resultados_8) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_8, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_8.csv", row.names = FALSE)
        
        
      } 
      
      # Noveno bloque de 500k
      if (z == 9) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_8[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_8) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_9 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_9 <- mcmc_results_9$theta
        
        resultados_9 <- as.data.frame(resultados_9)
        
        
        colnames(resultados_9) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_9, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_9.csv", row.names = FALSE)
        
        
      } 
      
      # Décimo bloque de 500k
      if (z == 10) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_9[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_9) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_10 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_10 <- mcmc_results_10$theta
        
        resultados_10 <- as.data.frame(resultados_10)
        
        
        colnames(resultados_10) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_10, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_10.csv", row.names = FALSE)
        
        
      } 
      
      # Undécimo bloque de 500k
      if (z == 11) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_10[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_10) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_11 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_11 <- mcmc_results_11$theta
        
        resultados_11 <- as.data.frame(resultados_11)
        
        
        colnames(resultados_11) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_11, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_11.csv", row.names = FALSE)
        
        
      } 
      
      # Doceavo bloque de 500k
      if (z == 12) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_11[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_11) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_12 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_12 <- mcmc_results_12$theta
        
        resultados_12 <- as.data.frame(resultados_12)
        
        
        colnames(resultados_12) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_12, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_12.csv", row.names = FALSE)
        
        
      } 
      
      # Treceavo bloque de 500k
      if (z == 13) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_12[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_12) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_13 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_13 <- mcmc_results_13$theta
        
        resultados_13 <- as.data.frame(resultados_13)
        
        
        colnames(resultados_13) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_13, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_13.csv", row.names = FALSE)
        
        
      } 
      
      # Catorceavo bloque de 500k
      if (z == 14) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_13[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_13) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_14 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_14 <- mcmc_results_14$theta
        
        resultados_14 <- as.data.frame(resultados_14)
        
        
        colnames(resultados_14) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_14, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_14.csv", row.names = FALSE)
        
        
      } 
      
      # Quinceavo bloque de 500k
      if (z == 15) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_14[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_14) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_15 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_15 <- mcmc_results_15$theta
        
        resultados_15 <- as.data.frame(resultados_15)
        
        
        colnames(resultados_15) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_15, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_15.csv", row.names = FALSE)
        
        
      } 
      
      
      # 16 bloque de 500k
      if (z == 16) { 
        
        #Se dan valores para iniciar el siguiente bloque 
        
        initial_values<-list(theta =resultados_15[500000, ]   ) #Tomamos el último renglón de los resultados anteriores
        rm(resultados_15) # Ahora que ya tomamos lo que ocupamos del anterior lo podemos borrar de memoria

        update_blocks <- list(
          1:12,        # alpha_t
          13:24,       # gamma_t
          25:40,       # u_i
          41:56,       # v_i
          57:(57 + 16 * 12 - 1),        # delta_it
          (57 + 16 * 12):(61 + 16 * 12)           # hypers
        )
        
        proposal_sds <- c(0.003, 0.003, 0.003, 0.003, 0.003)
        n_iter <- 500000
        
        set.seed(182435)
        mcmc_results_16 <- run_MwG(n_iter = n_iter, y, C, W, initial_values, proposal_sds)
        
        
        resultados_16 <- mcmc_results_16$theta
        
        resultados_16 <- as.data.frame(resultados_16)
        
        
        colnames(resultados_16) <- columnas_resultados
        
        
        # Guardar como csv
        write.csv(resultados_16, "/Users/emmagct/Desktop/Spatial_Temporal_SpatioTemporal/Resultados/archivo_16.csv", row.names = FALSE)
        
        
      } 
      

      
      
      
      
    
    }
    
    
    
  }
  
  
})


