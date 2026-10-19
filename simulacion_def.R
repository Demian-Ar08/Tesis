########################################################################################
setwd("C:/Users/Bernardo Aroca/OneDrive/Escritorio/tesis")

library(extrafrail)
library(survival)
library(dplyr)
library(peakRAM)
library(statmod) 
library(writexl) 
library(frailtypack)
library(parfm)

########################################
# Función de simulación 
#######################################
simular_datos_fragilidad <- function(n_pacientes_vec, dist_fragilidad, theta, lambda_w, rho_w, beta_cov1, beta_cov2, tasa_censura_objetivo) {
  n_clinicas <- length(n_pacientes_vec)
  
  if (dist_fragilidad == "gamma") {
    Z_clinicas <- rgamma(n_clinicas, shape = 1/theta, rate = 1/theta)
  } else if (dist_fragilidad == "WL") {
    Z_clinicas <- rWL(n_clinicas, theta = theta) 
  } else if (dist_fragilidad == "BS") {
    phi_theta <- (1 + sqrt(3 * theta + 1) - theta) / theta
    alpha_bs <- sqrt(2 / phi_theta)
    beta_bs <- phi_theta / (phi_theta + 1)
    Z_norm <- rnorm(n_clinicas)
    Z_clinicas <- beta_bs * (0.5 * alpha_bs * Z_norm + sqrt((0.5 * alpha_bs * Z_norm)^2 + 1))^2
  } else if (dist_fragilidad == "IG") {
    Z_clinicas <- statmod::rinvgauss(n_clinicas, mean = 1, shape = 1/theta)
  } else { stop("Distribución no soportada.") }
  
  Z_pacientes <- rep(Z_clinicas, times = n_pacientes_vec)
  id_clinica <- rep(1:n_clinicas, times = n_pacientes_vec)
  total_pacientes <- sum(n_pacientes_vec)
  
  X1 <- rnorm(total_pacientes, mean = 0, sd = 1)        
  X2 <- rbinom(total_pacientes, size = 1, prob = 0.5)   
  
  U <- runif(total_pacientes)
  tiempo_falla <- (-log(U) / (Z_pacientes * lambda_w * exp(X1 * beta_cov1 + X2 * beta_cov2)))^(1/rho_w)
  
  ################################################
  # RESOLUCIÓN NUMÉRICA DE CENSURA INDIVIDUAL 
  ####################################################
  eta_vec <- (Z_pacientes * lambda_w * exp(X1 * beta_cov1 + X2 * beta_cov2))^(-1/rho_w)
  
  lambda_c_vec <- sapply(1:total_pacientes, function(j) {
    eta_j <- eta_vec[j]
    
    f_obj_indiv <- function(lc) {
      integrand <- function(t) {
        (1 - exp(-lc * t)) * dweibull(t, shape = rho_w, scale = eta_j)
      }
      res <- tryCatch(integrate(integrand, lower = 0, upper = Inf)$value, error = function(e) NA)
      if(is.na(res)) return(100) 
      return(res - tasa_censura_objetivo)
    }
    
    tryCatch({
      uniroot(f_obj_indiv, lower = 1e-8, upper = 100)$root
    }, error = function(e) {
      0.05 
    })
  })
  ##
  tiempo_censura <- rexp(total_pacientes, rate = lambda_c_vec) 
  tiempo_observado <- pmin(tiempo_falla, tiempo_censura)
  evento <- ifelse(tiempo_falla <= tiempo_censura, 1, 0)
  
  datos_simulados <- data.frame(
    CLINICA = as.factor(id_clinica),
    Fragilidad_Z = Z_pacientes,
    Covariable_X1 = X1,
    Covariable_X2 = X2,
    Tiempo = tiempo_observado,
    Evento = evento
  )
  return(datos_simulados)
}

#######################################
# Parámetros, 8 Escenarios Paramétricos y Clústeres
#######################################
N_simulaciones <- 1000  
tasas_censura <- c(0.01, 0.05, 0.20) 

# Definición de los 8 casos exigidos
betas_opciones <- list(c(0.8, -0.5), c(-0.7, 0.9))
weibull_opciones <- list(c(0.1, 1.2), c(0.05, 1.5))
thetas_opciones <- c(0.5, 1.0)

grid_parametros <- expand.grid(
  id_beta = 1:2,
  theta = thetas_opciones,
  id_weibull = 1:2
)

escenario_1_base <- c(rep(22, 4), rep(20, 10), rep(19, 10), rep(18, 10), rep(17, 6), rep(15, 6))
escenario_2_doble_pac <- escenario_1_base * 2
escenario_3_doble_cli <- c(escenario_1_base, escenario_1_base)

lista_escenarios <- list(
  "Escenario 1 (Base)" = escenario_1_base,
  "Escenario 2 (Doble Pacientes)" = escenario_2_doble_pac,
  "Escenario 3 (Doble Clinicas)" = escenario_3_doble_cli
)

# Función auxiliar para calcular métricas
calcular_metricas <- function(est, se_est, valor_real, Nc) {
  bias <- mean(est) - valor_real
  ecm <- mean((est - valor_real)^2)
  se_empirico <- ifelse(Nc > 1, sd(est), 0) 
  
  if(any(is.na(se_est))) {
    pc <- NA
  } else {
    lim_inf <- est - 1.96 * se_est
    lim_sup <- est + 1.96 * se_est
    pc <- mean(valor_real >= lim_inf & valor_real <= lim_sup) * 100
  }
  return(c(Bias = bias, SE = se_empirico, ECM = ecm, PC = pc))
}

###################################
# 1. BUCLE EXTRAFRAIL
###################################

resultados_extrafrail <- data.frame()

for (caso in 1:nrow(grid_parametros)) {
  
  beta1_real <- betas_opciones[[grid_parametros$id_beta[caso]]][1]
  beta2_real <- betas_opciones[[grid_parametros$id_beta[caso]]][2]
  theta_real <- grid_parametros$theta[caso]
  lambda_real <- weibull_opciones[[grid_parametros$id_weibull[caso]]][1]
  rho_real <- weibull_opciones[[grid_parametros$id_weibull[caso]]][2]
  
  cat(sprintf("\n=======================================================\n"))
  cat(sprintf(" EXTRAFRAIL | CASO %d DE 8 | (B1=%.1f, B2=%.1f, Th=%.1f, L=%.2f, p=%.1f) \n", 
              caso, beta1_real, beta2_real, theta_real, lambda_real, rho_real))
  cat(sprintf("=======================================================\n"))
  
  for (nombre_escenario in names(lista_escenarios)) {
    tam_actuales <- lista_escenarios[[nombre_escenario]]
    cat(sprintf("\n*** %s | EXTRAFRAIL ***\n", toupper(nombre_escenario)))
    
    modelos_a_evaluar <- c("gamma", "IG", "WL", "BS")
    
    for (modelo in modelos_a_evaluar) {
      cat(sprintf("\n-> Modelo: %s\n", toupper(modelo)))
      
      for(tc in tasas_censura) {
        est_beta1 <- numeric(N_simulaciones); se_beta1 <- numeric(N_simulaciones)
        est_beta2 <- numeric(N_simulaciones); se_beta2 <- numeric(N_simulaciones)
        est_theta <- numeric(N_simulaciones); se_theta <- numeric(N_simulaciones)
        t_ejec <- numeric(N_simulaciones); mem_uso <- numeric(N_simulaciones)
        
        set.seed(2026) 
        
        for (i in 1:N_simulaciones) {
          datos_sim <- simular_datos_fragilidad(tam_actuales, modelo, theta_real, lambda_real, rho_real, beta1_real, beta2_real, tc)
          datos_sim$CLINICA <- as.numeric(as.character(datos_sim$CLINICA))
          
          ajuste_seguro <- NULL
          medicion <- peakRAM({
            temp <- tryCatch({ 
              frailty.fit(Surv(Tiempo, Evento) ~ Covariable_X1 + Covariable_X2 + cluster(CLINICA), 
                          dist.frail = modelo, dist = "weibull", data = datos_sim) 
            }, error = function(e) NULL)
            ajuste_seguro <<- temp 
          })
          
          if (!is.null(ajuste_seguro) && length(ajuste_seguro$coefficients) >= 5) { 
            coefs <- ajuste_seguro$coefficients
            err_est <- if(!is.null(ajuste_seguro$se)) ajuste_seguro$se else if(!is.null(ajuste_seguro$var)) sqrt(diag(ajuste_seguro$var)) else rep(NA, length(coefs))
            
            est_beta1[i] <- coefs[1]; se_beta1[i] <- err_est[1]
            est_beta2[i] <- coefs[2]; se_beta2[i] <- err_est[2]
            est_theta[i] <- coefs[5]; se_theta[i] <- err_est[5]
            t_ejec[i] <- medicion$Elapsed_Time_sec; mem_uso[i] <- medicion$Peak_RAM_Used_MiB
          } else { 
            est_beta1[i] <- NA; se_beta1[i] <- NA
            est_beta2[i] <- NA; se_beta2[i] <- NA
            est_theta[i] <- NA; se_theta[i] <- NA
            t_ejec[i] <- NA; mem_uso[i] <- NA 
          }
        }
        
        validos <- !is.na(est_theta)
        Nc <- sum(validos)
        TC <- (Nc / N_simulaciones) * 100
        
        if (Nc > 0) {
          m_b1 <- calcular_metricas(est_beta1[validos], se_beta1[validos], beta1_real, Nc)
          m_b2 <- calcular_metricas(est_beta2[validos], se_beta2[validos], beta2_real, Nc)
          m_th <- calcular_metricas(est_theta[validos], se_theta[validos], theta_real, Nc)
          
          cat(sprintf("  Cen: %.2f | Conv: %3.1f%% | T/iter: %.3fs | Mem: %.2fMB\n", tc, TC, mean(t_ejec[validos]), mean(mem_uso[validos])))
          
          nueva_fila <- data.frame(
            Caso = caso, Escenario = nombre_escenario, Paquete = "extrafrail", Modelo = modelo, Censura = tc,
            Beta1_True = beta1_real, Beta2_True = beta2_real, Theta_True = theta_real, Lambda_True = lambda_real, Rho_True = rho_real,
            Convergencia_Pct = TC, Tiempo_Segundos = mean(t_ejec[validos]), Memoria_MB = mean(mem_uso[validos]),
            Beta1_Bias = m_b1[1], Beta1_SE = m_b1[2], Beta1_ECM = m_b1[3], Beta1_PC = m_b1[4],
            Beta2_Bias = m_b2[1], Beta2_SE = m_b2[2], Beta2_ECM = m_b2[3], Beta2_PC = m_b2[4],
            Theta_Bias = m_th[1], Theta_SE = m_th[2], Theta_ECM = m_th[3], Theta_PC = m_th[4]
          )
        } else {
          cat(sprintf("  Cen: %.2f | Conv: 0.0%% | El modelo no convergió.\n", tc))
          nueva_fila <- data.frame(
            Caso = caso, Escenario = nombre_escenario, Paquete = "extrafrail", Modelo = modelo, Censura = tc,
            Beta1_True = beta1_real, Beta2_True = beta2_real, Theta_True = theta_real, Lambda_True = lambda_real, Rho_True = rho_real,
            Convergencia_Pct = 0, Tiempo_Segundos = NA, Memoria_MB = NA,
            Beta1_Bias = NA, Beta1_SE = NA, Beta1_ECM = NA, Beta1_PC = NA,
            Beta2_Bias = NA, Beta2_SE = NA, Beta2_ECM = NA, Beta2_PC = NA,
            Theta_Bias = NA, Theta_SE = NA, Theta_ECM = NA, Theta_PC = NA
          )
        }
        resultados_extrafrail <- rbind(resultados_extrafrail, nueva_fila)
      }
    }
  }
}

write_xlsx(resultados_extrafrail, "Resultados_Extrafrail_8Casos.xlsx")

##################################
# 2. BUCLE COMPETIDORES (SOLO PARFM)
##################################
resultados_competidores <- data.frame()

for (caso in 8:nrow(grid_parametros)) {
  
  beta1_real <- betas_opciones[[grid_parametros$id_beta[caso]]][1]
  beta2_real <- betas_opciones[[grid_parametros$id_beta[caso]]][2]
  theta_real <- grid_parametros$theta[caso]
  lambda_real <- weibull_opciones[[grid_parametros$id_weibull[caso]]][1]
  rho_real <- weibull_opciones[[grid_parametros$id_weibull[caso]]][2]
  
  cat(sprintf("\n=======================================================\n"))
  cat(sprintf(" COMPETIDORES | CASO %d DE 8 | (B1=%.1f, B2=%.1f, Th=%.1f, L=%.2f, p=%.1f) \n", 
              caso, beta1_real, beta2_real, theta_real, lambda_real, rho_real))
  cat(sprintf("=======================================================\n"))
  
  for (nombre_escenario in names(lista_escenarios)) {
    tam_actuales <- lista_escenarios[[nombre_escenario]]
    cat(sprintf("\n*** %s | COMPETIDORES ***\n", toupper(nombre_escenario)))
    
    # Solo parfm con Gamma e Inversa Gaussiana
    configuraciones_comp <- list(
      list(pkg = "parfm", modelo = "gamma"),
      list(pkg = "parfm", modelo = "ingau") 
    )
    
    for (config in configuraciones_comp) {
      paquete <- config$pkg
      modelo <- config$modelo
      cat(sprintf("\n-> Paquete: %s | Modelo: %s\n", toupper(paquete), toupper(modelo)))
      
      for(tc in tasas_censura) {
        est_beta1 <- numeric(N_simulaciones); se_beta1 <- numeric(N_simulaciones)
        est_beta2 <- numeric(N_simulaciones); se_beta2 <- numeric(N_simulaciones)
        est_theta <- numeric(N_simulaciones); se_theta <- numeric(N_simulaciones)
        t_ejec <- numeric(N_simulaciones); mem_uso <- numeric(N_simulaciones)
        
        set.seed(2026)
        
        for (i in 1:N_simulaciones) {
          dist_gen <- ifelse(modelo == "ingau", "IG", "gamma")
          datos_sim <- simular_datos_fragilidad(tam_actuales, dist_gen, theta_real, lambda_real, rho_real, beta1_real, beta2_real, tc)
          datos_sim$CLINICA <- as.numeric(as.character(datos_sim$CLINICA))
          
          ajuste_seguro <- NULL
          
          # Único bloque de estimación: parfm
          medicion <- peakRAM({
            temp <- tryCatch({ 
              parfm(Surv(Tiempo, Evento) ~ Covariable_X1 + Covariable_X2, 
                    cluster = "CLINICA", data = datos_sim, dist = "weibull", frailty = modelo) 
            }, error = function(e) NULL)
            ajuste_seguro <<- temp
          })
          
          if (!is.null(ajuste_seguro) && !inherits(ajuste_seguro, "try-error") && "theta" %in% rownames(ajuste_seguro)) {
            est_beta1[i] <- ajuste_seguro["Covariable_X1", "ESTIMATE"]; se_beta1[i] <- ajuste_seguro["Covariable_X1", "SE"]
            est_beta2[i] <- ajuste_seguro["Covariable_X2", "ESTIMATE"]; se_beta2[i] <- ajuste_seguro["Covariable_X2", "SE"]
            est_theta[i] <- ajuste_seguro["theta", "ESTIMATE"]; se_theta[i] <- ajuste_seguro["theta", "SE"]
            t_ejec[i] <- medicion$Elapsed_Time_sec; mem_uso[i] <- medicion$Peak_RAM_Used_MiB
          } else { 
            est_beta1[i] <- NA; se_beta1[i] <- NA; est_beta2[i] <- NA; se_beta2[i] <- NA
            est_theta[i] <- NA; se_theta[i] <- NA; t_ejec[i] <- NA; mem_uso[i] <- NA 
          }
        }
        
        validos <- !is.na(est_theta)
        Nc <- sum(validos)
        TC <- (Nc / N_simulaciones) * 100
        
        if (Nc > 0) {
          m_b1 <- calcular_metricas(est_beta1[validos], se_beta1[validos], beta1_real, Nc)
          m_b2 <- calcular_metricas(est_beta2[validos], se_beta2[validos], beta2_real, Nc)
          m_th <- calcular_metricas(est_theta[validos], se_theta[validos], theta_real, Nc)
          
          cat(sprintf("  Cen: %.2f | Conv: %3.1f%% | T/iter: %.3fs | Mem: %.2fMB\n", tc, TC, mean(t_ejec[validos]), mean(mem_uso[validos])))
          
          nueva_fila <- data.frame(
            Caso = caso, Escenario = nombre_escenario, Paquete = paquete, Modelo = modelo, Censura = tc,
            Beta1_True = beta1_real, Beta2_True = beta2_real, Theta_True = theta_real, Lambda_True = lambda_real, Rho_True = rho_real,
            Convergencia_Pct = TC, Tiempo_Segundos = mean(t_ejec[validos]), Memoria_MB = mean(mem_uso[validos]),
            Beta1_Bias = m_b1[1], Beta1_SE = m_b1[2], Beta1_ECM = m_b1[3], Beta1_PC = m_b1[4],
            Beta2_Bias = m_b2[1], Beta2_SE = m_b2[2], Beta2_ECM = m_b2[3], Beta2_PC = m_b2[4],
            Theta_Bias = m_th[1], Theta_SE = m_th[2], Theta_ECM = m_th[3], Theta_PC = m_th[4]
          )
        } else {
          cat(sprintf("  Cen: %.2f | Conv: 0.0%% | El modelo no convergió.\n", tc))
          nueva_fila <- data.frame(
            Caso = caso, Escenario = nombre_escenario, Paquete = paquete, Modelo = modelo, Censura = tc,
            Beta1_True = beta1_real, Beta2_True = beta2_real, Theta_True = theta_real, Lambda_True = lambda_real, Rho_True = rho_real,
            Convergencia_Pct = 0, Tiempo_Segundos = NA, Memoria_MB = NA,
            Beta1_Bias = NA, Beta1_SE = NA, Beta1_ECM = NA, Beta1_PC = NA,
            Beta2_Bias = NA, Beta2_SE = NA, Beta2_ECM = NA, Beta2_PC = NA,
            Theta_Bias = NA, Theta_SE = NA, Theta_ECM = NA, Theta_PC = NA
          )
        }
        resultados_competidores <- rbind(resultados_competidores, nueva_fila)
        
        # EL CHECKPOINT
        write.csv(resultados_competidores, "Resultados_Competidores_Backup_parte3.csv", row.names = FALSE)
      }
      gc()
    }
  }
}

write_xlsx(resultados_competidores, "Resultados_Competidores_8Casos.xlsx")


