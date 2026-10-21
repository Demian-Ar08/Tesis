##############################################
# Librerias y preparacion de datos 
############################################
setwd("C:/Users/Bernardo Aroca/OneDrive/Escritorio/tesis")
library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(extrafrail)
library(readxl)

df_principal <- read_excel("pacigeral (2).xlsx") 


datos_esofago <- df_principal %>%
  filter(TOPOGRUP == "C15") %>% 
  filter(ERRO == 0) %>% 
  mutate(
    fecha_diag = as.Date(DTDIAG, format="%d/%m/%Y"),
    fecha_ult  = as.Date(DTULTINFO, format="%d/%m/%Y"),
    tiempo_meses = as.numeric(difftime(fecha_ult, fecha_diag, units = "days")) / 30.44,
    evento = ifelse(ULTINFO >= 3, 1, 0),
    
    Cirugia_Total = ifelse(CIRURGIA == 1 | CIRURANT == 1 | CIRURAPOS == 1, 1, 0),
    Radio_Total   = ifelse(RADIO == 1 | RADIOANT == 1 | RADIOAPOS == 1, 1, 0),
    Quimio_Total  = ifelse(QUIMIO == 1 | QUIMIOANT == 1 | QUIMIOAPOS == 1, 1, 0),
    CLINICA_ID = as.numeric(as.factor(CLINICA))
  ) %>%
  mutate(across(c(Cirugia_Total, Radio_Total, Quimio_Total), ~ifelse(is.na(.x), 0, .x))) %>%
  filter(!is.na(tiempo_meses) & tiempo_meses > 0 & !is.na(evento))

##########################
## Kaplan- meires
######################
datos_km <- df_principal %>%
  filter(TOPOGRUP == "C15") %>% # Filtrar por Cancer de esofago
  mutate(
    f_diag = as.Date(DTDIAG, format="%d/%m/%Y"),
    f_ult  = as.Date(DTULTINFO, format="%d/%m/%Y"),
    
    # Calcular tiempo en años 
    tiempo_anios = as.numeric(difftime(f_ult, f_diag, units = "days")) / 365.25,
    
    # Evento: 1 si ULTINFO es 3 (Óbito por cáncer), 0 en otros casos
    evento = ifelse(ULTINFO == 3, 1, 0)
  ) %>%
  filter(!is.na(tiempo_anios) & tiempo_anios >= 0)

# Ajustar el modelo 
library(survival)
fit_sencillo <- survfit(Surv(tiempo_anios, evento) ~ 1, data = datos_km)

# Gráfico de Kaplan-Meier
plot(fit_sencillo,
     main = "Curva de Sobrevivencia de Kaplan-Meier\n Cáncer de esófago",
     xlab = "Tiempo (años)",
     ylab = "Probabilidad de supervivencia",
     col = "blue", 
     lwd = 2, 
     conf.int = TRUE,
     conf.int.col = "lightblue",
     mark.time = F,  
     pch = 3,           
     ylim = c(0, 1))
#Agregar una cuadrícula 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")


###################################
# Parte definitiva (Modelo)
#################################
# Preparacion de datos
library(readxl)
library(dplyr)
library(tidyr)
library(survival)
library(parfm)
library(peakRAM)
library(extrafrail) 

setwd("C:/Users/Bernardo Aroca/OneDrive/Escritorio/tesis")
#df_principal <- read_excel("pacigeral (2).xlsx") 


datos_esofago <- df_principal %>%
  filter(TOPOGRUP == "C15") %>% 
  filter(ERRO == 0) %>% 
  mutate(
    fecha_diag = as.Date(DTDIAG, format="%d/%m/%Y"),
    fecha_ult  = as.Date(DTULTINFO, format="%d/%m/%Y"),
    # Tiempo observado Y = min(T, C)
    tiempo_meses = as.numeric(difftime(fecha_ult, fecha_diag, units = "days")) / 30.44,
    # Variable indicadora de Censura/Falla (1 = Evento, 0 = Censura)
    evento = ifelse(ULTINFO >= 3, 1, 0), 
    
    # Covariables dicotómicas
    Cirugia_Total = ifelse(CIRURGIA == 1 | CIRURANT == 1 | CIRURAPOS == 1, 1, 0),
    Radio_Total   = ifelse(RADIO == 1 | RADIOANT == 1 | RADIOAPOS == 1, 1, 0),
    Quimio_Total  = ifelse(QUIMIO == 1 | QUIMIOANT == 1 | QUIMIOAPOS == 1, 1, 0),
    # El clúster será el hospital/clínica
    CLINICA_ID = as.numeric(as.factor(CLINICA))
  ) %>%
  mutate(across(c(Cirugia_Total, Radio_Total, Quimio_Total), ~ifelse(is.na(.x), 0, .x))) %>%
  filter(!is.na(tiempo_meses) & tiempo_meses > 0 & !is.na(evento))


# ###################################################
# Modelos y rendimiento computacional 
# ###################################################
# Se declaran las fórmulas especificando el clúster.
formula_extrafrail <- survival::Surv(tiempo_meses, evento) ~ Cirugia_Total + Radio_Total + Quimio_Total + Cirugia_Total:Radio_Total + Cirugia_Total:Quimio_Total + survival::cluster(CLINICA_ID)

formula_parfm <- survival::Surv(tiempo_meses, evento) ~ Cirugia_Total + Radio_Total + Quimio_Total + Cirugia_Total:Radio_Total + Cirugia_Total:Quimio_Total

# Se testean las 4 distribuciones de extrafrail (Gamma, IG, WL, BS) y 2 de parfm (Gamma, IG)
rendimiento <- peakRAM(
  Cox_Usual = { 
    mod_cox <<- coxph(formula_extrafrail, data = datos_esofago) 
  },
  Gamma_Ext = { 
    mod_gamma_ext <<- tryCatch(frailty.fit(formula_extrafrail, dist.frail = "gamma", dist = "weibull", data = datos_esofago), error = function(e) NULL) 
  },
  IG_Ext = { 
    mod_ig_ext <<- tryCatch(frailty.fit(formula_extrafrail, dist.frail = "invgauss", dist = "weibull", data = datos_esofago), error = function(e) NULL) 
  },
  WL_Ext = { 
    mod_wl_ext <<- tryCatch(frailty.fit(formula_extrafrail, dist.frail = "WL", dist = "weibull", data = datos_esofago), error = function(e) NULL) 
  },
  BS_Ext = { 
    mod_bs_ext <<- tryCatch(frailty.fit(formula_extrafrail, dist.frail = "BS", dist = "weibull", data = datos_esofago), error = function(e) NULL) 
  },
  Gamma_Parfm = { 
    mod_gamma_parfm <<- tryCatch(parfm(formula_parfm, cluster = "CLINICA_ID", data = datos_esofago, dist = "weibull", frailty = "gamma"), error = function(e) NULL) 
  },
  IG_Parfm = { 
    mod_ig_parfm <<- tryCatch(parfm(formula_parfm, cluster = "CLINICA_ID", data = datos_esofago, dist = "weibull", frailty = "ingau"), error = function(e) NULL) 
  }
)
print(rendimiento)

######################################################
# Valores de log-verosimilitud y p-valor  
#####################################################

extraer_loglik <- function(modelo, nombre) {
  if (is.null(modelo)) {
    cat(nombre, "-> Fallo de convergencia\n")
  } else {
    ll <- tryCatch(
      as.numeric(logLik(modelo)), 
      error = function(e) tryCatch(modelo$loglik[2], error = function(e) modelo$loglik)
    )
    cat(nombre, "->", ll, "\n")
  }
}

extraer_loglik(mod_cox, "Cox Usual")
extraer_loglik(mod_gamma_ext, "Gamma (Extrafrail)")
extraer_loglik(mod_ig_ext, "IG (Extrafrail)")
extraer_loglik(mod_wl_ext, "WL (Extrafrail)")
extraer_loglik(mod_bs_ext, "BS (Extrafrail)")
extraer_loglik(mod_gamma_parfm, "Gamma (Parfm)")
extraer_loglik(mod_ig_parfm, "IG (Parfm)")


imprimir_pvalores <- function(modelo, nombre_modelo) {
  if(!is.null(modelo)) {
    cat("\n---", nombre_modelo, "---\n")
    coeficientes <- modelo$coefficients
    
    errores_est <- NULL
    if(!is.null(modelo$se)) errores_est <- modelo$se
    else if(!is.null(modelo$var)) errores_est <- sqrt(diag(modelo$var))
    
    if(!is.null(errores_est)) {
      z_valores <- coeficientes / errores_est
      p_valores <- 2 * (1 - pnorm(abs(z_valores)))
      
      tabla_p <- data.frame(
        Coeficiente = round(coeficientes, 4),
        Error_Est = round(errores_est, 4),
        Valor_Z = round(z_valores, 4),
        P_Valor = round(p_valores, 5)
      )
      print(tabla_p)
    } else {
      cat("No se pudieron extraer los errores estándar.\n")
    }
  }
}

imprimir_pvalores(mod_gamma_ext, "Modelo Gamma (Extrafrail)")
imprimir_pvalores(mod_wl_ext, "Modelo WL (Extrafrail)")




####################################
# Cálculo de AIC y BIC 
######################################
calcular_criterios <- function(modelo, nombre, k_params = 8, n_obs = nrow(datos_esofago)) {
  # Extraer log-verosimilitud
  ll <- tryCatch({
    # Intento 1: logLik estándar de R
    val <- as.numeric(logLik(modelo))
    if(is.na(val) || length(val) == 0) stop()
    val
  }, error = function(e) {
    # Intento 2: Búsqueda directa en el objeto (usando logLik con L mayúscula)
    if (!is.null(modelo$logLik)) return(max(modelo$logLik))
    return(NA)
  })
  
  if (is.na(ll)) {
    return(data.frame(Modelo = nombre, LogLik = NA, AIC = NA, BIC = NA))
  }
  
  # Cálculo de AIC y BIC
  aic_val <- -2 * ll + 2 * k_params
  bic_val <- -2 * ll + k_params * log(n_obs)
  
  return(data.frame(Modelo = nombre, LogLik = ll, AIC = aic_val, BIC = bic_val))
}

resultados_aic_bic <- rbind(
  calcular_criterios(mod_gamma_parfm, "Gamma (parfm)"),
  calcular_criterios(mod_gamma_ext, "Gamma (extrafrail)"),
  calcular_criterios(mod_wl_ext, "WL (extrafrail)")
)
print(resultados_aic_bic)

###################################################
# Analisis Univariado (Log-Rank y Wilcoxon-Breslow)
###################################################
# Función para imprimir resultados más limpios
imprimir_pruebas <- function(variable, nombre_var) {
  cat(sprintf("\n--- %s ---\n", nombre_var))
  
  # Log-Rank (rho = 0)
  lr <- survdiff(Surv(tiempo_meses, evento) ~ eval(parse(text=variable)), data = datos_esofago, rho = 0)
  p_lr <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  
  # Wilcoxon-Breslow (rho = 1)
  wil <- survdiff(Surv(tiempo_meses, evento) ~ eval(parse(text=variable)), data = datos_esofago, rho = 1)
  p_wil <- 1 - pchisq(wil$chisq, length(wil$n) - 1)
  
  cat(sprintf("Log-Rank: Chi-cuadrado = %.2f, p-valor = %.5f\n", lr$chisq, p_lr))
  cat(sprintf("Wilcoxon: Chi-cuadrado = %.2f, p-valor = %.5f\n", wil$chisq, p_wil))
}

# Ejecutar para las 3 variables
imprimir_pruebas("Cirugia_Total", "CIRUGÍA TOTAL")
imprimir_pruebas("Radio_Total", "RADIOTERAPIA TOTAL")
imprimir_pruebas("Quimio_Total", "QUIMIOTERAPIA TOTAL")

######################################################################
# FRAGILIDADES (Zi) Y ANÁLISIS DESCRIPTIVO (VERSIÓN MEJORADA)
#################################################################
library(dplyr)

# Extraer nombres de clínicas y la fragilidad del modelo Gamma (extrafrail)
nombres_clinicas <- levels(as.factor(datos_esofago$CLINICA))
fragilidades <- mod_gamma_ext$frail 
if(is.null(fragilidades)) fragilidades <- mod_gamma_ext$z 

# Crear dataframe base con las fragilidades
df_fragilidades <- data.frame(
  CLINICA = as.character(nombres_clinicas),
  Zi = round(fragilidades, 4),
  Riesgo_Extra_pct = round((fragilidades - 1) * 100, 1) 
)

# Calcular estadísticas descriptivas generales para TODAS las clínicas cruzadas con su fragilidad (Z)
analisis_todas <- datos_esofago %>%
  mutate(CLINICA = as.character(CLINICA)) %>% # Evitar el error de cruce de datos
  group_by(CLINICA) %>%
  summarise(
    Pacientes = n(),
    Mortalidad_pct = round(mean(evento) * 100, 1),
    Mediana_Meses = round(median(tiempo_meses), 1),
    Cirugia_pct = round(mean(Cirugia_Total) * 100, 1),
    Quimio_pct = round(mean(Quimio_Total) * 100, 1),
    Radio_pct = round(mean(Radio_Total) * 100, 1),
    .groups = 'drop'
  ) %>%
  left_join(df_fragilidades, by = "CLINICA") %>%
  arrange(Zi) 

cat("\n--- EXTREMOS BRUTOS: TOP 5 MÁS SEGURAS (Con ruido de muestras pequeñas) ---\n")
print(head(analisis_todas, 5))

cat("\n--- EXTREMOS BRUTOS: TOP 5 MÁS RIESGOSAS (Con ruido de muestras pequeñas) ---\n")
print(tail(analisis_todas, 5))

################################################################################
# SELECCIÓN REPRESENTATIVA (Para caracterización y gráficos finales)
###############################################################################
#Clínicas con volumen intermedio-alto equivalente (500 a 800 pacientes)
clinicas_representativas <- analisis_todas %>% 
  filter(Pacientes >= 500 & Pacientes <= 800)

print("\n--- MAS SEGURA REPRESENTATIVA ---\n")
print(head(clinicas_representativas, 1))

print("\n--- MAS RIESGOSA REPRESENTATIVA ---\n")
print(tail(clinicas_representativas, 1))

# Crear el resumen final 
analisis_extremos_elegidos <- clinicas_representativas %>%
  filter(CLINICA %in% c(head(clinicas_representativas$CLINICA, 1), tail(clinicas_representativas$CLINICA, 1))) %>%
  mutate(Perfil = ifelse(Zi < 1, "Protectora (Z < 1)", "Riesgosa (Z > 1)")) %>%
  relocate(Perfil) %>%
  arrange(Zi)

print(as.data.frame(analisis_extremos_elegidos))

####################################
# Cluster
#####################################
if(length(nombres_clinicas) == length(fragilidades)) {
  # Preparar el DataFrame
  df_grafico <- tibble(
    Clinica = nombres_clinicas, 
    Fragilidad_Z = fragilidades,            # El valor bruto Z_i
    Riesgo_Extra = fragilidades - 1         # La desviación
  ) %>%
    arrange(Riesgo_Extra) %>% 
    mutate(
      Clinica = factor(Clinica, levels = Clinica),
      # Crear una variable para colorear fácilmente si es protector o riesgoso
      Perfil = ifelse(Riesgo_Extra > 0, "Riesgosa (Z > 1)", "Protectora (Z < 1)")
    )
  
  # Crear el Gráfico con mejoras visuales
  grafico_cluster <- ggplot(df_grafico, aes(x = Clinica, y = Riesgo_Extra, color = Perfil)) +
    geom_segment(aes(x = Clinica, xend = Clinica, y = 0, yend = Riesgo_Extra), color = "gray80") + # Líneas guía (Lollipop chart)
    geom_point(size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) + 
    coord_flip() +
    scale_color_manual(values = c("Riesgosa (Z > 1)" = "#d73027", "Protectora (Z < 1)" = "#4575b4")) +
    labs(
      title = "Estimación de Fragilidad Institucional en Cáncer de Esófago",
      subtitle = "Desviación del Riesgo Basal (Z - 1) por Centro de Atención (Modelo Gamma)",
      x = "Código del Centro de Atención (Clínica)", 
      y = "Desviación del Riesgo Relativo"
    ) +
    theme_minimal() + 
    theme(
      axis.text.y = element_text(size = 8, face = "bold"), # Tamaño legible para 40 clínicas
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  ggsave(filename = "cluster_esofago_alta_resolucion.pdf", 
         plot = grafico_cluster, 
         width = 10,       
         height = 12,      
         dpi = 300)        
  
  cat("\nEl gráfico de clúster se ha guardado en tu directorio de trabajo como 'cluster_esofago_alta_resolucion.pdf'\n")
  
} else {
  cat("\nError: La cantidad de clínicas no coincide con la cantidad de fragilidades estimadas.\n")
}


####################################################
# Graficos de sobrevivencia y riesgo (36 vs 9)
########################################################
library(ggplot2)
library(dplyr)

# Parámetros del modelo Gamma Extrafrail
rho <- 0.8496
lambda <- 0.2104
beta_cir <- -1.3488
beta_rad <- -0.5483

# Clinicas seleccionadas
z_segura <- 0.7593     # Clínica 36
z_riesgosa <- 1.0973   # Clínica 9

# Nombres descriptivos para los paneles del gráfico
nombre_segura <- "Clínica 36 (Protectora) | Z = 0.759 | n = 586"
nombre_riesgosa <- "Clínica 9 (Riesgosa) | Z = 1.097 | n = 760"

# Exponenciales de los tratamientos
exp_base <- exp(0)           # Sin Tratamiento
exp_cir <- exp(beta_cir)     # Solo Cirugía
exp_rad <- exp(beta_rad)     # Solo Radioterapia

# Secuencia de tiempo (0.1 a 60 meses)
t <- seq(0.1, 60, by = 0.5)

# Funciones matematicas de Riesgo y Sobrevivencia
calc_hazard <- function(t, z, exp_xb) { z * lambda * rho * (t^(rho - 1)) * exp_xb }
calc_surv   <- function(t, z, exp_xb) { exp(-z * lambda * (t^rho) * exp_xb) }

escenarios <- expand.grid(
  Tiempo = t,
  Cluster = c(nombre_segura, nombre_riesgosa),
  Tratamiento = c("Sin Tratamiento", "Solo Cirugía", "Solo Radioterapia")
)

df_curvas <- escenarios %>%
  mutate(
    Zi = ifelse(Cluster == nombre_segura, z_segura, z_riesgosa),
    Exp_XB = case_when(
      Tratamiento == "Sin Tratamiento" ~ exp_base,
      Tratamiento == "Solo Cirugía" ~ exp_cir,
      Tratamiento == "Solo Radioterapia" ~ exp_rad
    ),
    Riesgo = calc_hazard(Tiempo, Zi, Exp_XB),
    Sobrevivencia = calc_surv(Tiempo, Zi, Exp_XB)
  )
# orden de leyenda
df_curvas$Tratamiento <- factor(df_curvas$Tratamiento, levels = c("Sin Tratamiento", "Solo Radioterapia", "Solo Cirugía"))

##########################################
# GRÁFICO 1: SOBREVIVENCIA S(t)
#############################################
plot_surv <- ggplot(df_curvas, aes(x = Tiempo, y = Sobrevivencia, color = Tratamiento)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Cluster) + # Divide el gráfico en dos paneles
  scale_color_manual(values = c("Sin Tratamiento" = "#d73027", "Solo Cirugía" = "#4575b4", "Solo Radioterapia" = "#fdae61")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Probabilidad de Sobrevivencia S(t) por Perfil Institucional",
    subtitle = "Comparación de recintos representativos de volumen equivalente",
    x = "Tiempo de seguimiento (meses)", 
    y = "Probabilidad de Sobrevivencia S(t)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 13), 
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "gray95"))

print(plot_surv)
ggsave("Sobrevivencia_Representativos.pdf", plot = plot_surv, width = 10, height = 5, dpi = 300)

#############################################
# GRÁFICO 2: RIESGO CONDICIONAL h(t)
##############################################
plot_hazard <- ggplot(df_curvas, aes(x = Tiempo, y = Riesgo, color = Tratamiento)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Cluster) + # Divide el gráfico en dos paneles
  scale_color_manual(values = c("Sin Tratamiento" = "#d73027", "Solo Cirugía" = "#4575b4", "Solo Radioterapia" = "#fdae61")) +
  coord_cartesian(ylim = c(0, 0.4)) + 
  labs(
    title = "Función de Riesgo Condicional h(t) por Perfil Institucional",
    subtitle = expression(paste("Riesgo basal decreciente (Weibull, ", rho, " < 1)")),
    x = "Tiempo de seguimiento (meses)", 
    y = "Riesgo Condicional h(t)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 13), 
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "gray95"))

print(plot_hazard)
ggsave("Riesgo_Representativos.pdf", plot = plot_hazard, width = 10, height = 5, dpi = 300)








