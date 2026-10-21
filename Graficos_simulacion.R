#######################################
# Graficos para la simulación
#######################################
setwd("C:/Users/Bernardo Aroca/OneDrive/Escritorio/tesis")
library(readxl)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel) 

#############################################
# Limpieza de base de datos parfm
#############################################
limpiar_errores_parfm <- function(x) {
  x <- as.character(x)
  idx_e <- grepl("E|e", x)
  x[idx_e] <- gsub(",", ".", x[idx_e])
  idx_coma <- grepl(",", x) & !idx_e
  x[idx_coma] <- gsub("\\.", "", x[idx_coma])
  x[idx_coma] <- gsub(",", ".", x[idx_coma])  
  idx_puntos <- grepl("\\..*\\.", x) & !idx_e
  x[idx_puntos] <- sapply(x[idx_puntos], function(val) {
    partes <- unlist(strsplit(val, "\\."))
    paste0(partes[1], ".", paste0(partes[-1], collapse=""))
  })
  return(as.numeric(x))
}

format_num_grafico <- function(x) {
  sapply(x, function(val) {
    if(is.na(val)) return(NA)
    if(abs(val) < 0.0001 & val != 0 | abs(val) > 10000) {
      return(formatC(val, format = "e", digits = 1))
    } else {
      return(sprintf("%.4f", val))
    }
  })
}

############################################
# Union de base de datos 
##########################################
df_oficial <- read_excel("Resultados_Extrafrail_8Casos_oficial.xlsx")
df_comp <- read_delim("Resultados_Competidores_Backup_parte2.csv", delim = ";", col_types = cols(.default = "c"))

df_comp_limpio <- df_comp %>%
  mutate(
    Caso = as.numeric(Caso), # <-- SOLUCIÓN AL ERROR DE UNIÓN
    across(c(Censura, Tiempo_Segundos, Memoria_MB, contains("Bias"), contains("SE"), contains("ECM")), limpiar_errores_parfm),
    Modelo = ifelse(Modelo == "ingau", "IG", Modelo),
    Modelo = ifelse(Modelo == "gamma", "Gamma", Modelo)
  )

df_oficial_limpio <- df_oficial %>%
  mutate(
    Caso = as.numeric(Caso),
    Modelo = tools::toTitleCase(tolower(Modelo)),
    Modelo = ifelse(Modelo == "Ig", "IG", Modelo)
  )

vars_select <- c("Caso", "Escenario", "Censura", "Paquete", "Modelo", 
                 "Tiempo_Segundos", "Memoria_MB",
                 "Beta1_Bias", "Beta1_SE", "Beta1_ECM",
                 "Beta2_Bias", "Beta2_SE", "Beta2_ECM",
                 "Theta_Bias", "Theta_SE", "Theta_ECM")

df_total <- bind_rows(
  df_oficial_limpio %>% dplyr::select(all_of(vars_select)),
  df_comp_limpio %>% dplyr::select(all_of(vars_select))
)

###################################################
# Funcion para general los graficos de panel 
########################################################
generar_grafico_panel <- function(datos, y_var, titulo, subtitulo, colores) {
  pd <- position_dodge(width = 0.3)
  
  p <- ggplot(datos, aes(x = factor(Censura), y = .data[[y_var]], color = Paquete, group = Paquete)) +
    geom_line(linewidth = 0.8, position = pd, na.rm = TRUE) +
    geom_point(size = 2, position = pd, na.rm = TRUE) +
    geom_text_repel(aes(label = format_num_grafico(.data[[y_var]])), 
                    position = pd, 
                    size = 2.8, 
                    box.padding = 0.4, 
                    point.padding = 0.3,
                    segment.color = 'grey50',
                    na.rm = TRUE,
                    show.legend = FALSE) +
    facet_wrap(~ Caso, ncol = 4, scales = "free_y", labeller = labeller(Caso = function(x) paste("Caso", x))) +
    scale_color_manual(values = colores) +
    labs(title = titulo, subtitle = subtitulo, x = "Nivel de Censura", y = "Valor") +
    theme_bw() +
    theme(
      legend.position = "bottom", 
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 13)
    )
  
  return(p)
}

###############################################
# Buble para graficos de simulacion
##############################################
escenarios <- unique(df_total$Escenario)
modelos <- c("Gamma", "IG")

# Definir paletas de colores
col_b1 <- c("extrafrail" = "#E63946", "parfm" = "#1D3557") 
col_b2 <- c("extrafrail" = "#F4A261", "parfm" = "#2A9D8F") 
col_th <- c("extrafrail" = "#9C6644", "parfm" = "#5C677D") 
col_cp <- c("extrafrail" = "#8338EC", "parfm" = "#FF006E") 


for(esc in escenarios) {
  for(mod in modelos) {
    
    # Filtrar datos
    df_plot <- df_total %>% filter(Escenario == esc, Modelo == mod)
    if(nrow(df_plot) == 0) next
    
    # Nombre limpio para exportar el archivo
    nombre_limpio_esc <- gsub("[ ()]", "", esc)
    prefijo_archivo <- paste0(nombre_limpio_esc, "_", mod)
    subtitulo_graf <- paste0("Modelo ", mod, " | ", esc)
    
    #Graficos de Tiempo y Memoria RAM
    p_time <- generar_grafico_panel(df_plot, "Tiempo_Segundos", "Tiempo de Ejecución (s)", subtitulo_graf, col_cp)
    ggsave(paste0("Tiempo_", prefijo_archivo, ".pdf"), plot = p_time, width = 11, height = 7, dpi = 300)
    
    p_ram <- generar_grafico_panel(df_plot, "Memoria_MB", "Uso Máximo de Memoria RAM (MB)", subtitulo_graf, col_cp)
    ggsave(paste0("RAM_", prefijo_archivo, ".pdf"), plot = p_ram, width = 11, height = 7, dpi = 300)
    
    #Graficos Beta 1
    p_b1_b <- generar_grafico_panel(df_plot, "Beta1_Bias", "Sesgo (Bias) - Beta 1", subtitulo_graf, col_b1)
    ggsave(paste0("Beta1_Bias_", prefijo_archivo, ".pdf"), plot = p_b1_b, width = 11, height = 7, dpi = 300)
    
    p_b1_se <- generar_grafico_panel(df_plot, "Beta1_SE", "Error Estándar (SE) - Beta 1", subtitulo_graf, col_b1)
    ggsave(paste0("Beta1_SE_", prefijo_archivo, ".pdf"), plot = p_b1_se, width = 11, height = 7, dpi = 300)
    
    p_b1_ecm <- generar_grafico_panel(df_plot, "Beta1_ECM", "Error Cuadrático Medio (ECM) - Beta 1", subtitulo_graf, col_b1)
    ggsave(paste0("Beta1_ECM_", prefijo_archivo, ".pdf"), plot = p_b1_ecm, width = 11, height = 7, dpi = 300)
    
    # Graficos Beta 2
    p_b2_b <- generar_grafico_panel(df_plot, "Beta2_Bias", "Sesgo (Bias) - Beta 2", subtitulo_graf, col_b2)
    ggsave(paste0("Beta2_Bias_", prefijo_archivo, ".pdf"), plot = p_b2_b, width = 11, height = 7, dpi = 300)
    
    p_b2_se <- generar_grafico_panel(df_plot, "Beta2_SE", "Error Estándar (SE) - Beta 2", subtitulo_graf, col_b2)
    ggsave(paste0("Beta2_SE_", prefijo_archivo, ".pdf"), plot = p_b2_se, width = 11, height = 7, dpi = 300)
    
    p_b2_ecm <- generar_grafico_panel(df_plot, "Beta2_ECM", "Error Cuadrático Medio (ECM) - Beta 2", subtitulo_graf, col_b2)
    ggsave(paste0("Beta2_ECM_", prefijo_archivo, ".pdf"), plot = p_b2_ecm, width = 11, height = 7, dpi = 300)
    
    #Graficos Theta
    p_th_b <- generar_grafico_panel(df_plot, "Theta_Bias", "Sesgo (Bias) - Varianza Clúster (Theta)", subtitulo_graf, col_th)
    ggsave(paste0("Theta_Bias_", prefijo_archivo, ".pdf"), plot = p_th_b, width = 11, height = 7, dpi = 300)
    
    p_th_se <- generar_grafico_panel(df_plot, "Theta_SE", "Error Estándar (SE) - Varianza Clúster (Theta)", subtitulo_graf, col_th)
    ggsave(paste0("Theta_SE_", prefijo_archivo, ".pdf"), plot = p_th_se, width = 11, height = 7, dpi = 300)
    
    p_th_ecm <- generar_grafico_panel(df_plot, "Theta_ECM", "Error Cuadrático Medio (ECM) - Varianza Clúster (Theta)", subtitulo_graf, col_th)
    ggsave(paste0("Theta_ECM_", prefijo_archivo, ".pdf"), plot = p_th_ecm, width = 11, height = 7, dpi = 300)
    
  }
}






















