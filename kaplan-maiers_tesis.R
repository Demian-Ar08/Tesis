##Antiproyecto
install.packages("timereg")
library(timereg)
data("melanoma")
attach(melanoma)
t1<-days/365.25
d1<-ifelse(status==1,1,0)

Surv(t1,d1)
KM<-survfit(Surv(t1,d1) ~ 1,conf.type="log-log")
summary(KM)

# Gráfico limpio con IC 95% y línea de mediana
plot(KM,
     main = "Curva de Supervivencia de Kaplan-Meier\nPacientes con Melanoma",
     xlab = "Tiempo (años)",
     ylab = "Probabilidad de supervivencia",
     col = "black",
     lwd = 2,
     conf.int = TRUE,           # Intervalos de confianza al 95%
     mark.time = TRUE,          # Marcas de censura
     xlim = c(0, max(t1, na.rm = TRUE)),
     ylim = c(0, 1))# Línea horizontal en 0.5 (marca la mediana)
abline(h = 0.5, col = "red", lty = 2, lwd = 1.5)

# Cuadrícula suave (opcional, ayuda a la lectura)
grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")

# Leyenda para explicar los elementos
legend("bottomleft",
       legend = c("Curva de supervivencia", "IC 95%", "Mediana (0.5)"),
       col = c("black", "gray", "red"),
       lwd = c(2, 2, 1.5),
       lty = c(1, 1, 2),
       bty = "n")


KM_sex<-survfit(Surv(t1,d1) ~ sex,conf.type="log-log")
plot(KM_sex)
