setwd("C:/Users/Bernardo Aroca/OneDrive/Escritorio/tesis/Funciones_pdf_2")
##########
#Exporencial
###########
# Parámetros y etiquetas comunes
p1 <- 0.2; p2 <- 0.5; p3 <- 1.0; p4 <- 1.5
t <- seq(0, 5, length.out = 300)
aa <- expression(lambda==0.2, lambda==0.5, lambda==1.0, lambda==1.5)

# --- G1: Densidad f(t) ---
pdf("exp_f_t.pdf", width = 8, height = 7) # Ajustamos dimensiones del archivo
par(mai=c(0.9, 0.99, 0.05, 0.05)) 
plot(t, dexp(t, p1), type="l", col="blue", lwd=2, ylim=c(0, 1.6), 
     ylab="f(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, dexp(t, p2), col="red", lwd=2)
lines(t, dexp(t, p3), col="darkgreen", lwd=2)
lines(t, dexp(t, p4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), 
       lwd=2, bty="n", cex=2.5)
grid()
dev.off()

# --- G2: Sobrevivencia S(t) ---
pdf("exp_S_t.pdf", width = 8, height = 7)
par(mai=c(0.9, 0.99, 0.05, 0.05)) 
plot(t, pexp(t, p1, lower.tail=FALSE), type="l", col="blue", lwd=2, 
     ylim=c(0, 1), ylab="S(t)", xlab="Tiempo", 
     cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, pexp(t, p2, lower.tail=FALSE), col="red", lwd=2)
lines(t, pexp(t, p3, lower.tail=FALSE), col="darkgreen", lwd=2)
lines(t, pexp(t, p4, lower.tail=FALSE), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), 
       lwd=2, bty="n", cex=2.5)
grid()
dev.off()

# --- G3: Riesgo λ(t) ---
pdf("exp_h_t.pdf", width = 8, height = 7)
par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, rep(p1, length(t)), type="l", col="blue", lwd=2, 
     ylim=c(0, 2), ylab=expression(lambda(t)), xlab="Tiempo", 
     cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, rep(p2, length(t)), col="red", lwd=2)
lines(t, rep(p3, length(t)), col="darkgreen", lwd=2)
lines(t, rep(p4, length(t)), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), 
       lwd=2, bty="n", cex=2.5)
grid()
dev.off()
##########
#Weibull
###########
##########
# Weibull con dos parámetros
##########
# Parámetros (nu = forma, lambda = escala)
nu1 <- 0.5; l1 <- 1.0; nu2 <- 1.0; l2 <- 1.0
nu3 <- 1.5; l3 <- 1.0; nu4 <- 3.0; l4 <- 2.0
t <- seq(0.01, 5, length.out = 300)

aa <- expression(paste(nu==0.5, ", ", lambda==1.0), paste(nu==1.0, ", ", lambda==1.0),
                 paste(nu==1.5, ", ", lambda==1.0), paste(nu==3.0, ", ", lambda==2.0))

# --- G1: Densidad ---
pdf("weibull_f_t.pdf", width = 8, height = 7)
par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, dweibull(t, nu1, l1), type="l", col="blue", lwd=2, ylim=c(0, 1.5), 
     ylab="f(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, dweibull(t, nu2, l2), col="red", lwd=2)
lines(t, dweibull(t, nu3, l3), col="darkgreen", lwd=2)
lines(t, dweibull(t, nu4, l4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5)
grid(); dev.off()

# --- G2: Sobrevivencia ---
pdf("weibull_S_t.pdf", width = 8, height = 7)
par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, pweibull(t, nu1, l1, lower.tail=F), type="l", col="blue", lwd=2, ylim=c(0, 1), 
     ylab="S(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, pweibull(t, nu2, l2, lower.tail=F), col="red", lwd=2)
lines(t, pweibull(t, nu3, l3, lower.tail=F), col="darkgreen", lwd=2)
lines(t, pweibull(t, nu4, l4, lower.tail=F), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5)
grid(); dev.off()

# --- G3: Riesgo ---
h_weib <- function(t, nu, lam) { dweibull(t, nu, lam) / pweibull(t, nu, lam, lower.tail=F) }
pdf("weibull_h_t.pdf", width = 8, height = 7)
par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, h_weib(t, nu1, l1), type="l", col="blue", lwd=2, ylim=c(0, 4), 
     ylab="h(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, h_weib(t, nu2, l2), col="red", lwd=2)
lines(t, h_weib(t, nu3, l3), col="darkgreen", lwd=2)
lines(t, h_weib(t, nu4, l4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5)
grid(); dev.off()
##########
#Gompertz
###########
##########
# Gompertz con dos parámetros
##########
library(flexsurv)
e1 <- 0.2; r1 <- 0.1; e2 <- 0.4; r2 <- 0.1; e3 <- 0.7; r3 <- 0.1; e4 <- 1.0; r4 <- 0.3
t <- seq(0, 8, length.out = 300)
aa <- expression(paste(varphi==0.2, ", ", lambda==0.1), paste(varphi==0.4, ", ", lambda==0.1),
                 paste(varphi==0.7, ", ", lambda==0.1), paste(varphi==1.0, ", ", lambda==0.3))

# --- G1: Densidad ---
pdf("gompertz_f_t.pdf", width = 8, height = 7)
par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, dgompertz(t, e1, r1), type="l", col="blue", lwd=2, ylim=c(0, 0.6), ylab="f(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, dgompertz(t, e2, r2), col="red", lwd=2)
lines(t, dgompertz(t, e3, r3), col="darkgreen", lwd=2)
lines(t, dgompertz(t, e4, r4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5)
grid(); dev.off()

# --- G2: Sobrevivencia ---
pdf("gompertz_S_t.pdf", width = 8, height = 7)
par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, pgompertz(t, e1, r1, lower.tail=F), type="l", col="blue", lwd=2, ylim=c(0, 1), ylab="S(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, pgompertz(t, e2, r2, lower.tail=F), col="red", lwd=2)
lines(t, pgompertz(t, e3, r3, lower.tail=F), col="darkgreen", lwd=2)
lines(t, pgompertz(t, e4, r4, lower.tail=F), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5)
grid(); dev.off()

# --- G3: Riesgo ---
pdf("gompertz_h_t.pdf", width = 8, height = 7)
par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, hgompertz(t, e1, r1), type="l", col="blue", lwd=2, ylim=c(0, 8), ylab="h(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, hgompertz(t, e2, r2), col="red", lwd=2)
lines(t, hgompertz(t, e3, r3), col="darkgreen", lwd=2)
lines(t, hgompertz(t, e4, r4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5)
grid(); dev.off()
#############
# MODELOS PARAMETRICOS ALTERNATIVAS
##########
##########
#Gamma
###########
# Parámetros y etiquetas
s1 <- 0.5; r1 <- 1.0; s2 <- 1.5; r2 <- 1.0; s3 <- 3.0; r3 <- 1.0; s4 <- 5.0; r4 <- 2.0
t <- seq(0, 8, length.out = 300)
aa <- expression(paste(beta==0.5, ", ", lambda==1.0), paste(beta==1.5, ", ", lambda==1.0),
                 paste(beta==3.0, ", ", lambda==1.0), paste(beta==5.0, ", ", lambda==2.0))

# --- G1: Densidad ---
pdf("gamma_f_t.pdf", width = 8, height = 7); par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, dgamma(t, s1, r1), type="l", col="blue", lwd=2, ylim=c(0, 1), ylab="f(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, dgamma(t, s2, r2), col="red", lwd=2); lines(t, dgamma(t, s3, r3), col="darkgreen", lwd=2); lines(t, dgamma(t, s4, r4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid()
dev.off()

# --- G2: Sobrevivencia ---
pdf("gamma_S_t.pdf", width = 8, height = 7); par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, pgamma(t, s1, r1, lower.tail=F), type="l", col="blue", lwd=2, ylim=c(0, 1), ylab="S(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, pgamma(t, s2, r2, lower.tail=F), col="red", lwd=2); lines(t, pgamma(t, s3, r3, lower.tail=F), col="darkgreen", lwd=2); lines(t, pgamma(t, s4, r4, lower.tail=F), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid()
dev.off()

# --- G3: Riesgo ---
pdf("gamma_h_t.pdf", width = 8, height = 7); par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, dgamma(t, s1, r1)/pgamma(t, s1, r1, lower.tail=F), type="l", col="blue", lwd=2, ylim=c(0, 2.5), ylab=expression(lambda(t)), xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, dgamma(t, s2, r2)/pgamma(t, s2, r2, lower.tail=F), col="red", lwd=2); lines(t, dgamma(t, s3, r3)/pgamma(t, s3, r3, lower.tail=F), col="darkgreen", lwd=2); lines(t, dgamma(t, s4, r4)/pgamma(t, s4, r4, lower.tail=F), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid()
dev.off()

##########
#Distribución Gaussiana Inversa (Wald)
###########
if(!require(statmod)) install.packages("statmod")
library(statmod)

m1 <- 1.0; s1 <- 0.2; m2 <- 1.0; s2 <- 1.0; m3 <- 1.5; s3 <- 2.0; m4 <- 2.0; s4 <- 0.5
t <- seq(0.01, 5, length.out = 300)
aa <- expression(paste(mu==1.0, ", ", lambda==0.2), paste(mu==1.0, ", ", lambda==1.0),
                 paste(mu==1.5, ", ", lambda==2.0), paste(mu==2.0, ", ", lambda==0.5))

# --- G1: Densidad ---
pdf("wald_f_t.pdf", width = 8, height = 7); par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, dinvgauss(t, m1, s1), type="l", col="blue", lwd=2, ylim=c(0, 1.5), ylab="f(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, dinvgauss(t, m2, s2), col="red", lwd=2); lines(t, dinvgauss(t, m3, s3), col="darkgreen", lwd=2); lines(t, dinvgauss(t, m4, s4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid()
dev.off()

# --- G2: Sobrevivencia ---
pdf("wald_S_t.pdf", width = 8, height = 7)
par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, pinvgauss(t, m1, s1, lower.tail=F), type="l", col="blue", lwd=2, ylim=c(0, 1), ylab="S(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, pinvgauss(t, m2, s2, lower.tail=F), col="red", lwd=2); lines(t, pinvgauss(t, m3, s3, lower.tail=F), col="darkgreen", lwd=2); lines(t, pinvgauss(t, m4, s4, lower.tail=F), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid()
dev.off()

# --- G3: Riesgo ---
pdf("wald_h_t.pdf", width = 8, height = 7)
par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, dinvgauss(t, m1, s1)/pinvgauss(t, m1, s1, lower.tail=F), type="l", col="blue", lwd=2, ylim=c(0, 2), ylab=expression(lambda(t)), xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, dinvgauss(t, m2, s2)/pinvgauss(t, m2, s2, lower.tail=F), col="red", lwd=2); lines(t, dinvgauss(t, m3, s3)/pinvgauss(t, m3, s3, lower.tail=F), col="darkgreen", lwd=2); lines(t, dinvgauss(t, m4, s4)/pinvgauss(t, m4, s4, lower.tail=F), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid()
dev.off()

##########
#Distribución Birnbaum-Saunders
##########
f_bs <- function(t, a, b) { z <- (sqrt(t/b) - sqrt(b/t)) / a; return(((sqrt(t/b) + sqrt(b/t)) / (2 * a * t)) * dnorm(z)) }
s_bs <- function(t, a, b) { z <- (sqrt(t/b) - sqrt(b/t)) / a; return(1 - pnorm(z)) }
a1 <- 0.5; b1 <- 1.0; a2 <- 1.0; b2 <- 1.0; a3 <- 1.5; b3 <- 1.0; a4 <- 0.5; b4 <- 2.0
t <- seq(0.01, 6, length.out = 300)
aa <- expression(paste(alpha==0.5, ", ", beta==1.0), paste(alpha==1.0, ", ", beta==1.0),
                 paste(alpha==1.5, ", ", beta==1.0), paste(alpha==0.5, ", ", beta==2.0))

# --- G1, G2, G3 ---
pdf("bs_f_t.pdf", width = 8, height = 7); par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, f_bs(t, a1, b1), type="l", col="blue", lwd=2, ylim=c(0, 1.2), ylab="f(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, f_bs(t, a2, b2), col="red", lwd=2); lines(t, f_bs(t, a3, b3), col="darkgreen", lwd=2); lines(t, f_bs(t, a4, b4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid()
dev.off()

pdf("bs_S_t.pdf", width = 8, height = 7); par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, s_bs(t, a1, b1), type="l", col="blue", lwd=2, ylim=c(0, 1), ylab="S(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, s_bs(t, a2, b2), col="red", lwd=2); lines(t, s_bs(t, a3, b3), col="darkgreen", lwd=2); lines(t, s_bs(t, a4, b4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid(); dev.off()

pdf("bs_h_t.pdf", width = 8, height = 7); par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, f_bs(t, a1, b1)/s_bs(t, a1, b1), type="l", col="blue", lwd=2, ylim=c(0, 1.5), ylab=expression(lambda(t)), xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, f_bs(t, a2, b2)/s_bs(t, a2, b2), col="red", lwd=2); lines(t, f_bs(t, a3, b3)/s_bs(t, a3, b3), col="darkgreen", lwd=2); lines(t, f_bs(t, a4, b4)/s_bs(t, a4, b4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid()
dev.off()
##########
# Distribución Lindley Ponderada (Weighted Lindley)
##########
f_wl <- function(t, theta, alpha) {
  num <- theta^2 * (1 + t) * exp(-theta * t)
  den <- (1 + theta) * (1 - alpha + alpha * (1 + theta + theta * t))
  return(num / den)
}
s_wl <- function(t, theta, alpha) {
  return(1 - (1 - exp(-theta * t) * (1 + theta * t/(1 + theta))) * (1 - alpha))
}

a1 <- 0.5; w1 <- 0.2; a2 <- 1.0; w2 <- 0.5; a3 <- 0.5; w3 <- 0.8; a4 <- 1.5; w4 <- 0.3
t <- seq(0, 6, length.out = 300)
aa <- expression(paste(alpha==0.5, ", ", omega==0.2), paste(alpha==1.0, ", ", omega==0.5),
                 paste(alpha==0.5, ", ", omega==0.8), paste(alpha==1.5, ", ", omega==0.3))

# --- G1, G2, G3 ---
pdf("wl_f_t.pdf", width = 8, height = 7); par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, f_wl(t, a1, w1), type="l", col="blue", lwd=2, ylim=c(0, 0.8), ylab="f(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, f_wl(t, a2, w2), col="red", lwd=2); lines(t, f_wl(t, a3, w3), col="darkgreen", lwd=2); lines(t, f_wl(t, a4, w4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid(); dev.off()

pdf("wl_S_t.pdf", width = 8, height = 7); par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, s_wl(t, a1, w1), type="l", col="blue", lwd=2, ylim=c(0, 1), ylab="S(t)", xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, s_wl(t, a2, w2), col="red", lwd=2); lines(t, s_wl(t, a3, w3), col="darkgreen", lwd=2); lines(t, s_wl(t, a4, w4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid(); dev.off()

pdf("wl_h_t.pdf", width = 8, height = 7); par(mai=c(0.9, 0.99, 0.05, 0.05))
plot(t, f_wl(t, a1, w1)/s_wl(t, a1, w1), type="l", col="blue", lwd=2, ylim=c(0, 2), ylab=expression(lambda(t)), xlab="Tiempo", cex.lab=1.8, cex.axis=1.4, las=1)
lines(t, f_wl(t, a2, w2)/s_wl(t, a2, w2), col="red", lwd=2); lines(t, f_wl(t, a3, w3)/s_wl(t, a3, w3), col="darkgreen", lwd=2); lines(t, f_wl(t, a4, w4)/s_wl(t, a4, w4), col="orange", lwd=2)
legend("topright", legend=aa, col=c("blue","red","darkgreen","orange"), lwd=2, bty="n", cex=2.5); grid(); dev.off()

