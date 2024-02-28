library(survival)
library(km.ci) # pour les bandes de confiance
library(survminer) # pour les graphiques
library(flexsurv) # pour l'ajustement de modèle
library(lmtest) # pour les test
library(MASS) # pour la sélection de variable avec AIC


#### Données ####
data("diabetic")
str(diabetic)
head(diabetic)

dt <- Surv(diabetic$time, diabetic$status)

#### Ajustement de modèles ####

# Modèle de Weibull
modW1 <- survreg(dt ~ laser + age + eye + trt + risk, data = diabetic,
                dist = "weibull")
summary(modW1)

# test (modèle nul contre modèle complet)
anova(modW1)
lrtest(modW1) # test du rapport de vraisemblance
waldtest(modW1) # test de wald

# sélection de variable
stepAIC(modW1, k = 2) # AIC
stepAIC(modW1, k = log(nrow(diabetic))) # BIC

# modèle final
modW2 <- flexsurvreg(dt ~ trt + risk, dist = "weibull",
                        data = diabetic)
modW2

# comparaison à l'estimateur de Kaplan-Meier
plot(modW2, type = "survival", est = TRUE, ci = TRUE,
     xlab = "Temps", ylab = "Probabiltié de survie",
     main = "Estimation de la fonction de survie dans le modèle de Weibull")
legend(x = "bottomleft", legend = c("Modèle de Weibull", "Kaplan-Meier"),
       title = "Estimateur", col = c("red", "black"), lty = 1, lwd = c(2, 1),
       cex = 0.7)

# résidus de Cox-Snell
resW <- residuals(modW2, type = "coxsnell")
plot(survfit(Surv(resW, status) ~ 1,
             data = diabetic, ctype = 1),
     fun = "cumhaz", xlab = "Résidus", ylab = "Risque cumulé",
     main = expression(paste(
       "Estimateur de Nelson-Aalen de la fonction de risque cumulé")),
     lty = c(1, 2, 2)
)
abline(a = 0, b = 1, col = "red")



# Modèle de log-logistique
modLl1 <- survreg(dt ~ laser + age + eye + trt + risk, data = diabetic,
                dist = "loglogistic")
summary(modLl1)

# test (modèle nul contre modèle complet)
anova(modLl1)
lrtest(modLl1) # test du rapport de vraisemblance
waldtest(modLl1) # test de wald

# sélection de variable
stepAIC(modLl1, k = 2) # AIC
stepAIC(modLl1, k = log(nrow(diabetic))) # BIC

# modèle final
modLl2 <- flexsurvreg(dt ~ trt + risk, dist = "llogis",
                        data = diabetic)
modLl2

# comparaison à l'estimateur de Kaplan-Meier
plot(modLl2, type = "survival", est = TRUE, ci = TRUE,
     main = "Estimation de la fonction de survie dans le modèle log-logistique",
     xlab = "Temps", ylab = "Probabilité de survie", col = "green")
lines(modW2, type = "survival", est = TRUE, ci = FALSE, col = "red", lwd = 1)
legend(x = "bottomleft",legend = c("Modèle de log-logistique",
                                   "Modèle de Weibull",
                                   "Kaplan-Meier"),
       title = "Estimateur", col = c("green", "red", "black"), lty = 1,
       lwd = c(2, 1, 1), cex = 0.7)

# résidus de Cox-Snell
resll <- residuals(modLl2, type = "coxsnell")
plot(survfit(Surv(resll, status) ~ 1, data = diabetic, ctype = 1),
     fun = "cumhaz", xlab = "Résidus", ylab = "Risque cumulé",
     main = expression(paste(
       "Estimateur de Nelson-Aalen de la fonction de risque cumulé")),
     lty = c(1,2,2)
)
abline(a = 0, b = 1, col = "red")

# Modèle de log-normale
modLn1 <- survreg(dt ~ laser + age + eye + trt + risk, data = diabetic,
                  dist = "lognormal")
summary(modLn1)

# test (modèle nul contre modèle complet)
anova(modLn1)
lrtest(modLn1) # test du rapport de vraisemblance
waldtest(modLn1) # test de wald

# sélection de variable
stepAIC(modLn1, k = 2) # AIC
stepAIC(modLn1, k = log(nrow(diabetic))) # BIC

# modèle final
modLn2 <- flexsurvreg(dt ~ trt + risk, dist = "lnorm",
                      data = diabetic)
modLn2

# comparaison à l'estimateur de Kaplan-Meier
plot(modLn2, type = "survival", est = TRUE, ci = TRUE,
     main = "Estimation de la fonction de survie dans le modèle log-normale",
     xlab = "Temps", ylab = "Probabilité", col = "blue")
legend(x = "bottomleft",legend = c("Moddèle log-normale", "Kaplan-Meier"),
       title = "Estimateur", col = c("blue", "black"), lty = 1,
       lwd = c(2, 1), cex = 0.7)

# résidus de Cox-Snell
resln <- residuals(modLn2, type = "coxsnell")

plot(survfit(Surv(resln, status) ~ 1, data = diabetic, ctype = 1),
     fun = "cumhaz", xlab = "Temps", ylab = "Risque cumulé",
     main = expression(paste(
       "Estimateur de Nelson-Aalen de la fonction de risque cumulé")),
     sub = "Modèle de log-logistique", lty = c(1,2)
)
abline(a = 0, b = 1, col = "red")

# Modèle de gompertz
modGp <- flexsurvreg(dt ~ trt + risk, dist = "gompertz",
                      data = diabetic)
modGp

# comparaison à l'estimateur de Kaplan-Meier
plot(modLl2, type = "survival", est = TRUE, ci = TRUE,
     main = "Estimateur de la fonction de survie \n Modèle log-logistique",
     xlab = "Temps", ylab = "Probabilité", sub = "Modèle log-logistique")

# résidus de Cox-Snell
resgp <- residuals(modGp, type = "coxsnell")

plot(survfit(Surv(resgp, status) ~ 1, data = diabetic, ctype = 1),
     fun = "cumhaz", xlab = "Temps", ylab = "Risque cumulé",
     main = expression(paste(
       "Estimateur de Nelson-Aalen de la fonction de risque cumulé")),
     sub = "Modèle gompertz", lty = c(1,2)
)
abline(a = 0, b = 1, col = "red")
