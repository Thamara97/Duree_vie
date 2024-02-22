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

# test (modèle null contre modèle complet)
anova(modW1)
lrtest(modW1) # test du rapport de vraisemblance
waldtest(modW1) # test de wald

# sélection de variable
stepAIC(modW1, k = 2) # AIC
stepAIC(modW1, k = log(nrow(diabetic))) # BIC

# modèle final
modW2 <- flexsurvreg(dt ~ eye + trt + risk, dist = "weibull",
                        data = diabetic)
modW2

# comparaison à l'estimateur de Kaplan-Meier
plot(modW2, type = "survival", est = TRUE, ci = TRUE)

# résidus de Cox-Snell
gp <- ifelse(diabetic$eye == "right", 1, 0)

Rw <- (exp(-(log(466.4812) - 0.4287 * gp + 1.0166 * diabetic$trt -
               0.1743 * diabetic$risk)) * diabetic$time)^0.8181

plot(survfit(Surv(Rw, status) ~ 1, data = diabetic, ctype = 1),
     fun = "cumhaz", xlab = "Temps", ylab = "Risque cumulé",
     main = expression(paste(
       "Estimateur de Nelson-Aalen de la fonction de risque cumulé")),
     sub = "Modèle de Weibull", lty = c(1,2)
)
abline(a = 0, b = 1, col = "red")

# Modèle de log-logistique
modLl1 <- survreg(dt ~ laser + age + eye + trt + risk, data = diabetic,
                dist = "loglogistic")
summary(modLl1)

# test (modèle null contre modèle complet)
anova(modLl1)
lrtest(modLl1) # test du rapport de vraisemblance
waldtest(modLl1) # test de wald

# sélection de variable
stepAIC(modLl1, k = 2) # AIC
stepAIC(modLl1, k = log(nrow(diabetic))) # BIC

# modèle final
modLl2 <- flexsurvreg(dt ~ eye + trt + risk, dist = "llogis",
                        data = diabetic)
modLl2

# comparaison à l'estimateur de Kaplan-Meier
plot(modLl2, type = "survival", est = TRUE, ci = TRUE,
     main = "Estimateur de la fonction de survie \n Modèle log-logistique",
     xlab = "Temps", ylab = "Probabilité", sub = "Modèle log-logistique")

# résidus de Cox-Snell
Rll <- log(1 + (exp(-(log(380.6968) -0.5072 * gp + 1.0777 * diabetic$trt
          -0.1979 * diabetic$risk)) * diabetic$time)^0.9600)

plot(survfit(Surv(Rll, status) ~ 1, data = diabetic, ctype = 1),
     fun = "cumhaz", xlab = "Temps", ylab = "Risque cumulé",
     main = expression(paste(
       "Estimateur de Nelson-Aalen de la fonction de risque cumulé")),
     sub = "Modèle de log-logistique", lty = c(1,2)
)
abline(a = 0, b = 1, col = "red")

# Modèle de log-normale
modLn1 <- survreg(dt ~ laser + age + eye + trt + risk, data = diabetic,
                  dist = "lognormal")
summary(modLn1)

# test (modèle null contre modèle complet)
anova(modLn1)
lrtest(modLn1) # test du rapport de vraisemblance
waldtest(modLn1) # test de wald

# sélection de variable
stepAIC(modLn1, k = 2) # AIC
stepAIC(modLn1, k = log(nrow(diabetic))) # BIC

# modèle final
modLn2 <- flexsurvreg(dt ~ eye + trt + risk, dist = "lnorm",
                      data = diabetic)
modLn2

# comparaison à l'estimateur de Kaplan-Meier
plot(modLn2, type = "survival", est = TRUE, ci = TRUE,
     main = "Estimateur de la fonction de survie \n Modèle log-normale",
     xlab = "Temps", ylab = "Probabilité")

# résidus de Cox-Snell
modLn3 <- survreg(dt ~ eye + trt + risk, data = diabetic, dist = "lognormal")
Rln <- - log(1 - pnorm((log(diabetic$time) - modLn3$coefficients[1] -
                          modLn3$coefficients[2] * gp -
                          modLn3$coefficients[3] * diabetic$trt -
                          modLn3$coefficients[4] * diabetic$risk) /
                         modLn3$scale))

plot(survfit(Surv(Rln, status) ~ 1, data = diabetic, ctype = 1),
     fun = "cumhaz", xlab = "Temps", ylab = "Risque cumulé",
     main = expression(paste(
       "Estimateur de Nelson-Aalen de la fonction de risque cumulé")),
     sub = "Modèle de log-logistique", lty = c(1,2)
)
abline(a = 0, b = 1, col = "red")
