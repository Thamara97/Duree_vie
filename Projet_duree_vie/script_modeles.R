library(survival)
library(survminer) # graphiques
library(flexsurv) # ajustement de modèle
library(lmtest) # test
library(MASS) # sélection de variable


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

# modèle final (BIC)
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
     main = "Estimateur de Nelson-Aalen de la fonction de risque cumulé",
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
legend(x = "bottomleft",legend = c("Modèle de log-logistique", "Kaplan-Meier"),
       title = "Estimateur", col = c("green", "black"), lty = 1,
       lwd = c(2, 1, 1), cex = 0.7)

# résidus de Cox-Snell
resll <- residuals(modLl2, type = "coxsnell")
plot(survfit(Surv(resll, status) ~ 1, data = diabetic, ctype = 1),
     fun = "cumhaz", xlab = "Résidus", ylab = "Risque cumulé",
     main = "Estimateur de Nelson-Aalen de la fonction de risque cumulé",
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

# modèle final (BIC)
modLn2 <- flexsurvreg(dt ~ trt + risk, dist = "lnorm",
                      data = diabetic)
modLn2

# comparaison à l'estimateur de Kaplan-Meier
plot(modLn2, type = "survival", est = TRUE, ci = FALSE,
     main = "Estimation de la fonction de survie dans le modèle log-normale",
     xlab = "Temps", ylab = "Probabilité", col = "blue")
legend(x = "bottomleft",legend = c("Moddèle log-normale","Kaplan-Meier"),
       title = "Estimateur", col = c("blue", "black"), lty = 1,
       lwd = c(2, 1, 1, 1), cex = 0.7)

# résidus de Cox-Snell
resln <- residuals(modLn2, type = "coxsnell")
plot(survfit(Surv(resln, status) ~ 1, data = diabetic, ctype = 1),
     fun = "cumhaz", xlab = "Temps", ylab = "Risque cumulé",
     main = "Estimateur de Nelson-Aalen de la fonction de risque cumulé",
     lty = c(1,2)
)
abline(a = 0, b = 1, col = "red")


#### Comparaison des modèles ####

# Kaplan-Meier
plot(modLn2, type = "survival", est = TRUE, ci = TRUE,
     main = "Estimation de la fonction dans les différents modèles",
     xlab = "Temps", ylab = "Probabilité", col = "blue")
lines(modLl2, type = "survival", est = TRUE, ci = T, col = "green")
lines(modW2, type = "survival", est = TRUE, ci = T)
legend(x = "bottomleft",legend = c("Moddèle log-normale",
                                   "Modèle de log-logistique",
                                   "Modèle de Weibull", "Kaplan-Meier"),
       title = "Estimateur", col = c("blue", "green", "red", "black"), lty = 1,
       lwd = c(2, 2, 2, 1), cex = 0.7)

#### Ajustement du modèle de Cox ####
data <- diabetic
data$trt <- relevel(as.factor(data$trt), ref = "1")

modCox <- coxph(dt ~ trt + risk, data = data)
summary(modCox)

ggforest(modCox, data = data)

#### Vérification de l'hypothèse des risques proportionnels ####

# Résidus de Schöenfeld
residus <- cox.zph(modCox)
ggcoxzph(residus)

# Graphique "LML" pour trt
indiv_new <- data.frame(trt = c("1", "0"), risk = c(9, 9))
survie_pred <- survfit(modCox, newdata = indiv_new, conf.type = "none")
plot(survie_pred, fun = "cloglog", col = c("black", "red"),
     xlab = "Temps", ylab = "ln(-ln(S(t|Z)))",
     main = "Courbes \"LML\" pour les modalités de la covariable trt
     quand risk = 9")

#### Prédiction de la fonction de survie ####
indiv_new <- data.frame(trt = rep(c("0", "1"), each = 6),
                        risk = rep(6:12, times = 2)[-c(2, 9)])
surv.new <- survfit(modCox, newdata = indiv_new)
leng <- paste("trt =", rep(c("0", "1"), each = 6), "risk =", c(6:12)[-2])
med <- c(63.3, 48.9, 42.4, 34.4, 26.2)

plot(surv.new, col = 1:14, lty = rep(c(1,2), each = 6),
     xlab = "Temps", ylab = "Probabilité de survie",
     main = "Prédiction de la fonction de survie pour plusieurs individus")
legend(x = "bottomleft", legend = leng, col = 1:12, lty = rep(c(1,2), each = 6),
       lwd = 1, cex = 0.7)
for (i in 2:6) {
  lines(x = c(-2.5, med[i-1]), y = c(0.5, 0.5), lty = 3, col = i, lwd = 2)
  lines(x = c(med[i-1], med[i-1]), y = c(-0.5, 0.5), lty = 3, col = i, lwd = 2)
}
axis(side = 2, at = seq(0, 1, by = 0.1))
