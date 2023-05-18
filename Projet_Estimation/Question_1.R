library(smoothmest)

# Taille des échantillions
n <- 8 ; m <- 13

# Vecteurs des puissances empirique
ps <- c() ; pw <- c() ; pks <- c()

# Test numéro j
j <- 0

for (theta in seq(0, 4, by = 0.25)) {
  
  j <- j + 1
 
  # Nombre de p-value < 0.05
  nb_ps <- 0 ; nb_pw <- 0 ; nb_pks <- 0
  
  # Nombre de tests
  u <- 0

  for (i in 1:10000) {
    
    # Echantillons
    X <- rdoublex(n, mu = 0, lambda = 1)
    Y <- rdoublex(m, mu = theta, lambda = 1)
      
    # Normalité des données
    if (shapiro.test(X)$p.value > 0.05 | shapiro.test(Y)$p.value > 0.05) {
      
      # Test de Student
      if (bartlett.test(c(X,Y)~c(rep(1,n), rep(2,m)))$p.value > 0.05) {
        Stud_test <- t.test(X, Y, alternative = "less", var.equal = TRUE)
      } else {
          Stud_test <- t.test(X, Y, alternative = "less")
          } 
      
      if (Stud_test$p.value < 0.05) {nb_ps <- nb_ps + 1} 
        
      # Test de Wilcoxon
      Wil_test <- wilcox.test(X, Y, alternative = "less")
      if (Wil_test$p.value < 0.05) {nb_pw <- nb_pw + 1}
      
      # Test de Kolmogorov-Smirnov
      KS_test <- ks.test(X, Y, alternative = "two.sided")
      if (KS_test$p.value < 0.05) {nb_pks <- nb_pks + 1}
      
      u <- u + 1
    }
    
  }
  
  ps[j] <- nb_ps/u
  pw[j] <- nb_pw/u
  pks[j] <- nb_pks/u
}

plot(seq(0, 4, by = 0.25), ps, type = "o", col = "blue", xlab = expression(theta), ylab = "Puissance")
points(seq(0, 4, by = 0.25), pw, type = "o", col = "red")
points(seq(0, 4, by = 0.25), pks, type = "o", col = "green")


#QUESTION 2 

library(wmwpow)
library(smoothmest)

# déterminer la valeur de theta*

p = 0.7    # p = f(theta)
N = 10000
C = c()

for (theta in seq(0, 5, by = 0.5)) {
  c = 0
    for(i in (1:100)){
       X <- rdoublex(N, mu = 0, lambda = 1) # échantillon X
       Y <- rdoublex(N, mu = theta, lambda = 1) # échantillon Y
       c = c(c,sum(X<=Y)/N)
    }
  C = c(C,mean(c))
}

plot(seq(0, 5, by = 0.5), C, type = "o")
abline(h = p)

# maintenant on sait que theta est entre 0.5 et 1, on refait alors le même calcul pour theta allant de 0 à 1

C = c()

for (theta in seq(0, 1, by = 0.1)) {
  c = 0
    for(i in (1:100)){
       X <- rdoublex(N, mu = 0, lambda = 1) # échantillon X
       Y <- rdoublex(N, mu = theta, lambda = 1) # échantillon Y
       c = c(c,sum(X<=Y)/N)
    }
  C = c(C,mean(c))
}

C

plot(seq(0, 1, by = 0.1), C, type = "o")
abline(h = p)

#maintenant on sait que theta est entre 0.9 et 1

C = c()

for (theta in seq(0.9, 1, by = 0.01)) {
  c = 0
    for(i in (1:100)){
       X <- rdoublex(N, mu = 0, lambda = 1) # échantillon X
       Y <- rdoublex(N, mu = theta, lambda = 1) # échantillon Y
       c = c(c,sum(X<=Y)/N)
    }
  C = c(C,mean(c))
}

C

plot(seq(0.9, 1, by = 0.01), C, type = "o")
abline(h = p)

# on sait maintenant que theta est entre 0.90 et 0.91
C = c()

for (theta in seq(0.9, 0.91, by = 0.001)) {
  c = 0
    for(i in (1:100)){
       X <- rdoublex(N, mu = 0, lambda = 1) # échantillon X
       Y <- rdoublex(N, mu = theta, lambda = 1) # échantillon Y
       c = c(c,sum(X<=Y)/N)
    }
  C = c(C,mean(c))
}

C

plot(seq(0.9, 0.91, by = 0.001), C, type = "o")
abline(h = p)

# on approxime alors la valeur de theta à 0.91

theta_e = 0.91

# tests pour trouver n = m

pow = c()
for(n in (1:10)){
  powpow = wmwpowp(n,n,"doublex", p= 0.7, sides = "less", alpha = 0.05)
  pow = c(pow, powpow$empirical_power)
}

plot(1:10, pow, type = "o")
abline(h = 0.8)

# on sait que n> 10

pow = c()
for(n in (11:20)){
  powpow = wmwpowp(n,n,"doublex", p= 0.7, sides = "less", alpha = 0.05)
  pow = c(pow, powpow$empirical_power)
}

plot(11:20, pow, type = "o")
abline(h = 0.8)

# d'après nos tests n > 20

pow = c()
for(n in (21:30)){
  powpow = wmwpowp(n,n,"doublex", p= 0.7, sides = "less", alpha = 0.05)
  pow = c(pow, powpow$empirical_power)
}

# on a alors n = 25

n = 25

# QUESTION 3

# Données
X <- c(-0.5, 0., 3.1, 0.8, 1.1, -0.4, -0.5, -0.2, 0.3, -0.2, -0.2, 0.5, 1.3, 0.7, 1.6, -0.2, 1.7, 1., 0.5, -0.2, -0.7, 1.4, -0.6, 2.1, 0.6, -0.2)

Y <- c(1.4, 3.9, 1.1, 1.6, 1.2, 0.2, 1.1, -0.5, 1.8, 0.6, -0.7, 3.7, -1.7, 0., 1.9, -0.2, 1., 0.4, 2.7, 1.1, 1.2, 1.5, -1.6, 0.6, 0.7, -0.1)

bartlett.test(c(X, Y) ~ c(rep(1, 26), rep(2, 26)))
# même variance

# Test
wilcox.test(X, Y, alternative = "less")
# Non rejet de H_0, on l'accepte au risque de 20%
