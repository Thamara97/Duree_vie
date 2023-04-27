library(smoothmest)

n <- 8 ; m <- 13

ps <- c() ; pw <- c() ; pks <- c()

j <- 0

# theta
for (theta in seq(0, 5, by = 0.5)) {
  
  j <- j + 1
 
  # nbr p-value < 0.05
  nb_ps <- 0 ; nb_pw <- 0 ; nb_pks <- 0

  for (i in 1:10000) {
    
    # Echantillon 
    X <- rdoublex(n, mu = 0, lambda = 1) # échantillon X
    Y <- rdoublex(m, mu = theta, lambda = 1) # échantillon Y
      
    if (shapiro.test(X)$p.value > 0.05 | shapiro.test(Y)$p.value > 0.05) {
      # Test de Student
      Stud_test <- t.test(X, Y, alternative = "less")
      if (Stud_test$p.value < 0.05) {nb_ps <- nb_ps + 1} 
        
      # Test de Wilcoxon
      Wil_test <- wilcox.test(X, Y, alternative = "less")
      if (Wil_test$p.value < 0.05) {nb_pw <- nb_pw + 1}
      
      # Test de Kolmogorov-Smirnov
      KS_test <- ks.test(X, Y, alternative = "two.sided")
      if (KS_test$p.value < 0.05) {nb_pks <- nb_pks + 1}
    }
    
  }
  
  ps[j] <- nb_ps/10000 # 
  pw[j] <- nb_pw/10000
  pks[j] <- nb_pks/10000
}

plot(seq(0, 5, by = 0.5), ps, type = "o", col = "blue", xlab = "", ylab = "Puissance")
points(seq(0, 5, by = 0.5), pw, type = "o", col = "red")
points(seq(0, 5, by = 0.5), pks, type = "o", col = "green")


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

plot(seq(0, 5, by = 0.5), C, type = "o")
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

plot(seq(0, 5, by = 0.5), C, type = "o")
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

plot(seq(0, 5, by = 0.5), C, type = "o")
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

# d'après nos tests on a n = 19

n = 19
