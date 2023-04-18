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
  
  ps[j] <- nb_ps/10000 # proportion de p-value Student
  pw[j] <- nb_pw/10000 # proportion de p-value Wilcoxon
  pks[j] <- nb_pks/10000 # proportion de p-value Kolmogorov
}

plot(seq(0, 5, by = 0.5), ps, type = "o", col = "blue", xlab = "", ylab = "Puissance")
points(seq(0, 5, by = 0.5), pw, type = "o", col = "red")
points(seq(0, 5, by = 0.5), pks, type = "o", col = "green")
