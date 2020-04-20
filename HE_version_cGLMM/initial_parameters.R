# initial parameters

num_fe <- 10
epsilon <- 1e10
threshold <- 0.05
beta <- rep(1,10)
sigma1 <- 0.5
K <- 2
burnin <- 198
n_user <- 2
n = 2 # sampling size for Metropolis Hasting

Add_Big_Number = 10000
Precision = 10000
save(n_user, Add_Big_Number, Precision, num_fe, epsilon, threshold, beta, sigma1, K, burnin, n, file="initial_parameters.RData")