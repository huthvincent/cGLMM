setwd("/home/huthvincent/Desktop/glmm/150/Fed_GLMM")
load("/home/huthvincent/Desktop/glmm/150/Fed_GLMM/X.Rdata")
load("/home/huthvincent/Desktop/glmm/150/Fed_GLMM/Y.Rdata")
load("/home/huthvincent/Desktop/glmm/150/Fed_GLMM/Z.Rdata")

# calcuate the real_parameters by taking the Z1 as know in advance
options(digits = 15)
real_parameters <- function(X, Y, Z1, initial_beta = rnorm(50,0,0.01)) {
  num_fe = dim(X[[1]])[2]
  
  K = length(Z1)
  #real_sigma1 就是 Z1的
  real_sigma1 <- sqrt(mean(Z1 ^ 2))
  real_beta = initial_beta
  epi <- 2
  
  #这里的f2就是likelihood l 对 beta的二阶导
  while (epi > 1e-5) {
    H = matrix(0, num_fe, num_fe)
    for (i in c(1:num_fe)) {
      for (j in c(1:num_fe)) {
        f2 = rep(0, K)
        for (k in c(1:K)) {
          num_pat = dim(X[[k]])[1]
          f2[k] = sum(X[[k]][, i] * X[[k]][, j] * exp(X[[k]] %*% real_beta  + rep(Z1[k], num_pat)) /
                        (1 + exp(
                          X[[k]] %*% real_beta  + rep(Z1[k], num_pat)
                        )) ^ 2)
        }
        H[i, j] = sum(f2)
      }
    }
    
    H_inverse = solve(H)
    
    f1 = rep(0, num_fe)
    for (i in c(1:num_fe)) {
      temp = rep(0, K)
      for (k in c(1:K)) {
        num_pat = dim(X[[k]])[1]
        temp[k] = sum(X[[k]][, i] * (Y[[k]] - 1 + 1 / (1 + exp(
          X[[k]] %*% real_beta  + rep(Z1[k], num_pat)
        ))))
      }
      f1[i] = sum(temp)
    }
    
    updata = H_inverse %*% f1
    
    new_beta = real_beta + updata
    epi = max(abs(new_beta - real_beta))
    real_beta <- new_beta
    
  }
  
  return(list(real_beta, real_sigma1))
}

true_value = real_parameters(X, Y, Z)
true_value
true_beta = true_value[[1]]
true_sigma = true_value[[2]]
save.image(file="cur_session.RData")
save(true_beta,file="ground_truth_beta.RData")
save(true_sigma,file="ground_truth_sigma.RData")
# method 1
# cut data into two parts both has 10 random effects 
K = 10
m = 50
X_1 = rep(list(matrix(0, 1, m)), K)
X_2 = rep(list(matrix(0, 1, m)), K)

Y_1 = rep(list(matrix(0, 1, 1)), K)
Y_2 = rep(list(matrix(0, 1, 1)), K)


for (i in c(1:K)) {
  rows = length(unlist(X[i]))/m
  tmp_X = matrix(unlist(X[i]), rows , m)
  tmp_Y = matrix(unlist(Y[i]), rows , 1)
  first_part = rows%/%2
  second_part = first_part + 1
  X_1[[i]] = tmp_X[1:first_part,]
  Y_1[[i]] = tmp_Y[1:first_part,]
  #cat("first_part",first_part, "rows", rows)
  X_2[[i]] = tmp_X[second_part:rows,]
  Y_2[[i]] = tmp_Y[second_part:rows,]
}
save(X_1, Y_1, file="client1_latest.RData")
save(X_2, Y_2, file="client2_latest.RData")


# method 2
# cut data into two parts
# K = 10
# m = 4
# part_size = K/2
# Part1_X = rep(list(matrix(0, 1, m)), part_size)
# Part2_X = rep(list(matrix(0, 1, m)), part_size)
# 
# Part1_Y = rep(list(matrix(0, 1, 1)), part_size)
# Part2_Y = rep(list(matrix(0, 1, 1)), part_size)
# for (i in c(1:part_size)) {
#   Part1_X[i] = X[i]
#   Part2_X[i] = X[i+part_size]
#   Part1_Y[i] = Y[i]
#   Part2_Y[i] = Y[i+part_size]
# }
# 
# save(Part1_X, Part1_Y, file="client1.RData")
# save(Part2_X, Part2_Y, file="client2.RData")
# 
# 
# # add zero for part_1 and part2
# tmp = matrix(0,1,1)
# for (i in 6:10) {
#   Part1_X[[i]] = tmp
#   Part1_Y[[i]] = tmp
# }
# for (i in 1:5) {
#   Part2_X[[i+5]] = Part2_X[[i]]
#   Part2_Y[[i+5]] = Part2_Y[[i]]
#   Part2_X[[i]] = tmp
#   Part2_Y[[i]] = tmp 
# }
# 
# save(Part1_X, Part1_Y, file="client1_padding.RData")
# save(Part2_X, Part2_Y, file="client2_padding.RData")


