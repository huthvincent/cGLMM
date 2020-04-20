setwd("~/Desktop/Fed_GLMM 1024")
require(Rlab)
## creating data
#n个病人一共
#m个feature
#r.e有K种
#生成数据函数
generate_data <- function(n = 300,
                          m = 100,
                          beta = c(-49:50)/100,
                          K = 10,
                          sigma1 = 3,
                          p = 0.5) {
  
  pre_X1 = rep(0,n*m)
  for(v in c(1:(n*m))){
    genome_seq = rbern(2, p)
    if(genome_seq[1] == 0){
      if(genome_seq[2] == 0){
        pre_X1[v] = 0
      }
    }
    if(genome_seq[1] == 1){
      if(genome_seq[2] == 1){
        pre_X1[v] = 2
      }
    }
    else{
      pre_X1[v] = 1
    }
  }
  X1 = matrix(pre_X1,n,m)
  
  r_e = round(runif(n,1,K))
  
  X = rep(list(matrix(0,1,m)),K)
  Y = rep(list(matrix(0,1,1)),K)
  
  for(i in c(1:n)){
    
    X[[r_e[i]]] = rbind(X[[r_e[i]]],unlist(X1[i,]))
    Y[[r_e[i]]] = rbind(Y[[r_e[i]]],c(0))
  }    
  
  for(i in c(1:K)){
    X[[i]] = X[[i]][-1,] 
    Y[[i]] = Y[[i]][-1,]
  }
  
  #制造一个与X同结构的 y放label
  
  Z1 = rnorm(K, 0, sigma1)
  
  #matrix(U,n,T_)是一个每一列是U，一共T_列的矩阵
  #x是10列100行的矩阵 每个元素加上对应的Z1和Z2
  for (i in c(1:K)) {
    patient_amount = dim(X[[i]])[1]
    x = X[[i]] %*% beta + Z1[i]
    P <- exp(x) / (1 + exp(x))
    Y[[i]] <-
      ifelse(matrix(runif(patient_amount * 1), patient_amount, 1) < P, 1, 0)
  }
  
  return(list(X, Y, Z1, r_e))
}
alldata <- generate_data()
X = alldata[[1]]
Y = alldata[[2]]
Z = alldata[[3]]
Z_index = alldata[[4]]
Z_index

# calcuate the real_parameters by taking the Z1 as know in advance

real_parameters <- function(X, Y, Z1, initial_beta = c(-49:50)/100) {
  num_fe = dim(X[[1]])[2]
  
  K = length(Z1)
  #real_sigma1 就是 Z1的
  real_sigma1 <- sqrt(mean(Z1 ^ 2))
  real_beta = initial_beta
  epi <- 2
  
  #这里的f2就是likelihood l 对 beta的二阶导
  while (epi > 1e-2) {
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

save.image(file="cur_session.RData")


# method 1
# cut data into two parts both has 10 random effects 
K = 10
m = 100
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


