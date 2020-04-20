setwd("/home/huthvincent/Desktop/cGLMM_latest/Fed_GLMM_Latest")

## creating data
#n个病人一共
#m个feature
#r.e有K种
#生成数据函数
n = 1000
K = 2
m = 10

generate_data <- function(n = 1000,
                          m = 10,
                          beta = c(1:10),
                          K = 2,
                          sigma1 = 0.5) {
  X1 = matrix(rnorm(n * m), n, m)
  r_e = round(runif(n, 1, K))
  
  X = rep(list(matrix(0, 1, m)), K)
  Y = rep(list(matrix(0, 1, 1)), K)
  
  for (i in c(1:n)) {
    X[[r_e[i]]] = rbind(X[[r_e[i]]], unlist(X1[i,]))
    Y[[r_e[i]]] = rbind(Y[[r_e[i]]], c(0))
  }
  
  for (i in c(1:K)) {
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

real_parameters <- function(X, Y, Z1, initial_beta = c(1, 1, 3, 2)) {
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

