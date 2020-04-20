require(tcltk)
require(svMisc)
require(svSocket)
library(gmp)
library(homomorpheR)
load("/home/huthvincent/Desktop/cGLMM_latest/Fed_GLMM_Latest/pub_key.RData")
# parameters or computing coefficient loading from file
load("/home/huthvincent/Desktop/cGLMM_latest/Fed_GLMM_Latest/initial_parameters.RData")

# debug 
load("/home/huthvincent/Desktop/cGLMM_latest/Fed_GLMM_Latest/priv_key.RData")

BigzMean <- function(BigZ) {
  sumZ <- as.bigz(0)
  len <- length(BigZ)
  for (idx in c(1:len)) {
    sumZ <- add.bigz(sumZ, BigZ[idx])
  }
  mean_val <- sumZ/length(BigZ)
  return (mean_val)
}

BigzMatrixMean <- function(BigZ) {
  sumZ = as.bigz(0)
  len = dim(BigZ)[1]
  for (i in c(1:len)){
    for (j in c(1:len)) {
      sumZ <- add.bigz(sumZ, BigZ[i,j])
    }
  }
  mean_val <- sumZ/(len*len)
  return (mean_val)
}

server2Client <- function(BigZ, n) {
  decrypted = priv_key$decrypt(BigZ)
  decrypted = decrypted / Precision
  decrypted = decrypted - n * Add_Big_Number
  return (decrypted)
}

client2Server <- function(number){
  number = number + Add_Big_Number
  number = number * Precision
  return (pub_key$encrypt(number))
}

MH <- function(A1, Z1, newZ1, K, burnin, n) {
    # A1 is encrypted, Z1 and newZ1 is not enctyped
    # update Z1
    Z1 <- ifelse(client2Server(log(runif(n)) + (n_user+1) * Add_Big_Number) > A1, Z1, newZ1)
    return(Z1)
}


# authenticat two clients with task ID and password
u_username = ""
usernames <- c("client1", "client2")
n_user <- length(usernames)
taskid <- "GLMMtest"
passwords <- c("", "")
portNumber <- 8999

beta_list <- matrix(beta, num_fe, 1)
sigma1_list <- c(sigma1)
step <- 1

# for receiving Hessian matrix f1 and beta from clients
encrypted_H <- list(as.bigz(matrix(rep(0,num_fe*num_fe), num_fe, num_fe)), as.bigz(matrix(rep(0, num_fe*num_fe), num_fe, num_fe)))
A1 <- list(as.bigz(rep(0, n)), as.bigz(rep(0, n)))
encrypted_f1 <- list(as.bigz(matrix(rep(0,num_fe), num_fe)), as.bigz(matrix(rep(0,num_fe), num_fe)))
beta_receiver <- list(rep(0, num_fe), rep(0, num_fe))

# control the connections and the ordre of clients communication
receive_ind <- 0
stop_connection <- rep(0, n_user)
startSocketServer(
  port = portNumber,
  server.name = "R_yw",
  procfun = processSocket,
  secure = FALSE,
  local = FALSE
)

# algorithm iterations
iter <- 0

# epsilon controls whether the beta and sigma are converged
while (epsilon > threshold) {
  
  # epi controls whether the beta is converged
  epi <- 1
  while (epi > 0.1) {
    
    # Step 1: Sampling Z1_list
    Z1_MH = array(0, dim = c(n, K + burnin))
    Z1 = rnorm(n, 0, sigma1)
    Z1_list_flag <- 0
    H_sum_flag <- 0
    
    for (i in 1:(K + burnin)){
      cat("server: i = ", i, "\n")
      
      #if (i%%200 == 0){
      #    cat("server: i = ", i, "\n")
      #}
      newZ1 = rnorm(n, mean = 0, sd = sigma1)
      
      # each Z1_MH sampling need gather A1 from all clients
      for (receive_ind in 1:n_user) {
        cat("get A1 from client ", receive_ind, "\n")
        while (BigzMean(abs(A1[[receive_ind]])) == 0) {
          Sys.sleep(0.1)
        }
        
      }
      ready <- 0
      # get encrypted A1 from clients, added elementwise and pass to MH function
      A_sum = client2Server(as.bigz(rep(0,n)))
      for (p in c(1:n)) {
        for (q in c(1:n_user)){
          A_sum[p] = A_sum[p] * A1[[q]][p]
        }
      }
      # get Z1_MH[,i]
      Z1_MH[,i] = MH(A_sum, Z1, newZ1, K, burnin, n)
      Z1 = Z1_MH[,i]
      
      #reset A1
      A1 <- list(as.bigz(rep(0,n)), as.bigz(rep(0,n)))
      ready <- 1
      receive_ind <- 0
    }
    # Z1_list and sigma1 not encrypted
    print("\nZ1_MH\n")
    print(Z1_MH)
    Z1_list = Z1_MH[,-c(1:burnin)]
    print("Server : Z1_list ")
    print(Z1_list)
    
    Z1_list_flag <- 1
    # update sigma1
    sigma1 = sqrt(mean(Z1_list ^ 2))
    
    receive_ind <- 0
    # get from clients Hessian Matrix and first derivative
    for (receive_ind in 1:n_user) {
      cat("get H and f1 from client ", receive_ind, "\n")
      
      while (BigzMatrixMean(encrypted_H[[receive_ind]]) == 0 || BigzMean(encrypted_f1[[receive_ind]]) == 0) {
        #cat("waiting")
        Sys.sleep(0.1)
      }
      cat("decrypt Received H is",as.numeric(server2Client(encrypted_H[[receive_ind]], 1)) )
    }
    
    # H_sum is not ready
    sum_f1_flag = 0
    H_sum_flag = 0
    H_sum = client2Server(as.bigz(matrix(rep(0, num_fe*num_fe), num_fe, num_fe)))
    
    for (H_i in c(1:num_fe)) {
      for (H_j in c(1:num_fe)){
        for (H_k in c(1:n_user)){
          # here * == + for homomorphic operation
          H_sum[H_i,H_j] = H_sum[H_i,H_j] * encrypted_H[[H_k]][H_i,H_j]
        }
      }
    }
    print("before send out, check H_sum")
    print(as.numeric(server2Client(H_sum, n_user+1)))
    
    # H_sum is ready
    H_sum_flag = 1
    
    sum_f1 = client2Server(as.bigz(matrix(rep(0,num_fe), num_fe)))
    for (i in c(1:n_user)){
      # here * == + 
      sum_f1 = sum_f1 * encrypted_f1[[i]]
    }
    # sum_f1 ready
    sum_f1_flag = 1
    
    # Waiting client pull H_sum and sum_f1
  
    # get new_beta from client
    for (receive_ind in 1:n_user) {
      cat("get new_beta from client ", receive_ind, "\n")
      while (mean(abs(beta_receiver[[receive_ind]])) == 0) {
        Sys.sleep(0.1)
      }
      
    }
    print("server computing epi")
    new_beta = beta_receiver[[1]]
    epi = max(abs(new_beta - beta))
    beta <- new_beta
    
    print("server reset receiver variables")
    # reset encrypted H and f1 for receiving the next iteration 
    encrypted_H <- list(as.bigz(matrix(rep(0,num_fe*num_fe), num_fe, num_fe)), as.bigz(matrix(rep(0, num_fe*num_fe), num_fe, num_fe)))
    encrypted_f1 <- list(as.bigz(matrix(rep(0,num_fe), num_fe)), as.bigz(matrix(rep(0,num_fe), num_fe)))
    # reset beta_receiver
    beta_receiver <- list(rep(0, num_fe), rep(0, num_fe))
    
    cat("Iteration=", iter, "finish \n")
    iter <- iter + 1
    cat("beta=")
    print(beta)
    cat("sigma=", sigma1, "\n")
  }
  
  epsilon <-
    max(max(abs(beta - beta_list[,step])), abs(sigma1 - sigma1_list[step]))
  
  beta_list = cbind(beta_list, beta)
  sigma1_list[step + 1] <- sigma1
  step <- step + 1
  cat("Step is", step, "\n")
}

# print the final result
print("the final beta is")
print(beta)
print("the final sigma is")
print(sigma1)

# break the connections
receive_ind <- 0
for (receive_ind in 1:n_user)
{
  while (stop_connection[receive_ind] != 1)
  {
    Sys.sleep(0.1)
  }
}

Sys.sleep(1)
stopSocketServer(port = portNumber)
# end of Fed_GLMM_server
