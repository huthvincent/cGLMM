require(tcltk)
require(svMisc)
require(svSocket)
require(Metrics)
options(digits = 15)
MH <- function(A1, Z1, newZ1, K, burnin, n) {
 
    # update Z1
    Z1 <- ifelse(log(runif(n)) > A1, Z1, newZ1)
    #Z1_MH[, i] <- Z1
    #cat("Z1_MH[,",i,"]=", Z1)
    sum_t = 0
    sum_t = sum(ifelse(Z1 == newZ1, yes = 1, 0))
    # print("z1")
    # print(Z1)
    # print("newz1")
    # print(newZ1)
    # print("sum_t")
    # print(sum_t)
    return(list(Z1, sum_t))
}

# authenticat two clients with task ID and password
u_username = ""
usernames <- c("client1", "client2")
n_user <- length(usernames)
taskid <- "GLMMtest"
passwords <- c("", "")
portNumber <- 8999

load("/home/huthvincent/Desktop/glmm/150/Fed_GLMM/ground_truth_beta.RData")

ground_truth = true_beta


epsilon <- 1e10
threshold <- 0.01
sleep_time = 0.1
# parameters or computing coefficient
num_fe <- 150
load("/home/huthvincent/Desktop/glmm/150/Fed_GLMM/init_beta.Rdata")
beta = init_beta
sigma1 <- 5.6

beta_list <- matrix(beta, num_fe, 1)
sigma1_list <- c(sigma1)
mse_list = c()
epi_list = c()


step <- 1
K <- 2
burnin <- 998
#K <- 1000
# burnin <- 300
n = 10 # sampling size for Metropolis Hasting

# for receiving Hessian matrix from clients
H <- array(10^(-20), c(num_fe, num_fe, n_user))
A1 <- array(10^(-20), c(1, n, n_user))
f1 <- array(10^(-20), c(1, num_fe, n_user))

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
timestart<-Sys.time()
# epsilon controls whether the beta and sigma are converged
while (epsilon > threshold) {
  
  # epi controls whether the beta is converged
  epi <- 1
#  while (epi > 0.1) {
    
    # Step 1: Sampling Z1_list
    Z1_MH = array(0, dim = c(n, K + burnin))
    Z1 = rnorm(n, 0, sigma1)
    Z1_list_flag <- 0
    sum_sum_t = 0
    for (i in 1:(K + burnin)){
      #cat("server i = ", i, "\n")
      if (i%%200 == 0){
          cat("server: i = ", i, "\n")
      }
      newZ1 = rnorm(n, mean = 0, sd = sigma1)
      # each Z1_MH sampling need gather A1 from all clients
      for (receive_ind in 1:n_user) {
        #cat("Waiting for client ID is", receive_ind)
        while (mean(abs(A1[,,receive_ind]))==10^(-20)) {
          #cat("waiting")
          Sys.sleep(sleep_time)
        }
        #cat("A1[,,",receive_ind,"] is ", A1[,,receive_ind], "\n")
      }
      ready <- 0
      # get A1 from clients, added elementwise and pass to MH function
      A_sum = rowSums(A1,dims = 2)
      # get Z1_MH[,i]
      Z1_sum_t = MH(A_sum, Z1, newZ1, K, burnin, n)
      # print("z1")
      # print(Z1_sum_t[[1]])
      # print("sum_t")
      # print(Z1_sum_t[[2]])
      sum_t = Z1_sum_t[[2]]
      sum_sum_t = sum_sum_t + sum_t
      Z1 = Z1_sum_t[[1]]
      Z1_MH[,i] = Z1
      A1 <- array(10^ (-20), c(1, n, n_user))
      ready <- 1
      #cat("server finish i=", i, "\n")
      receive_ind <- 0
    }
    cat("SUM_SUM_t is", "\n")
    print(sum_sum_t)
    Z1_list = Z1_MH[,-c(1:burnin)]
    #cat("Z1_list in Server is ", Z1_list, "\n")
    Z1_list_flag <- 1
    # update sigma1
    sigma1 = sqrt(mean(Z1_list ^ 2))
    
    receive_ind <- 0
    # get from clients Hessian Matrix and first derivative
    for (receive_ind in 1:n_user) {
      #cat("get from client ", receive_ind, "\n")
      #cat("H is", H[,,receive_ind])
      #cat("f1 is", f1[,,receive_ind])
      while (mean(abs(H[,,receive_ind]))==10^(-20) || mean(abs(f1[,,receive_ind]))==10^(-20)) {
        #cat("waiting")
        Sys.sleep(sleep_time)
      }
      #cat("Received H is", H[,,receive_ind])
      #cat("Received f1 is ", f1[,,receive_ind])
    }
    
    H_sum = rowSums(H,dims = 2)
    #cat("H_Sum is", H_sum)
    H_sum = matrix(H_sum, num_fe, num_fe)
    H_inverse = solve(H_sum)
    #cat("H_inverse", H_inverse)
    if (det(H_inverse) == 0) {
      cat("H_inverse is 0", "\n")
      break
    }
   
    f1 = rowSums(f1, dims = 2)
    f1 = matrix(f1, nrow=num_fe, ncol=1)
    # update beta
    updata = H_inverse %*% f1
    new_beta = beta + updata
    epi = mean(abs(new_beta - beta))
    beta <- new_beta
    
    # reset H and f1 for receiving the next iteration 
    H <- array(10^(-20), c(num_fe, num_fe, n_user))
    f1 <- array(10^(-20), c(1, num_fe, n_user))
    iter <- iter + 1
    cat("Iteration=", iter, "\n")
    #cat("beta=", beta,"\n")
    cat("MSE = ", "\n")
    print(mean(abs((ground_truth- beta)/ground_truth)))
    mse_list = append(mse_list,mean(abs((ground_truth- beta)/ground_truth)))  
    
    cat("sigma=", sigma1, "\n")
    flush.console()
#  }
  
  epsilon <-
    max(mean(abs(beta - beta_list[,step])), abs(sigma1 - sigma1_list[step])/20)
  print("epsilon")
  print(epsilon)
  epi_list = append(epi_list,epsilon)
  beta_list = cbind(beta_list, beta)
  sigma1_list[step + 1] <- sigma1
  step <- step + 1
  cat("Step is", step, "\n")
  flush.console()
}

# print the final result
#print("the final beta is")
#print(beta)
print("the final sigma is")
print(sigma1)

# break the connections
receive_ind <- 0
for (receive_ind in 1:n_user)
{
  while (stop_connection[receive_ind] != 1)
  {
    Sys.sleep(sleep_time)
  }
}

Sys.sleep(sleep_time)
stopSocketServer(port = portNumber)
# endof Fed_GLMM_server
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)  

save(epi_list,file='/home/huthvincent/Desktop/cGLMM/draw_graph/cGLMM_epi.Rdata')
save(mse_list,file='/home/huthvincent/Desktop/glmm/150/Fed_GLMM/dataANA/cGLMM_mse.Rdata')
save(sigma1_list,file='/home/huthvincent/Desktop/glmm/150/Fed_GLMM/dataANA/cGLMM_sigma.Rdata')
