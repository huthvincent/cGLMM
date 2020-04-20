require(tcltk)
require(svMisc)
require(svSocket)
library(gmp)
library(homomorpheR)
source("evalServer6.R")
load("/home/huthvincent/Desktop/cGLMM_latest/Fed_GLMM_Latest/pub_key.RData")
load("/home/huthvincent/Desktop/cGLMM_latest/Fed_GLMM_Latest/priv_key.RData")
# load initial parameters from file
load("/home/huthvincent/Desktop/cGLMM_latest/Fed_GLMM_Latest/initial_parameters.RData")

# receive encrypted data from server
server2Client <- function(BigZ, ser_num) {
  decrypted = priv_key$decrypt(BigZ)
  decrypted = decrypted / Precision
  decrypted = decrypted - ser_num * Add_Big_Number
  return (decrypted)
}
# sending encrypted data to server
client2Server <- function(number){
  number = number + Add_Big_Number
  number = number * Precision
  return (pub_key$encrypt(number))
}

Fed_GLMM_client <-
  function(username,
           taskid,
           password,
           portNumber,
           localFileName)
  {
    con = socketConnection(host = "127.0.0.1",
                           port = portNumber,
                           blocking = FALSE)
    # local attestation with server
    authenticationFlag = FALSE
    server_taskid <- evalServer6(con, "temp=taskid")
    cat(server_taskid)
    evalServer6(con, u_username, username)
    if (server_taskid == taskid)
    {
      clientID = evalServer6(con, "clientID=which(usernames==u_username)")
      cat("clientID=", clientID, "\n")
      if (length(clientID))
        authenticationFlag = TRUE
    }
    
    if (authenticationFlag == FALSE)
    {
      close(con)
      stop("Please check your taskid and username and rerun your local code")
    }
    
    # load generated data for testing only
    cat("loading data from file")
    cat(localFileName)
    data = load(localFileName)
    X = get(data[1])
    Y = get(data[2])
    
    # synchronize with server for iteration
    iter <- 0
    step <- 1
    
    epsilon <- 1e10
    threshold <- 0.05
    
    H_sum_receiver = as.bigz(matrix(rep(0, num_fe*num_fe), num_fe, num_fe))
    sum_f1_receiver = as.bigz(matrix(rep(0,num_fe), num_fe))
    
    # whether the beta and sigma are converged
    while (epsilon > threshold) {
      epi <- 1
      while (epi > 0.1) {
        
        # get beta and sigma from server
        sigma1 = evalServer6(con, "sigma1=sigma1")
        Z1 = evalServer6(con, "Z1=Z1")
        
        for (i in 1:(K + burnin)) {
          # computing A1 on local
          cat("\nclient: i = ", i, "\n")
          #if (i%%200 == 0){
          # cat("client: i = ", i, "\n")
          #}
          # processing client one by one
          repeat {
            Sys.sleep(0.1)
            if (evalServer6(con, "i=i") == i) {
              # waiting for server processing
              while (evalServer6(con, "receive_ind=receive_ind") < clientID) {
                Sys.sleep(0.1)
              }
              newZ1 = evalServer6(con, "newZ1=newZ1")
              A1 = as.bigz(rep(0, n))
              for (j in c(1:n)) {
                num_pat = dim(X[[j]])[1]
                if (num_pat > 1) {
                  A1[j] <-
                    sum(Y[[j]] * (newZ1[j] - Z1[j])) + sum(log(1 + exp(
                      X[[j]] %*% beta + rep(Z1[j], num_pat)
                    )) - log(1 + exp(
                      X[[j]] %*% beta + rep(newZ1[j], num_pat)
                    )))
                }
              }
              # sending encrypted A1 to the server
              print("sending A1 to server")
              encrypted_A1 = client2Server(A1)
              evalServer6(con, A1[[receive_ind]], encrypted_A1)
              
              print("wating server get All A1")
              ready = evalServer6(con, "ready=ready")
              while (ready == 0) {
                Sys.sleep(0.1)
                #evalServer6(con, A1[[receive_ind]], encrypted_A1)
                ready = evalServer6(con, "ready=ready")
              }
              #cat("server finished i=", i, "\n")
              break
            }
          }
        }
        
        Z1_list_flag = evalServer6(con, "Z1_list_flag=Z1_list_flag")
        while (Z1_list_flag == 0) {
          Sys.sleep(0.1)
          Z1_list_flag = evalServer6(con, "Z1_list_flag=Z1_list_flag")
        }
        Z1_list = evalServer6(con, "Z1_list=Z1_list")
        
        print("-------Z1_list received from server -------")
        print(Z1_list)
        
        # computing Hessian matrix locally
        
        print("computing Hessian")
        H = matrix(0, num_fe, num_fe)
        for (H_i in c(1:num_fe)) {
          for (H_j in c(1:num_fe)) {
            f2 = rep(0, K)
            for (mc in c(1:K)) {
              temp = rep(0, n)
              for (k in c(1:n)) {
                num_pat = dim(X[[k]])[1]
                if (num_pat > 1) {
                  temp[k] = sum(X[[k]][, H_i] * X[[k]][, H_j] * exp(X[[k]] %*% beta  + rep(Z1_list[k, mc], num_pat)) /
                                  (1 + exp(
                                    X[[k]] %*% beta  + rep(Z1_list[k, mc], num_pat)
                                  )) ^ 2)
                }
              }
              f2[mc] = sum(temp)
            }
            H[H_i, H_j] = sum(f2)
          }
        }
        print("computing f1")
        # computing first derivative locally
        f1 = rep(0, num_fe)
        for (f_i in c(1:num_fe)) {
          temp2 = rep(0, K)
          for (mc in c(1:K)) {
            temp = rep(0, n)
            for (k in c(1:n)) {
              num_pat = dim(X[[k]])[1]
              if (num_pat > 1) {
                temp[k] = sum(X[[k]][, f_i] * (Y[[k]] - 1 + 1 / (1 + exp(
                  X[[k]] %*% beta  + rep(Z1_list[k, mc], num_pat)
                ))))
              }
            }
            temp2[mc] = sum(temp)
          }
          f1[f_i] = sum(temp2)
        }
        
        # synchronize the order of client
        while (evalServer6(con, "receive_ind=receive_ind") < clientID) {
          Sys.sleep(0.1)
        }
        
        # encrypting H and f1
        encrypted_H = client2Server(H)
        encrypted_f1 = client2Server(f1)
        print("send out H and f")
        evalServer6(con, encrypted_H[[receive_ind]], encrypted_H)
        evalServer6(con, encrypted_f1[[receive_ind]],encrypted_f1)
        
        # wait until H_sum is ready
        H_sum_flag = evalServer6(con, "H_sum_flag=H_sum_flag")
        while (H_sum_flag == 0) {
          Sys.sleep(0.1)
          H_sum_flag = evalServer6(con, "H_sum_flag=H_sum_flag")
        }
        cat("Get H_sum from server ", H_sum_flag, "\n")
        # receive H_sum and decrypte and encrypt send out
        
        H_sum_receiver = evalServer6(con, "H_sum=H_sum")
        
        # wait until sum_f1_flag is ready
        sum_f1_flag = evalServer6(con, "sum_f1_flag=sum_f1_flag")
        while (sum_f1_flag == 0) {
          Sys.sleep(0.1)
          sum_f1_flag = evalServer6(con, "sum_f1_flag=sum_f1_flag")
        }
        
        sum_f1_receiver = evalServer6(con, "sum_f1=sum_f1")
        
        print("decrypted H_sum")
        print(as.numeric(server2Client(H_sum_receiver, n_user+1)))
        H_inverse = solve(matrix(as.numeric(server2Client(H_sum_receiver, n_user+1)),num_fe, num_fe))
        if (det(H_inverse) == 0) {
          print("det == 0")
          break
        }
        
        print("decrypted sum_f1")
        print(as.numeric(server2Client(sum_f1_receiver, n_user+1)))
        sum_f1 = as.numeric(server2Client(sum_f1_receiver, n_user+1))
        
        # update beta
        print("H_inverse")
        print(H_inverse)
       
        sum_f1 = matrix(sum_f1, nrow=num_fe, ncol=1)
        updata = H_inverse %*% sum_f1
        
        beta = beta + updata
        
        while (evalServer6(con, "receive_ind=receive_ind") < clientID) {
          Sys.sleep(0.1)
        }
        
        # send out beta
        evalServer6(con, beta_receiver[[receive_ind]],beta)
        cat("old iter =", iter, "\n")
        # Get the epi
        newiter = evalServer6(con, "iter=iter")
        cat("newiter = ", newiter, "\n")
        while (iter == newiter) {
          Sys.sleep(0.1)
          newiter = evalServer6(con, "iter=iter")
        }
        epi = evalServer6(con, "epi=epi")
        iter = newiter
        cat("Iteration=", iter, "finish \n")
        print("beta")
        print(beta)
      }
      
      # Get the step
      newstep = evalServer6(con, "step=step")
      while (step == newstep) {
        newstep = evalServer6(con, "step=step")
      }
      step = newstep
      cat("Step=", step, "\n")
    }
    
    while (evalServer6(con, "receive_ind=receive_ind") < clientID)
    {
      Sys.sleep(0.1)
    }
    
    evalServer6(con, stop_connection[receive_ind], 1)
    close(con)
    res <- list(beta, sigma)
    names(res) <- c("coefficients", "sigma")
    cat("result is", res)
    
  }

#username <- "dbmiclient3"
#	taskid <- "dbmiGLROEtest"
#	password <- ""
#	portNumber <- 8971
#	localFileName <- "C:/Users/yuan/Dropbox/distributed_lr_paper/edin_part3.txt"
