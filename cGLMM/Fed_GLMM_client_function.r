require(tcltk)
require(svMisc)
require(svSocket)
source("evalServer6.R")


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
    cat("print file name")
    cat(localFileName)
    data = load(localFileName)
    X = get(data[1])
    Y = get(data[2])
    
    # Random initial parameters like beta and sigma1 etc.
    #beta <- c(1,2,4,2)
    #K <- 100
    #burnin <- 30
    #sigma1 <- 3
    #n <- 10
    
    # synchronize with server for iteration
    iter <- 0
    step <- 1
    
    sleep_time = 0.1
    epsilon <- 1e10
    threshold <- 0.1
    
    # get parameters from server
    num_fe = evalServer6(con, "num_fe=num_fe")
    K = evalServer6(con, "K=K")
    burnin = evalServer6(con, "burnin=burnin")
    n = evalServer6(con, "n=n")
    
    # whether the beta and sigma are converged
    while (epsilon > threshold) {
      epi <- 1
      while (epi > 0.1) {
        # get beta and sigma from server
        sigma1 = evalServer6(con, "sigma1=sigma1")
        beta = evalServer6(con, "beta=beta")
        Z1 = evalServer6(con, "Z1=Z1")
        
        for (i in 1:(K + burnin)) {
          # computing A1 on local
          #cat("client i = ", i, "\n")
          if (i%%200 == 0){
           cat("client: i = ", i, "\n")
          }
          # processing client one by one
          repeat {
            Sys.sleep(sleep_time)
            if (evalServer6(con, "i=i") == i) {
              #cat("start i = ", i, "\n")
              # waiting for server processing
              while (evalServer6(con, "receive_ind=receive_ind") < clientID) {
                Sys.sleep(sleep_time)
              }
              newZ1 = evalServer6(con, "newZ1=newZ1")
              Z1 = evalServer6(con, "Z1=Z1")
              A1 = rep(0, n)
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
              # sending A1 to the server
              evalServer6(con, A1[, , receive_ind], A1)
              
              ready = evalServer6(con, "ready=ready")
              while (ready == 0) {
                Sys.sleep(sleep_time)
                ready = evalServer6(con, "ready=ready")
              }
              #cat("server finished i=", i, "\n")
              break
            }
          }
        }
        
        Z1_list_flag = evalServer6(con, "Z1_list_flag=Z1_list_flag")
        while (Z1_list_flag == 0) {
          Sys.sleep(sleep_time)
          Z1_list_flag = evalServer6(con, "Z1_list_flag=Z1_list_flag")
        }
        #cat("Z1_list_flag = ", Z1_list_flag, "\n")
        Z1_list = evalServer6(con, "Z1_list=Z1_list")
        #cat("Z1_list is", Z1_list,"\n")
        # computing Hessian matrix locally
        H = matrix(0, num_fe, num_fe)
        for (i in c(1:num_fe)) {
          for (j in c(i:num_fe)) {
            f2 = rep(0, K)
            for (mc in c(1:K)) {
              temp = rep(0, n)
              for (k in c(1:n)) {
                num_pat = dim(X[[k]])[1]
                if (num_pat > 1) {
                  temp[k] = sum(X[[k]][, i] * X[[k]][, j] * exp(X[[k]] %*% beta  + rep(Z1_list[k, mc], num_pat)) /
                                  (1 + exp(
                                    X[[k]] %*% beta  + rep(Z1_list[k, mc], num_pat)
                                  )) ^ 2)
                }
              }
              f2[mc] = sum(temp)
            }
            H[i, j] = sum(f2)
            H[j, i] = H[i, j]
            #cat("H[",i,j,"]=",sum(f2))
          }
        }
        #cat("finish computing Hessian matrix", H)
        # computing first derivative locally
        f1 = rep(0, num_fe)
        for (i in c(1:num_fe)) {
          temp2 = rep(0, K)
          for (mc in c(1:K)) {
            temp = rep(0, n)
            for (k in c(1:n)) {
              num_pat = dim(X[[k]])[1]
              if (num_pat > 1) {
                temp[k] = sum(X[[k]][, i] * (Y[[k]] - 1 + 1 / (1 + exp(
                  X[[k]] %*% beta  + rep(Z1_list[k, mc], num_pat)
                ))))
              }
            }
            temp2[mc] = sum(temp)
          }
          f1[i] = sum(temp2)
         # cat("f1[",i,"]=",sum(temp2))
        }
        #cat("finish computing f1 ", f1)
        while (evalServer6(con, "receive_ind=receive_ind") < clientID) {
          Sys.sleep(sleep_time)
        }
        evalServer6(con, H[, , receive_ind], H)
        evalServer6(con, f1[, , receive_ind], f1)
        #cat("sending H and f1 to server")
        # Get the epi
        newiter = evalServer6(con, "iter=iter")
        while (iter == newiter) {
          newiter = evalServer6(con, "iter=iter")
        }
        epi = evalServer6(con, "epi=epi")
        iter = newiter
        cat("Iteration=", iter, "\n")
        flush.console()
      }
      
      # Get the step
      newstep = evalServer6(con, "step=step")
      while (step == newstep) {
        newstep = evalServer6(con, "step=step")
      }
      step = newstep
      cat("Step=", step, "\n")
      flush.console()
    }
    
    while (evalServer6(con, "receive_ind=receive_ind") < clientID)
    {
      Sys.sleep(sleep_time)
    }
    
    evalServer6(con, stop_connection[receive_ind], 1)
    close(con)
    res <- list(beta, sigma)
    names(res) <- c("coefficients", "sigma")
    cat("result is", res)
    flush.console()
  }

#username <- "dbmiclient3"
#	taskid <- "dbmiGLROEtest"
#	password <- ""
#	portNumber <- 8971
#	localFileName <- "C:/Users/yuan/Dropbox/distributed_lr_paper/edin_part3.txt"
