#client 1
setwd("/home/huthvincent/Desktop/glmm/Fed_GLMM")
source("Fed_GLMM_client_function.r")
Fed_GLMM_client("client1","GLMMtest","",8999,"client1_latest.RData")
#Fed_GLMM_client("client1","GLMMtest","",8999,"client1_padding.RData")
# function(username,taskid,password,portNumber,localFileName)
