
# server
options(digits = 15)
setwd("/home/huthvincent/Desktop/glmm/150/Fed_GLMM")
source("Fed_GLMM_server.r")

#client 1
setwd("/home/huthvincent/Desktop/glmm/150/Fed_GLMM")
source("Fed_GLMM_client_function.r")
Fed_GLMM_client("client1","GLMMtest","",8999,"client1_latest.RData")
#Fed_GLMM_client("client1","GLMMtest","",8999,"client1_padding.RData")
# function(username,taskid,password,portNumber,localFileName)

#client 2
setwd("/home/huthvincent/Desktop/glmm/150/Fed_GLMM")
source("Fed_GLMM_client_function.r")
Fed_GLMM_client("client2","GLMMtest","",8999,"client2_latest.RData")
#Fed_GLMM_client("client2","GLMMtest","",8999,"client2_padding.RData")

#clinet 3   
#setwd("~/Desktop/Fed_GLMM")
#source("GLOREclient_function.r")
#GLOREclient("dbmiclient3","dbmiGLROEtest","",8999,"edin_part3.txt")