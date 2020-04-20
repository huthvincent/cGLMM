# Fed_GLMM

Implementing the distributed version of GLMM in R

### Origin Version
muti_fe_latest.Rmd

### Generate simulated data
Generate_test_data.r is used for generating simulated data, two files saved at current fold named as "client1_latest.RData", "client2_latest.RData"

### Distributed Version
evalServer6.R, Fed_GLMM_client_function.r, Fed_GLMM_server.r, runFed_GLMM.r. The first file used for socket communication; two file represent client and server

### Run instruction
1. set the same fold as working direcotry for both client and server source code.
2. lanch at least 3 session stand for server, client1, and client2 in Rstudio
3. at each session, run the corresponding scripts at runFed_GLMM.r

### Part of Experiment report
https://docs.google.com/document/d/1gqOd-wv0CKdjdBndndxDSni4QSv6UKzjUTV-0JnDhAI/edit?usp=sharing
