# Fed_GLMM

Implementation of the distributed version of GLMM with Homomorphic encryption in R

## Non-distributed Version of GLMM
muti_fe_latest.Rmd

## Distributed Version
### Generate simulated data
Run Generate_test_data.r to generate two files stands for two clients
### Install requried packages
require(tcltk)
require(svMisc)
require(svSocket)
library(gmp)
library(homomorpheR)
### Run instruction
1. Run initial_parameters.R to save initial parameters as initial_parameters.RData
2. Run generate_keyPair.R to save Homomorphic Key Pair as priv_key.RData, pub_key.RData
3. Lanch Rstudio and set the same fold as working direcotry for both client (Fed_GLMM_client_function.r) and server (Fed_GLMM_server.r).
4. Lanch at least 3 session stand for server, client1, and client2 in Rstudio
5. At each session, run the corresponding scripts at runFed_GLMM.r

## Ref:
[Introduction to Homomorphic Computation in R](
http://cran.r-project.org.icopy.site/web/packages/homomorpheR/vignettes/introduction.html)
