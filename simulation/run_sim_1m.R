# The first two lines MUST be the comments.This Template if for Type I error simulation
#This is template file used by RunmefileCreator.R to create a number of R files, and a command bash file to run the sime. Multiple R files are created to split the simulation run to be run on each node/machine
#rm(list = ls())
run <-  1 ; numFile <-  1
nSim <-  5000
rep <- paste("rep");assign(rep, seq( (1+(1-1)*ceiling(5000/numFile)), ((1-1)*ceiling(5000/numFile)+5000 /numFile)) )
library(lavaan)
library(RMediation)
library(OpenMx)
library(doParallel)
library(dplyr)
source("mbcoSimFuns.R", echo = FALSE)     # modified 4/7/18
source("profLikSimFuns.R", echo = FALSE)  # modified 4/7/18
source("simFunctions_within_lav.R", echo = FALSE) # main simulation function, # modified 4/28/18

Quant <- c(~ b1 * b2 *b3)
type <- c("mbco", "joint", "prof", "mc", "perc", "bc") ###Methods
b1 <- b2 <- b3 <- c(0, 0.14, 0.39, 0.59)
samp <- c(50, 75, 100, 200, 500) #sample Sizes to generate mediation model
R = 1000L # The number of Bootstrap samples

dfParm <-
  expand.grid(rep=rep,
              quant = Quant,
              samp = samp,
              b1 = b1,
              b2 = b2,
              b3 = b3,
              s2x=1,
              s2em1=1,
              s2em2=1,
              s2ey=1)


dfParm1 <- dfParm %>% dplyr::filter(b1==0, b2==0, b3==0)
dfParm2 <- dfParm %>% dplyr::filter(b2!=0, b1>=b2 ,b3==0)
dfParm3 <- dfParm %>% dplyr::filter(b1!=0,b2==0,b3==0)
dfParm4 <- dfParm %>% dplyr::filter(b3!=0, b1>=b2, b2>=b3, b1>=b3) %>% arrange(b1,b2,b3)


dfParmn <- bind_rows(dfParm1,dfParm2,dfParm3,dfParm4)%>% arrange(b1,b2,b3)

mc <- detectCores() - 1
cl <- makeCluster(mc)
registerDoParallel(cl)

it <- iter(dfParm, by = 'row')
resSim <-
  foreach(
    x = it,
    .combine = "rbind",
    .packages = c("RMediation", "lavaan", "OpenMx"),
    .inorder = FALSE,
    .errorhandling = "remove",
    .verbose = F
  ) %dopar% sim(x, type = type)

stopCluster(cl)

write.table(resSim,
            file = file.path("res2M",paste("run", run, ".csv", sep = "")),
            row.names = FALSE,
            sep = ",")
