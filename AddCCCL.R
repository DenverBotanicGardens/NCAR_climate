library(ncdf4)
library(extRemes)
library(fields)
library(parallel)
library(doParallel)
library(foreach)
library(abind)

#point to P:/hackathon/Simulations/AddCCCL.r
# scp AddCCCL.r deprengm@yellowstone.ucar.edu:/glade/u/home/deprengm
#In logged in cygwin
# R --no-save < AddCCCL.r



sims <- sprintf('%0.3d', 1:30)

precc2006_2080 <- lapply(sims, function(sims){
  nc_open(paste("/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/PRECC/b.e11.BRCP85C5CNBDRD.f09_g16.",
                sims,".cam.h0.PRECC.200601-208012.nc", collapse="", sep=""))
})

precl2006_2080 <- lapply(sims, function(sims){
  nc_open(paste("/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/PRECL/b.e11.BRCP85C5CNBDRD.f09_g16.",
                sims,".cam.h0.PRECL.200601-208012.nc", collapse="", sep=""))
})

save(precc2006_2080, precl2006_2080, file = "precip.RData")


abind3 <- function(...) { abind(along = 3, ...) }

time <- ncvar_get(precl2006_2080[[1]], "date")

cl <- makeCluster(6)
registerDoParallel(cl)

#use foreach for loop instead of for or lapply
#by default results are returned in a list
#objects from parent environment not available
cc_cl_30 <- lapply(1:30, function(model){
  Ycl <- ncvar_get(precl2006_2080[[model]], "PRECL")
  Ycc <- ncvar_get(precc2006_2080[[model]], "PRECC")
  cccl <- foreach(timeindex = 1:dim(time), 
                  .combine = "abind3",
                  .packages = "ncdf4") %dopar% {
                    cc <- Ycc[,,timeindex]
                    cl <- Ycl[,,timeindex]
                    cc+cl
                  }
  rm(Ycl,Ycc)
  gc()
})

stopCluster(cl)


####Conductive
precc_annual <- lapply(1:30, function(model){
  precc1 <- ncvar_get(precc2006_2080[[model]], "PRECC")
  cl <- makeCluster(6)
  registerDoParallel(cl)
  precc <- foreach(timeindex = seq(1,900,12),
                   .combine = "abind3") %dopar% 
                   {
                     apply(precc1[,,timeindex:(timeindex+11)], MARGIN = c(1,2), sum)
                   }
  stopCluster(cl)
  gc()
  rm(precc1)
  precc
})


####Large Scale precipitation
precl_annual <- lapply(1:30, function(model){
  precl1 <- ncvar_get(precl2006_2080[[model]], "PRECL")
  cl <- makeCluster(6)
  registerDoParallel(cl)
  precl <- foreach(timeindex = seq(1,900,12),
                   .combine = "abind3") %dopar% 
                   {
                     apply(precl1[,,timeindex:(timeindex+11)], MARGIN = c(1,2), sum)
                   }
  stopCluster(cl)
  gc()
  rm(precl1)
  precl
})





# make an average over 30 year periods
#mean_2006_


# make a SD over 30 year periods




