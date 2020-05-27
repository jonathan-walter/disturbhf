## Example for using disturbhf package
## 2020-05-28

library(devtools)
devtools::install_github("jonathan-walter/disturbhf") #only needs to be run the first time or after package update
library(disturbhf)

# Simulate some fake data to play with
refy<-data.frame(tt=1:365, yy=rnorm(365))
testy<-data.frame(tt=1:365, yy=rnorm(365))


test1<-mwdistdiffz()