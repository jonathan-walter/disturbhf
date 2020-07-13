#rm(list=ls())
library(roxygen2)
library(devtools)

setwd("~/GitHub")#/disturbhf/")

#create("disturbhf")
build("disturbhf")   # need a folder named 'wtCNN, in which a file named DESCRIPTION and a sub-folder named 'R'. This will creat a zip file
document("disturbhf") # this will create a file 'NAMESPACE' and a folder 'man'
check("disturbhf")


#install("disturbhf")
#library("disturbhf")
install_github("jonathan-walter/disturbhf")
