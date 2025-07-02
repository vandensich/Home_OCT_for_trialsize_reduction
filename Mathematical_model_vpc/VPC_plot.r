rm(list=ls())
library(dplyr)
library(nonmem2R)

model <- 'HomeOCT_CST'

if (!file.exists("figures")) {dir.create ("figures")}
  
fname <- paste("figures/",model,"_VPCs.pdf", sep="")
if (file.exists (fname)){ file.remove (fname) }

pdf (file = fname, width=7, height = 5)

p1 <- vpcfig2(vpctab="vpctab",vpcresult="vpc_results.csv")  + xlab("Time (days)") + ylab("Prediction corrected PD")

print(p1)


dev.off()

# open created file
# if (file.exists(fname)) {
#   if (Sys.info()['sysname'] == 'Windows') { shell.exec(paste(getwd(),"/",fname,sep="")) }  # windows
#   else if (Sys.info()['sysname'] == 'Darwin') { system(paste ("open ",fname, sep="")) } # mac
#   else { system(paste ("xdg-open ",fname, sep="")) } # linux
# }