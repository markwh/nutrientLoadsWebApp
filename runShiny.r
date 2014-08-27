

# setwd("shiny")
if(!require(shiny)) {install.packages("shiny"); library(shiny)}

options(shiny.trace = T)
runApp(display.mode = "showcase")


## To deploy to shinyapps:

library(devtools)
if(!require(shinyapps)){install_github("rstudio/shinyapps"); libary(shinyapps)}

shinyapps::accountInfo(name = "markwh")

deployApp()
authorizedUsers()
applications()
addAuthorizedUser()
shinyapps::terminateApp()

## Authorizing users:
if(!require(scrypt)){devtools::install_github('rstudio/rscrypt'); library(scrypt)}

addAuthorizedUser("mhpark")


deployApp()
