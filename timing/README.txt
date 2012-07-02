##note: timing.Rnw is NOT a Sweave file - it's a knitr file. Calling Sweave
##on it probably won't work and knitr should be used instead. The following
##script will create the tex file, which can be compiled normally like any
##tex file:

library(knitr) ##install from CRAN using install.packages("knitr") if necessary
knit("timing.Rnw") ##creates tex file
