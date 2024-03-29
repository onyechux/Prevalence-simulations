#Reset the environment
rm(list = ls())

#Packages to be used
packages<-c("here","tidyverse","ggplot2","gridExtra","lme4","lmtest","readxl", "DT",
            "ggridges","viridis","hrbrthemes","tidyr","dplyr","forcats","ggpubr","merTools", 
            "RColorBrewer", "grid", "gtable" )

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))  

#Creating dir

dir.create(here("Figures"),showWarnings = F)
dir.create(here("Output"),showWarnings = F)


######################
# Optimizing cluster #
######################
source(here("Scripts","Ssq_stochastic_gradient.R"))
cluster(nrpl=10000,grad=1000)


##################
# Baseline model #
##################
source(here("Scripts","Baseline.R"))
baseline(nrpl=5000,n=56,c=0.61)

#############
# Scenarios #
#############
source(here("Scripts","Scenarios.R"))
scenarios(nrpl=5000,sce_n=c(10,33,56,79,102),sce_c=c(0.05,0.33,0.61,0.81,1.0))
