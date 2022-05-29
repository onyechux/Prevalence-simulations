#Reset the environment
rm(list = ls())

#Packages to be used
packages<-c("here","tidyverse","ggplot2","gridExtra","lme4","lmtest","readxl", "DT",
            "ggridges","viridis","hrbrthemes","tidyr","dplyr","forcats")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))  

#Creating dir

dir.create(here("Figures"))
dir.create(here("Output"))

######################
# Optimizing cluster #
######################
source(here("Scripts","Ssq_stochastic_gradient.R"))
cluster(rpl=5000,grad=300)


##################
# Baseline model #
##################
source(here("Scripts","Baseline.R"))
baseline(nrpl=50)

#############
# Scenarios #
#############
source(here("Scripts","Scenarios.R"))
scenarios(nrpl=5000)
