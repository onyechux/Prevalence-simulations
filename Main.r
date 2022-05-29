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

#Creating dir.

dir.create(here("Figures"))
dir.create(here("Output"))

######################
# Optimizing cluster #
######################


##################
# Baseline model #
##################
source(here("scripts","2Baseline.R"))
baseline(nrpl=50)

#############
# Scenarios #
#############
