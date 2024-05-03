# Definition du repertoire de fichiers
WD = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(WD)

# Charger les bibliotheques necessaires
if(!require("readxl")){install.packages("readxl")} ; library("readxl")
if(!require("tidyverse")){install.packages("tidyverse")} ; library("tidyverse")



HABREF_70 <- read.csv("../HABREF_70.csv", sep=";")
HABREF_CORRESP_HAB_70 <- read.csv("../HABREF_CORRESP_HAB_70.csv", sep=";")
HABREF_CORRESP_TAXON_70 <- read.csv("../HABREF_CORRESP_TAXON_70.csv", sep=";")
HABREF_DESCRIPTION_70 <- read.csv("../HABREF_DESCRIPTION_70.csv", sep=";")


TYPOREF_70 <- read.csv("../TYPOREF_70.csv", sep=";")


PVF2 = HABREF_70 %>% filter(CD_TYPO =="28")

PVF2_TAXON = left_join(PVF2,HABREF_CORRESP_TAXON_70,by=c("CD_HAB"="CD_HAB_ENTRE"))
