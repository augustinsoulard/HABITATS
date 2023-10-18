#Choisir l'environnement de tavail
WD = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(WD)

# Librairies
if (!require("readxl")) {install.packages("readxl")}+library("readxl") # API GBIF
if (!require("xlsx")) {install.packages("xlsx")}+library("xlsx")
if (!require("stringr")) {install.packages("stringr")}+library("stringr")




# Chargement des donnees

#PVF2

TYPO_PVF2_70 <- read_excel("../TYPO_PVF2_70.xlsx", sheet = "TYPO_PVF2_70")

#TAXREFv16
TAXREFv16_FLORE_FR = read.csv("TAXREFv16_FLORE_FR.csv",fileEncoding = "UTF-8")

# On conserve que les RANG especes, sous especes et variete
TAXREFv16_FLORE_FR = TAXREFv16_FLORE_FR[TAXREFv16_FLORE_FR$RANG %in% c("ES", "SSES", "VAR"),]

# Amorcage du tableau ? remplir
TAB_ESP_PVF2 = data.frame(CD_HAB = character(0),LB_CODE= character(0),LB_HAB_FR=character(0),CD_NOM=character(0),LB_NOM=character(0))


# Boucle de rattachement TAXREF - PVF2 avec affichage
for(i in 1:nrow(TYPO_PVF2_70)){
  resultats <- sapply(TAXREFv16_FLORE_FR$LB_NOM, function(x) grepl(x,TYPO_PVF2_70$COMBINAISON_ESPECES[i]))
  TAXREFv16_FLORE_FR[resultats,]$LB_NOM
  
    cat("PVF2 ligne : ",i,"\n")
    if(nrow(data.frame(TAXREFv16_FLORE_FR[resultats,]$LB_NOM))>0){
      
      TAB_TEMP = data.frame(CD_HAB = TYPO_PVF2_70$CD_HAB[i],LB_CODE= TYPO_PVF2_70$LB_CODE[i],LB_HAB_FR=TYPO_PVF2_70$LB_HAB_FR[i],CD_NOM=TAXREFv16_FLORE_FR[resultats,]$LB_NOM,LB_NOM=TAXREFv16_FLORE_FR[resultats,]$CD_REF)
      
      TAB_ESP_PVF2 = rbind(TAB_ESP_PVF2,TAB_TEMP)
 
  }
}

# Enregistrement en XLSX
write.xlsx(TAB_ESP_PVF2,"TAB_ESP_PVF2.xlsx")


########WORKFLOW####################

any(str_detect(TAXREFv16_FLORE_FR$LB_NOM, TYPO_PVF2_70$COMBINAISON_ESPECES[14]))

