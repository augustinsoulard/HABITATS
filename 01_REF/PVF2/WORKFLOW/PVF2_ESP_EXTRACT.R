#Choisir l'environnement de tavail
WD = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(WD)

# Librairies
if (!require("readxl")) {install.packages("readxl")}+library("readxl") # API GBIF
if (!require("xlsx")) {install.packages("xlsx")}+library("xlsx")
if (!require("stringr")) {install.packages("stringr")}+library("stringr")


###############Creation TYPO_PVF2_70_corrige
TYPO_PVF2_70_corrige <- read_excel("D:/Github/HABITATS/01_REF/PVF2/TYPO_PVF2_70_corrige.xlsx", 
                                   sheet = "TYPO_PVF2_70")

#Inition la typologie corrige
TYPO_PVF2_70_corrige$COMBINAISON_ESPECES_CORR = NA

#Boucle de correction
for(j in 1:nrow(TYPO_PVF2_70_corrige)){
  cat("Ligne : ",j,"\n")
  if(is.na(TYPO_PVF2_70_corrige$COMBINAISON_ESPECES[j])){
    next
  }
  mots = strsplit(TYPO_PVF2_70_corrige$COMBINAISON_ESPECES[j], " |,|\n")[[1]]
  
  lastGenre = "ERREUR"
  lastAdjectif = "ERREUR"
  
  for (i in 1:length(mots)) {
    #Traitement des genres
    if(str_detect(mots[i],"^[A-Z]\\.$") == TRUE){
      mots[i] = lastGenre
    }
    if(str_detect(mots[i],"^[A-Z]") == TRUE){
      lastGenre = mots[i]
    }
    #Traitement des infrataxons
    if(str_detect(mots[i],"\\.$") == FALSE){
      lastAdjectif = mots[i]
    }
    if(str_detect(mots[i],"^[a-z]\\.$") == TRUE){
      mots[i] = lastAdjectif
    }
  }
  TYPO_PVF2_70_corrige$COMBINAISON_ESPECES_CORR[j] = paste(mots,collapse = " ")
}
# Enregistrement en XLSX
write.csv(TYPO_PVF2_70_corrige,"TYPO_PVF2_70_corrige.csv",fileEncoding = "UTF-8",row.names = F)

# Chargement des donnees

#PVF2

#TYPO_PVF2_70 <- read_excel("../TYPO_PVF2_70_corrige.xlsx", sheet = "TYPO_PVF2_70")
TYPO_PVF2_70 = TYPO_PVF2_70_corrige

#TAXREFv16
TAXREFv16_FLORE_FR = read.csv("TAXREFv16_FLORE_FR.csv",fileEncoding = "UTF-8")

# On conserve que les RANG especes, sous especes et variete
TAXREFv16_FLORE_FR = TAXREFv16_FLORE_FR[TAXREFv16_FLORE_FR$RANG %in% c("ES", "SSES", "VAR"),]

# Amorcage du tableau ? remplir
TAB_ESP_PVF2 = data.frame(CD_HAB = character(0),LB_CODE= character(0),LB_HAB_FR=character(0),CD_NOM=character(0),LB_NOM=character(0))


# Boucle de rattachement TAXREF - PVF2 avec affichage
for(i in 1:nrow(TYPO_PVF2_70)){
  resultats <- sapply(TAXREFv16_FLORE_FR$LB_NOM, function(x) grepl(x,TYPO_PVF2_70$COMBINAISON_ESPECES_CORR[i]))
  TAXREFv16_FLORE_FR[resultats,]$LB_NOM
  
    cat("PVF2 ligne : ",i,"\n")
    if(nrow(data.frame(TAXREFv16_FLORE_FR[resultats,]$LB_NOM))>0){
      
      TAB_TEMP = data.frame(CD_HAB = TYPO_PVF2_70$CD_HAB[i],LB_CODE= TYPO_PVF2_70$LB_CODE[i],LB_HAB_FR=TYPO_PVF2_70$LB_HAB_FR[i],CD_NOM=TAXREFv16_FLORE_FR[resultats,]$CD_REF,LB_NOM=TAXREFv16_FLORE_FR[resultats,]$LB_NOM)
      
      TAB_ESP_PVF2 = rbind(TAB_ESP_PVF2,TAB_TEMP)
 
  }
}

# Enregistrement en XLSX
write.csv(TAB_ESP_PVF2,"TAB_ESP_PVF2.csv",fileEncoding = "UTF-8",row.names = F)



