setwd("C:/Users/Augustin Soulard/Documents/Programmation/Github/HABITATS/01_REF/PVF2/WORKFLOW")

# Chargement des données

#PVF2
library(readxl)
TYPO_PVF2_70 <- read_excel("../TYPO_PVF2_70.xlsx", sheet = "TYPO_PVF2_70")

#TAXREFv16
TAXREFv16_FLORE_FR = read.csv("../../../../BDD_FLORE_CONSTRUCT/TAXREF/TAXREFv16_FLORE_FR.csv",fileEncoding = "UTF-8")

# On conserve que les RANG especes, sous especes et variete
TAXREFv16_FLORE_FR = TAXREFv16_FLORE_FR[TAXREFv16_FLORE_FR$RANG %in% c("ES", "SSES", "VAR"),]

# Amorcage du tableau à remplir
TAB_ESP_PVF2 = data.frame(CD_HAB = "X",LB_CODE= "X",LB_HAB_FR="X",CD_NOM="X",LB_NOM="X")
TAB_ESP_PVF2 = TAB_ESP_PVF2[-1,]

# Boucle de rattachement TAXREF - PVF2
for(i in 1:nrow(TYPO_PVF2_70)){
  for( j in 1:nrow(TAXREFv16_FLORE_FR)){
    cat("PVF2 ligne : ",i," - TAXREF ligne : ",j,"\n")
    if(agrepl(TAXREFv16_FLORE_FR$LB_NOM[j],max.distance = 0.05,TYPO_PVF2_70$COMBINAISON_ESPECES[i],ignore.case = TRUE)){
      
      TAB_ESP_PVF2 = rbind(TAB_ESP_PVF2,c(TYPO_PVF2_70$CD_HAB[i],TYPO_PVF2_70$LB_CODE[i],TYPO_PVF2_70$LB_HAB_FR[i],TAXREFv16_FLORE_FR$CD_NOM[j],TAXREFv16_FLORE_FR$LB_NOM[j]))
      
    }
  }
}


 