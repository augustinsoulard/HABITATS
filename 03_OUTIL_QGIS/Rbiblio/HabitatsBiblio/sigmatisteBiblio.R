#Choisir l'environnement de tavail
WD = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(WD)

# Librairies
if (!require("sf")) {install.packages("sf")}+library("sf") # gestion carto
if (!require("tidyverse")) {install.packages("tidyverse")}+library("tidyverse")
if (!require("readxl")) {install.packages("readxl")}+library("readxl")

# Définition des variables
nomShp = "CLC2018_CUT"
nomPointFlore = "POINT_2024_01_25_16h17m44"


#Fonction de traitement des données


#Charger le fichier de TAB_ESP_PVF2
TAB_ESP_PVF2 = read.csv("TAB_ESP_PVF2.csv")

# Charger les fichiers Shapefile
OCCSOL <- st_read(paste0("../../Occsol/",nomShp,".shp"))
Point_Biblio_Flore <- st_read(paste0("../../FloreBiblio/",nomPointFlore,".shp"))

# Effectuer l'intersection
intersection_resultat <- st_intersection(Point_Biblio_Flore, OCCSOL)

#Garder uniquement les colonnes utiles a l'analyse
OBS_DATA_POLYGONE = intersection_resultat %>% select(cd_ref,nom_valide,ID,Code_18)

# Creation du tableau avec les habitats répété par ID de polygone et par occsol
OBS_DATA_ESP_PHYTO = merge(OBS_DATA_POLYGONE,TAB_ESP_PVF2,by.x = "cd_ref",by.y = "CD_NOM",all = TRUE)
OBS_DATA_ESP_PHYTO = OBS_DATA_ESP_PHYTO[!is.na(OBS_DATA_ESP_PHYTO$ID),]
OBS_DATA_ESP_PHYTO$COUNT = 1

# Calcule des points pour la typologie OCCSOL
POINT_TYPO_OCCSOL = aggregate(COUNT~Code_18+CD_HAB+LB_CODE+LB_HAB_FR,data = OBS_DATA_ESP_PHYTO,FUN = sum)

# Calcule des points pour l'ID OCCSOL
POINT_ID_OCCSOL = aggregate(COUNT~ID+Code_18+CD_HAB+LB_CODE+LB_HAB_FR,data = OBS_DATA_ESP_PHYTO,FUN = sum)

#Jointure pour avoir des points d'occsol et d'ID de typologie
POINT_TYPO_OCCSOL$KEY = paste0(POINT_TYPO_OCCSOL$Code_18,POINT_TYPO_OCCSOL$CD_HAB)
POINT_ID_OCCSOL$KEY = paste0(POINT_ID_OCCSOL$Code_18,POINT_ID_OCCSOL$CD_HAB)

POINT_ID_TYPO_OCCSOL = left_join(POINT_ID_OCCSOL,POINT_TYPO_OCCSOL,by="KEY")

#CALCULE des points par habitats
POINT_ID_TYPO_OCCSOL$TOTAL = POINT_ID_TYPO_OCCSOL$COUNT.x*5 +POINT_ID_TYPO_OCCSOL$COUNT.y

#Garder les 3 habitats les plus probables
POINT_ID_TYPO_OCCSOL_FILTRE <- POINT_ID_TYPO_OCCSOL %>%
  group_by(ID) %>%
  mutate(Rang = dense_rank(desc(TOTAL))) %>%
  filter(Rang <= 3) %>%
  select(-Rang) %>%
  ungroup()

# Faire le rattachement avec les autres typologies d'habitats
# Chargement des autres typologies
CRSP_PVF2_EUNIS_2012_70 <- read_excel("CRSP_PVF2_EUNIS_2012_70.xlsx", sheet = "CRSP_PVF2_EUNIS_2012_70")
CRSP_PVF2_CORINE_BIOTOPES_70 <- read_excel("CRSP_PVF2_CORINE_BIOTOPES_70.xlsx", sheet = "CRSP_PVF2_CORINE_BIOTOPES_70")
CRSP_PVF2_HIC_70 <- read_excel("CRSP_PVF2_HIC_70.xlsx", sheet = "CRSP_PVF2_HIC_70")

#Selection des colonnes
CRSP_PVF2_EUNIS_2012_70 = CRSP_PVF2_EUNIS_2012_70 %>% select(CD_HAB_ENTRE,CD_HAB_SORTIE,LB_CODE_SORTIE,LB_HAB_SORTIE)
CRSP_PVF2_CORINE_BIOTOPES_70 = CRSP_PVF2_CORINE_BIOTOPES_70 %>% select(CD_HAB_ENTRE,CD_HAB_SORTIE,LB_CODE_SORTIE,LB_HAB_SORTIE)
CRSP_PVF2_HIC_70 = CRSP_PVF2_HIC_70 %>% select(CD_HAB_ENTRE,CD_HAB_SORTIE,LB_CODE_SORTIE,LB_HAB_SORTIE)

# Jointure des typologies
POINT_ID_TYPO_OCCSOL_FILTRE$CD_HAB = as.character(POINT_ID_TYPO_OCCSOL_FILTRE$CD_HAB.y)
JOINTURE_BRUT = left_join(POINT_ID_TYPO_OCCSOL_FILTRE,CRSP_PVF2_EUNIS_2012_70,by=c("CD_HAB"="CD_HAB_ENTRE"))
JOINTURE_BRUT = left_join(JOINTURE_BRUT,CRSP_PVF2_CORINE_BIOTOPES_70,by=c("CD_HAB"="CD_HAB_ENTRE"))
JOINTURE_BRUT = left_join(JOINTURE_BRUT,CRSP_PVF2_HIC_70,by=c("CD_HAB"="CD_HAB_ENTRE"))

#Selection et renommage des colonnes
JOINTURE_BRUT = JOINTURE_BRUT %>%
  select(ID,OCCSOL = Code_18.x,CD_HAB,CD_PVF2 = LB_CODE.x,LB_PVF2 = LB_HAB_FR.x,
         CD_EUNIS2012=LB_CODE_SORTIE.x,LB_EUNIS2012=LB_HAB_SORTIE.x,
         CD_CORINE_BIOTOPE=LB_CODE_SORTIE.y,LB_CORINE_BIOTOPE=LB_HAB_SORTIE.y,
         CD_HIC=LB_CODE_SORTIE,LB_HIC=LB_HAB_SORTIE,
         POINT=TOTAL)

#Enregistrement de la sorte en csv
write.csv(JOINTURE_BRUT,"JOINTURE_BRUT.csv",fileEncoding = "UTF-8",row.names = F)

#Selection des colonnes de OCCSOL
OCCSOL = OCCSOL %>% select(ID,geometry)

#Jointure à OCCSOL
OCCSOL_COMPLET = merge(OCCSOL,JOINTURE_BRUT,by = "ID",all = TRUE)
OCCSOL_COMPLET_2D  = st_zm(OCCSOL_COMPLET, drop = TRUE)

# Ecriture de la couche finale
st_write(OCCSOL_COMPLET_2D, "../../Occsol/OCCSOL_COMPLET3.shp")
