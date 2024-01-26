#Choisir l'environnement de tavail
WD = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(WD)

# Librairies
if (!require("readxl")) {install.packages("readxl")}+library("readxl") # API GBIF
if (!require("xlsx")) {install.packages("xlsx")}+library("xlsx")
if (!require("stringr")) {install.packages("stringr")}+library("stringr")




# Chargement des donnees

#PVF2

TYPO_PVF2_70 <- read_excel("../TYPO_PVF2_70_corrige.xlsx", sheet = "TYPO_PVF2_70")

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



##############TEST EXPRESSION REGULIERE


# Fonction pour détecter et corriger les abréviations de genre

phrase = "Nardus stricta, Achillea ptarmica subsp. pyrenaica, Agrostis rupestris subsp. rupestris, 
Carex umbrosa subsp. huetiana, Epikeros pyrenaeus (= Selinum pyrenaeum), Euphrasia minima subsp. m., 
Gentiana alpina, G. pyrenaica, Hieracium lactucella, Pedicularis pyrenaica, Ranunculus pyrenaeus, Trifolium alpinum, 
Selaginella selaginoides"


mots = strsplit(phrase, " |,|\n")[[1]]

test_vide <- function(x) {
  return(nchar(x) != 0)
}

mot_Genre = mots[str_detect(mots,"^[A-Z]")]


#Liste de mot dans ma phrase
mots = mots[unlist(lapply(mot, function(x) test_vide(x)))]

grepl( mots[6])

sub("^([A-Z]+)\\.$", "\\1", mots[6])

gsub("^([A-Z]+)\\s+.*", "\\1", mots[6])

# Fonction pour détecter et corriger les abréviations de genre dans une phrase
detecter_et_corriger_abreviation_phrase <- function(phrase) {
  mots <- strsplit(phrase, " ")[[1]]
  genre_precedent <- NULL
  
  for (i in seq_along(mots)) {
    mot <- mots[i]
    
    # Utilisation d'une expression régulière pour détecter les abréviations de genre
    if (grepl("^\\w+\\.$", mot)) {
      # Extraire l'abréviation de genre
      abrev_genre <- gsub("^([A-Z]+)\\.$", "\\1", mot)
      
      # Corriger l'abréviation en utilisant le genre précédent
      if (!is.null(genre_precedent)) {
        mots[i] <- paste(genre_precedent, gsub("^[A-Z]+\\.$", "", mot), sep = " ")
      }
    } else {
      # Mettre à jour le genre précédent si le mot n'est pas une abréviation
      genre_precedent <- gsub("^([A-Z]+)\\s+.*", "\\1", mot)
    }
  }
  
  return(paste(mots, collapse = " "))
}

# Liste de phrases avec des noms d'espèces
phrases <- c("Nardus stricta est une plante.", 
             "Achillea ptarmica subsp. pyrenaica est une autre espèce.",
             "Gentiana alpina, G. pyrenaica, et Hieracium lactucella sont des exemples.")

# Appliquer la fonction sur la liste de phrases
phrases_corrigees <- sapply(phrases, detecter_et_corriger_abreviation_phrase)

# Afficher les résultats
cat("Liste de phrases initiale:\n", paste(phrases, collapse = "\n"), "\n\n")
cat("Liste de phrases corrigée:\n", paste(phrases_corrigees, collapse = "\n"))
