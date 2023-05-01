library(foreign)
DATAPHYTO = read.csv(choose.files())

library(reshape2)

#Décroiser le tableau
MELT_DATA <- melt(DATAPHYTO, id.vars = c("CD_NOM","NOM_VALIDE","RELEVE"), measure.vars = c("ARBO", "ARBU", "HERB", "MUCINALE"))

# Attribution d'un poid pour chaque classe phyto
library(tidyverse)
MELT_DATA$value = str_replace_all(MELT_DATA$value,c("5"="10","4"="8","3"="6","2"="4","1"="1","\\+"="0.5","r"="0.25","i"="0.125"))
MELT_DATA$value = as.numeric(MELT_DATA$value)

#Somme des scores par espèce
AGG_DATA = aggregate(value~CD_NOM+NOM_VALIDE+RELEVE,data=MELT_DATA,sum)

# Joindre les informations de baseflor
TAXREF_baseflor <- read.csv("~/R/Rphytosocio/TAXREF_baseflor.csv", sep=";")

JOIN_DATA  = left_join(AGG_DATA,TAXREF_baseflor,by ="CD_NOM")

JOIN_DATA_CALL = select(JOIN_DATA,RELEVE,CD_NOM,NOM_VALIDE,Nom.scientifique, value,code_CATMINAT)


DATARELEVE = data.frame(RELEVE = levels(factor(JOIN_DATA_CALL$RELEVE)),CATMINAT = "")
#Boucle qui tourne pour chaque relevé
for(i in 1:nrow(DATARELEVE)){
  JOIN_DATA_CALL = JOIN_DATA_CALL[JOIN_DATA_CALL$RELEVE == DATARELEVE$RELEVE[i],]
  CATMINAT_SEP = separate(JOIN_DATA_CALL, code_CATMINAT, into = c("code1", "code2", "code3", "code4", "code5","code6","code7"), sep = "/|\\.")
  JOIN_DATA_CALL = left_join(CATMINAT_SEP,JOIN_DATA_CALL,by="CD_NOM",keep = FALSE)[,-c(13:16)]
  # On analyse les code CATMINAT dans l'ordre d'importance
  for(j in 1:7){
    cat(j)
    col = paste0("code",j)
    if(all(is.na(JOIN_DATA_CALL[col]))){break}
    AGG = aggregate(value.x~get(paste0("code",j)),data=JOIN_DATA_CALL,sum)
    sort_value = sort(AGG$value.x,decreasing = T)
    if(sort_value[1]<sort_value[2]*1.2){
      cat("PROXIMITE OU EGALITE DES VALEURS AU NIVEAU ",j)
      break}
    
    #Attribution du code CATMINAT selon l'avancement
    CATMIN = as.character(AGG[AGG$value.x==max(AGG$value.x),][1])
    if(j == 1){sep1 = "/" ; sep2 = ""} 
    else if(j == 2){sep1 = "." ; sep2 = ""} 
    else if(j == 3){sep1 = "" ; sep2 = ""} 
    else{sep1 = "" ; sep2 = "."}
    DATARELEVE$CATMINAT[i] = paste0(DATARELEVE$CATMINAT[i],sep2,CATMIN,sep1)
    JOIN_DATA_CALL = JOIN_DATA_CALL[JOIN_DATA_CALL[col]== CATMIN,]

    }
  }

# Une fois un code attribué à chaque relevé, on y joint les informations de baseveg
baseveg <- read.csv("~/R/Rphytosocio/baseveg.csv", sep=";")
DATARELEVE_JOIN = left_join(DATARELEVE,baseveg,by=c("CATMINAT"="CODE..CATMINAT"))

