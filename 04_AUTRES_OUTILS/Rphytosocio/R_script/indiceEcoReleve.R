# Chargement de las base de donnees base flore avec CD_NOM
baseflorTAXREFv16 = read.csv("baseflorTAXREFv16.csv", sep=";")


# Jointure de baseflor a DATAOPHYTO
DATAbaseflorJOIN = left_join(DATAPHYTO,baseflorTAXREFv16,by="CD_NOM")

#Verification de la jointure et NA
if(any(is.na(DATAbaseflorJOIN$N._Nomenclatural_BDNFF)==TRUE)){
  cat("!!!ATTENTION JOINTURE INCOMPLETE !!!")
  
}

DATAbaseflorJOIN[355,]
