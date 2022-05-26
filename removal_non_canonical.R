
## Con este script removemos los cromosomas no canonicos
## Creamos el ambiente
## asumimos que se ha descargado desde UCSC
version <- "38"
genoma <- paste0("genomaHg",version)
archivo<- paste0("../datos/",genoma,"/genoma/hg",version,".fa")
carpeta_salida <- paste0("../datos/",genoma,"/genoma/canonical")
archivo_salida <- file.path(carpeta_salida,paste0(genoma,".fa"))
comando <- paste("./remove_non_canonical.sh",archivo,">",archivo_salida)

if(!dir.exists(carpeta_salida)){
  
  dir.create(carpeta_salida)
  print(comando)
  system(comando)
  
}else if(dir.exists(carpeta_salida) && !file.exists(archivo_salida)){
  print(comando)
  system(comando)
  
}else{
  message("Ya se ha canonizado")
}


