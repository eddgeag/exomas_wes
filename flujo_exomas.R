
## Flujo bioinfomartico para exomas
## Basado en subprocesos en R

## Primero nos descargamos la versio  hg38 de UCSC 

## Y corremos en el terminal el script sh para remover cromosomas no canonicos

## Creamos las funciones de calidad
## OJO  no creamos trimmomatic, en caso que sea necesario, pero por 
## el momento descartamos

## Creamos el ambiente

library(tools)

pipeline<-"../pipeline"
muestra <- "muestraToy"
folder_fq <- "fastq_files"
folder_fasta <-"../datos/cromosoma8"

files_folder <- file.path(pipeline,muestra,folder_fq)



fastqc_R <-
  function(input_directory=files_folder) {
    
    exit_status1 <-1
    exit_status0 <-0
    
    files_fastq <- list.files(input_directory,full.names = T)
    
    pattern_files <- unique(file_ext(files_fastq))
    
    output_dir <- file.path(paste(unlist(strsplit(file.path(input_directory),"/"))[1:3],collapse = "/")
                            ,"output_dir")
    output_dir_fqc <- file.path(output_dir,"output_QC")
    command <- paste("fastqc -t 4 ",paste0(input_directory,"/*.",pattern_files), "-o", output_dir_fqc)
    print(command)
    
    if(!dir.exists(output_dir) && !dir.exists(output_dir_fqc)){
      dir.create(output_dir)
      dir.create(output_dir_fqc)
      
      system(command)
      
      return(exit_status1)
      
    }else if(!dir.exists(output_dir_fqc) && dir.exists(output_dir)){
      
      dir.create(output_dir_fqc)
      system(command)
      return(exit_status1)
      
    }else if(dir.exists(output_dir_fqc) && dir.exists(output_dir)){
      
      if(length(list.files(output_dir_fqc))==0){
        
        system(command)
        
      }else{
        
        message("Ya se a hecho el control de calidad 1")
        
      }
      
    }
    
    
  }


index_fasta_samtools <- function(input_directory=folder_fasta){
  
  extension = unlist(lapply(list.files(input_directory,pattern="fa"),function(x) file_ext(x)))
  if(extension[1]=="fasta"||extension[1]=="fa"){
    extension<-extension
  }else{
    errorCondition("No existe archivo fasta")
  }
  fasta_file <- list.files(input_directory,pattern = paste0(".",extension),full.names = T)

  if(length(fasta_file)>1){
    stop("there are more than one file")
  }else if(length(fasta_file)==0){
    stop("there are not fasta files!")
  }
  
  fasta_file <- fasta_file[1]

  command <- paste("samtools faidx",fasta_file)
  print(command)
  system(command)
  
  
  
}

index_fasta_bwa <- function(input_directory=folder_fasta){
  
  extension = "fasta$"
  extension = unlist(lapply(list.files(input_directory,pattern="fa"),function(x) file_ext(x)))
  if(extension[1]=="fasta"||extension[1]=="fa"){
    extension<-extension
  }else{
    errorCondition("No existe archivo fasta")
  }
  fasta_file <- list.files(input_directory,pattern = paste0(".",extension),full.names = T)

  
  fasta_file <- fasta_file[1]

  command <- paste("bwa index",fasta_file)

  print(command)
  
  system(command)
  
}



bwa_mem2 <- function(fastq_folder=files_folder ,reference_genome =folder_fasta){
  
  ## Comprobamos que existe archivo fasta
  
  extension = "fasta$"
  extension = unlist(lapply(list.files(reference_genome,pattern="fa"),function(x) file_ext(x)))
  if(extension[1]=="fasta"||extension[1]=="fa"){
    extension<-extension
  }else{
    errorCondition("No existe archivo fasta")
  }
  fasta_file <- list.files(reference_genome,pattern = paste0(".",extension),full.names = T)
  
  ##### fai index exist ?######
  
  if(!length(file.exists(file.path(reference_genome,list.files(reference_genome,"fai"))))){
    
    print("### generating fai index...")
    
    index_fasta_samtools(input_directory=reference_genome)
    
  }
  
  ## Comprobamos que esiste bwa index files
  extension <- c("amb","ann","bwt","pac","sa" )
  
  if(!length(file.exists(list.files(reference_genome,extension[1]))) |
     !length(file.exists(list.files(reference_genome,extension[2])) )|
     !length(file.exists(list.files(reference_genome,extension[3]))) |
     !length(file.exists(list.files(reference_genome,extension[4]))) |
     !length(file.exists(list.files(reference_genome,extension[4])) )){
    ## no existen bwa index
    print("Creando ficheros índices para bwa mem...")
    index_fasta_bwa(input_directory=reference_genome)
    
  }
  

  
  ## Buscamos refenrecia fasta
  extensiones <- unlist(lapply(fasta_file,function(x) file_ext(x)))
  tmp.ext <- grep("fa$",extensiones)
  tmp.ext2 <- grep("fasta$",extensiones)
  extensiones_which<- c(tmp.ext,tmp.ext2)
  extension <- extensiones[extensiones_which]
  fasta_file <- list.files(reference_genome,full.names=T,pattern=paste0(extension,"$"))
  fastq_files<- list.files(fastq_folder,full.names=T)
  ## Archivos fastq
  fastq_full_path_files <- list(fastq_1 = fastq_files[1],fastq_2 = fastq_files[2])
  ## Creamos directorio de mapeo
  output_file_name<-unlist(strsplit(gsub("R[12]","map",fastq_files[1]),"/"))[5]
  
  output_dir <- file.path(paste(unlist(strsplit(file.path(fastq_folder),"/"))[1:3],collapse = "/")
                          ,"output_dir")
  output_folder<- file.path(output_dir,"mapping_output")
  if(!dir.exists(output_folder)){
    dir.create(output_folder)
    if(length(list.files(output_folder))==0)

      output_file_sam <- file.path(output_folder,paste0(output_file_name,".sam"))
      output_file_bam <- file.path(output_folder,paste0(output_file_name,".bam"))
      output_file_sorted_bam <- file.path(output_folder,paste0(output_file_name,".sorted.bam"))

 
      print("#### MAPPING...#####")
      
      comando_mapper <- paste("./bwa-mapper.sh ",fastq_full_path_files$fastq_1,
                              fastq_full_path_files$fastq_2,
                              fasta_file,output_file_bam)
      print(comando_mapper)
      system(comando_mapper)
          ## crea el mapeo

        command_out_bam_sorted <- paste("samtools sort ",output_file_bam ,"-o",output_file_sorted_bam)
        print(command_out_bam_sorted)
  
        system(command = command_out_bam_sorted)
  
        return(0)

  }else if(dir.exists(output_folder) && !dir.exists(output_folder)){
    dir.create(output_folder)
    output_file_sam <- file.path(output_folder,paste0(output_file_name,".sam"))
    output_file_bam <- file.path(output_folder,paste0(output_file_name,".bam"))
    output_file_sorted_bam <- file.path(output_folder,paste0(output_file_name,".sorted.bam"))
    
    
    print("#### MAPPING...#####")
    
    comando_mapper <- paste("./bwa-mapper.sh ",fastq_full_path_files$fastq_1,
                            fastq_full_path_files$fastq_2,
                            fasta_file,output_file_bam)
    print(comando_mapper)
    system(comando_mapper)
    ## crea el mapeo
    
    command_out_bam_sorted <- paste("samtools sort ",output_file_bam ,"-o",output_file_sorted_bam)
    print(command_out_bam_sorted)
    
    system(command = command_out_bam_sorted)
    
    
  }else if(dir.exists(output_folder) && dir.exists(output_folder)){
    
    if(length(list.files(output_folder))==0){
      
      output_file_sam <- file.path(output_folder,paste0(output_file_name,".sam"))
      output_file_bam <- file.path(output_folder,paste0(output_file_name,".bam"))
      output_file_sorted_bam <- file.path(output_folder,paste0(output_file_name,".sorted.bam"))
      
      
      print("#### MAPPING...#####")
      
      comando_mapper <- paste("./bwa-mapper.sh ",fastq_full_path_files$fastq_1,
                              fastq_full_path_files$fastq_2,
                              fasta_file,">",output_file_bam)
      print(comando_mapper)
      system(comando_mapper)
      ## crea el mapeo
      
      command_out_bam_sorted <- paste("samtools sort ",output_file_bam ,"-o",output_file_sorted_bam)
      print(command_out_bam_sorted)
      
      system(command = command_out_bam_sorted)
      
      
    }else{
      
      message("Ya se ha mapeado")
    }
    
    
  }

  }
  




bwa_mem2()
