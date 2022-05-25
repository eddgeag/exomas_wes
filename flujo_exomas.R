

## Flujo bioinfomartico para exomas enfermedades geneticas (1)
## Basado en subprocesos en R

## Primero nos descargamos la versio  hg38 de UCSC

## Y corremos en el terminal el script sh para remover cromosomas no canonicos

## Creamos las funciones de calidad
## OJO  no creamos trimmomatic, en caso que sea necesario, pero por
## el momento descartamos

## Creamos el ambiente

#### FATLA CREAR QUE SEA UN PIPELINE



library(tools)

pipeline <- "../pipeline"
muestra <- "muestraToy"
folder_fq <- "fastq_files"
folder_fasta <- "../datos/cromosoma8"

files_folder <- file.path(pipeline, muestra, folder_fq)
output_directory <- file.path(pipeline, muestra, "output_dir")


fastqc_R <-
  function(input_directory = files_folder) {
    exit_status1 <- 1
    exit_status0 <- 0
    
    files_fastq <- list.files(input_directory, full.names = T)
    
    pattern_files <- unique(file_ext(files_fastq))
    
    output_dir <-
      file.path(paste(unlist(strsplit(
        file.path(input_directory), "/"
      ))[1:3], collapse = "/")
      , "output_dir")
    output_dir_fqc <- file.path(output_dir, "output_QC")
    command <-
      paste("fastqc -t 4 ",
            paste0(input_directory, "/*.", pattern_files),
            "-o",
            output_dir_fqc)
    print(command)
    
    if (!dir.exists(output_dir) && !dir.exists(output_dir_fqc)) {
      dir.create(output_dir)
      dir.create(output_dir_fqc)
      
      system(command)
      
      return(exit_status1)
      
    } else if (!dir.exists(output_dir_fqc) &&
               dir.exists(output_dir)) {
      dir.create(output_dir_fqc)
      system(command)
      return(exit_status1)
      
    } else if (dir.exists(output_dir_fqc) &&
               dir.exists(output_dir)) {
      if (length(list.files(output_dir_fqc)) == 0) {
        system(command)
        
      } else{
        message("Ya se a hecho el control de calidad 1")
        
      }
      
    }
    
    
  }


index_fasta_samtools <- function(input_directory = folder_fasta) {
  extension = unlist(lapply(list.files(input_directory, pattern = "fa"), function(x)
    file_ext(x)))
  if (extension[1] == "fasta" || extension[1] == "fa") {
    extension <- extension
  } else{
    errorCondition("No existe archivo fasta")
  }
  fasta_file <-
    list.files(input_directory,
               pattern = paste0(".", extension),
               full.names = T)
  
  if (length(fasta_file) > 1) {
    stop("there are more than one file")
  } else if (length(fasta_file) == 0) {
    stop("there are not fasta files!")
  }
  
  fasta_file <- fasta_file[1]
  
  command <- paste("samtools faidx", fasta_file)
  print(command)
  system(command)
  
  
  
}

index_fasta_bwa <- function(input_directory = folder_fasta) {
  extension = "fasta$"
  extension = unlist(lapply(list.files(input_directory, pattern = "fa"), function(x)
    file_ext(x)))
  if (extension[1] == "fasta" || extension[1] == "fa") {
    extension <- extension
  } else{
    errorCondition("No existe archivo fasta")
  }
  fasta_file <-
    list.files(input_directory,
               pattern = paste0(".", extension),
               full.names = T)
  
  
  fasta_file <- fasta_file[1]
  
  command <- paste("bwa index", fasta_file)
  
  print(command)
  
  system(command)
  
}



bwa_mem2 <-
  function(fastq_folder = files_folder ,
           reference_genome = folder_fasta) {
    ## Comprobamos que existe archivo fasta
    
    extension = "fasta$"
    extension = unlist(lapply(list.files(reference_genome, pattern = "fa"), function(x)
      file_ext(x)))
    if (extension[1] == "fasta" || extension[1] == "fa") {
      extension <- extension
    } else{
      errorCondition("No existe archivo fasta")
    }
    fasta_file <-
      list.files(reference_genome,
                 pattern = paste0(".", extension),
                 full.names = T)
    
    ##### fai index exist ?######
    
    if (!length(file.exists(file.path(
      reference_genome, list.files(reference_genome, "fai")
    )))) {
      print("### generating fai index...")
      
      index_fasta_samtools(input_directory = reference_genome)
      
    }
    
    ## Comprobamos que esiste bwa index files
    extension <- c("amb", "ann", "bwt", "pac", "sa")
    
    if (!length(file.exists(list.files(reference_genome, extension[1]))) |
        !length(file.exists(list.files(reference_genome, extension[2]))) |
        !length(file.exists(list.files(reference_genome, extension[3]))) |
        !length(file.exists(list.files(reference_genome, extension[4]))) |
        !length(file.exists(list.files(reference_genome, extension[4])))) {
      ## no existen bwa index
      print("Creando ficheros índices para bwa mem...")
      index_fasta_bwa(input_directory = reference_genome)
      
    }
    
    
    
    ## Buscamos refenrecia fasta
    extensiones <- unlist(lapply(fasta_file, function(x)
      file_ext(x)))
    tmp.ext <- grep("fa$", extensiones)
    tmp.ext2 <- grep("fasta$", extensiones)
    extensiones_which <- c(tmp.ext, tmp.ext2)
    extension <- extensiones[extensiones_which]
    fasta_file <-
      list.files(reference_genome,
                 full.names = T,
                 pattern = paste0(extension, "$"))
    fastq_files <- list.files(fastq_folder, full.names = T)
    ## Archivos fastq
    fastq_full_path_files <-
      list(fastq_1 = fastq_files[1], fastq_2 = fastq_files[2])
    ## Creamos directorio de mapeo
    output_file_name <-
      unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))[5]
    
    output_dir <-
      file.path(paste(unlist(strsplit(
        file.path(fastq_folder), "/"
      ))[1:3], collapse = "/")
      , "output_dir")
    output_folder <- file.path(output_dir, "mapping_output")
    output_file_sam <-
      file.path(output_folder, paste0(output_file_name, ".sam"))
    output_file_bam <-
      file.path(output_folder, paste0(output_file_name, ".bam"))
    output_file_sorted_bam <-
      file.path(output_folder, paste0(output_file_name, ".sorted.bam"))
    if (!dir.exists(output_folder)) {
      dir.create(output_folder)
      if (length(list.files(output_folder)) == 0)
        
        
        print("#### MAPPING...#####")
      
      comando_mapper <-
        paste(
          "./bwa-mapper.sh ",
          fastq_full_path_files$fastq_1,
          fastq_full_path_files$fastq_2,
          fasta_file,
          output_file_sam
        )
      print(comando_mapper)
      system(comando_mapper)
      
      print("#### SAM TO BAM")
      
      command_sam_to_bam <-
        paste("samtools view -S -b",
              output_file_sam,
              "-o",
              output_file_bam)
      
      print(command_sam_to_bam)
      system(command_sam_to_bam)
      
      print("### BAM to sorted BAM")
      command_out_bam_sorted <-
        paste("samtools sort ",
              output_file_bam ,
              "-o",
              output_file_sorted_bam)
      
      print(command_out_bam_sorted)
      system(command_out_bam_sorted)
      
      
    } else if (dir.exists(output_folder) &&
               !dir.exists(output_folder)) {
      dir.create(output_folder)
      
      
      
      
      print("#### MAPPING...#####")
      
      comando_mapper <-
        paste(
          "./bwa-mapper.sh ",
          fastq_full_path_files$fastq_1,
          fastq_full_path_files$fastq_2,
          fasta_file,
          output_file_sam
        )
      print(comando_mapper)
      system(comando_mapper)
      
      print("#### SAM TO BAM")
      
      command_sam_to_bam <-
        paste("samtools view -S -b",
              output_file_sam,
              "-o",
              output_file_bam)
      
      print(command_sam_to_bam)
      system(command_sam_to_bam)
      
      print("### BAM to sorted BAM")
      command_out_bam_sorted <-
        paste("samtools sort ",
              output_file_bam ,
              "-o",
              output_file_sorted_bam)
      
      print(command_out_bam_sorted)
      system(command_out_bam_sorted)
      
      
    } else if (dir.exists(output_folder) &&
               dir.exists(output_folder)) {
      if (length(list.files(output_folder)) == 0) {
        print("#### MAPPING...#####")
        
        comando_mapper <-
          paste(
            "./bwa-mapper.sh ",
            fastq_full_path_files$fastq_1,
            fastq_full_path_files$fastq_2,
            fasta_file,
            output_file_sam
          )
        print(comando_mapper)
        system(comando_mapper)
        
        print("#### SAM TO BAM")
        
        command_sam_to_bam <-
          paste("samtools view -S -b",
                output_file_sam,
                "-o",
                output_file_bam)
        
        print(command_sam_to_bam)
        system(command_sam_to_bam)
        
        print("### BAM to sorted BAM")
        command_out_bam_sorted <-
          paste("samtools sort ",
                output_file_bam ,
                "-o",
                output_file_sorted_bam)
        
        print(command_out_bam_sorted)
        system(command_out_bam_sorted)
        
        
      } else{
        message("Ya se ha mapeado")
      }
      
      
    }
    
  }


filter1 <- function(input_directory = output_directory) {
  carpeta <- "mapping_output"
  output_folder <- file.path(input_directory, carpeta)
  input_bam <-
    list.files(output_folder, pattern = ".sorted.bam$", full.names = T)
  input_bam.name <-
    list.files(output_folder, pattern = ".sorted.bam$", full.names = F)
  input_bam.name.sans_ext <- file_path_sans_ext(input_bam.name)
  out_file <-
    file.path(output_folder,
              paste(input_bam.name.sans_ext, "filtraje_1.bam", sep = "_"))
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        "bamtools filter -in",
        input_bam,
        "-out",
        out_file,
        " -isMapped true  -isMateMapped true"
      )
    print(command)
    system(command)
  } else{
    message("Ya esta filtrado 1")
  }
  
  
  
  
}


RmDup <- function(input_directory = output_directory) {
  carpeta <- "mapping_output"
  output_folder <- file.path(input_directory, carpeta)
  input_bam_filt <-
    list.files(output_folder, pattern = "filtraje_1.bam$", full.names = T)
  input_bam_filt_sans_ext <- file_path_sans_ext(input_bam_filt)
  
  output_file_sort_bam <-
    paste(input_bam_filt_sans_ext, "sort.bam", sep = "_")
  output_file_sort_bam_sans_ext <-
    file_path_sans_ext(output_file_sort_bam)
  
  output_file_collate_bam <-
    paste(output_file_sort_bam_sans_ext, "collate.bam", sep = "_")
  output_file_collate_bam_sans_ext <-
    file_path_sans_ext(output_file_collate_bam)
  
  
  output_file_fixmate <-
    paste(output_file_collate_bam_sans_ext, "fixmate.bam", sep = "_")
  output_file_fixmate_sans_ext <-
    file_path_sans_ext(output_file_fixmate)
  
  
  output_file_fixmate_sort <-
    paste(output_file_collate_bam_sans_ext, "sort.bam", sep = "_")
  output_file_fixmate_sort_sans_ext <-
    file_path_sans_ext(output_file_fixmate_sort)
  
  
  output_file_rmdup <-
    paste(output_file_fixmate_sort_sans_ext, "rmdup.bam", sep = "_")
  output_file_rmdup_sans_ext <-
    file_path_sans_ext(output_file_rmdup)
  
  
  output_file_rmdup_sort <-
    paste(output_file_rmdup_sans_ext, "sorted.bam", sep = "_")
  output_file_rmdup_sort_sans_ext <-
    file_path_sans_ext(output_file_rmdup_sort)
  
  
  if (!file.exists(output_file_sort_bam) |
      !file.exists(output_file_collate_bam) |
      !file.exists(output_file_fixmate) |
      !file.exists(output_file_fixmate_sort) |
      !file.exists(output_file_rmdup) |
      !file.exists(output_file_rmdup_sort)) {
    ##### removing duplicates#####
    print("removing duplicates...")
    print("...Primero sorteamos")
    
    if (!file.exists(output_file_sort_bam)) {
      comando1 <-
        paste("samtools sort -n -o ",
              output_file_sort_bam,
              input_bam_filt)
      
      print(comando1)
      system(comando1)
      
    } else{
      message("Ya esta sorteado")
      print(output_file_collate_bam)
      print(output_file_sort_bam)
    }
    
    if (!file.exists(output_file_collate_bam)) {
      print("... Segundo collate")
      comando2 <-
        paste("samtools collate -o ",
              output_file_collate_bam,
              output_file_sort_bam)
      
      print(comando2)
      system(comando2)
    } else{
      message("Ya esta collatated")
    }
    
    if (!file.exists(output_file_fixmate)) {
      print("...Tercero anyadimos ms and MC tags para usar con markdup")
      comando3 <-
        paste("samtools fixmate -m",
              output_file_collate_bam,
              output_file_fixmate)
      
      print(comando3)
      system(comando3)
    } else{
      message("Ya esta fixeado")
    }
    
    if (!file.exists(output_file_fixmate_sort)) {
      print("... markdup necesita posicion orden")
      comando4 <-
        paste("samtools sort -o ",
              output_file_fixmate_sort,
              output_file_fixmate)
      
      print(comando4)
      system(comando4)
      
      
    } else{
      message("Ya esta fixeado y sorteado")
    }
    if (!file.exists(output_file_rmdup)) {
      print("... Llamamos a markdup")
      comando5 <-
        paste("samtools markdup ",
              output_file_fixmate_sort,
              output_file_rmdup)
      
      print(comando5)
      system(comando5)
      
    } else{
      message("ya esta marcado")
    }
    if (!file.exists(output_file_rmdup_sort)) {
      print("... Sorteamos para que podamos utilizar el archivo luego")
      comando6 <-
        paste("samtools sort -o ",
              output_file_rmdup_sort,
              output_file_rmdup)
      
      print(comando6)
      system(comando6)
      
    } else{
      message("Ya esta marcado y sorteado")
    }
    
    
    
    
  }
  
  else{
    message("Ya se han removido los duplicados")
  }
  
}










freebayes<-function(input_directory=output_directory,reference_genome =folder_fasta){

  
  extension = c(unlist(lapply(list.files(reference_genome, pattern = "fa$"), function(x)
    file_ext(x))),
  unlist(lapply(list.files(reference_genome, pattern = "fasta$"), function(x)
    file_ext(x))))
  print(extension)
  if (extension== "fasta" | extension == "fa") {
    extension <- extension
  } else{
    errorCondition("No existe archivo fasta")
  }
  fasta_file <-
    list.files(reference_genome,
               pattern = paste0(".", extension,"$"),
               full.names = T)
  print(fasta_file)


  output_folder<- file.path(input_directory,"calling_variants")
  input_folder_bams <-file.path(input_directory,"mapping_output")
  output_file.vcf <- file.path(output_folder,"freebayes.vcf")
  rmdup_file <- list.files(input_folder_bams,pattern="rmdup_sorted.bam$",full.names = T)



  if(!dir.exists(output_folder)){
    dir.create(output_folder)
    comando.index <- paste("samtools index ",rmdup_file)

    print(comando.index)
    system(comando.index)

    comando <- paste("freebayes -f",fasta_file,rmdup_file,">",output_file.vcf )
    print("LLAMANDO A LAS VARIANTES")
    print(comando)
    system(comando)





  }else if(dir.exists(output_folder) && !file.exists(output_file.vcf)){

    comando.index <- paste("samtools index ",rmdup_file)

    print(comando.index)
    system(comando.index)

    comando <- paste("freebayes -f",fasta_file,rmdup_file,">",output_file.vcf )
    print("LLAMANDO A LAS VARIANTES")
    print(comando)
    system(comando)




  }else if(dir.exists(output_folder) && file.exists(output_file.vcf)){

    message("Ya se han llamado a las variantes")
    }




}


freebayes()








