




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










freebayes <-
  function(input_directory = output_directory,
           reference_genome = folder_fasta) {
    extension = c(unlist(lapply(list.files(reference_genome, pattern = "fa$"), function(x)
      file_ext(x))),
      unlist(lapply(list.files(reference_genome, pattern = "fasta$"), function(x)
        file_ext(x))))
    print(extension)
    if (extension == "fasta" | extension == "fa") {
      extension <- extension
    } else{
      errorCondition("No existe archivo fasta")
    }
    fasta_file <-
      list.files(reference_genome,
                 pattern = paste0(".", extension, "$"),
                 full.names = T)
    
    
    output_folder <- file.path(input_directory, "calling_variants")
    input_folder_bams <-
      file.path(input_directory, "mapping_output")
    output_file.vcf <- file.path(output_folder, "freebayes.vcf")
    rmdup_file <-
      list.files(input_folder_bams,
                 pattern = "rmdup_sorted.bam$",
                 full.names = T)
    
    
    
    if (!dir.exists(output_folder)) {
      dir.create(output_folder)
      comando.index <- paste("samtools index ", rmdup_file)
      
      print(comando.index)
      system(comando.index)
      
      comando <-
        paste("freebayes -f",
              fasta_file,
              rmdup_file,
              ">",
              output_file.vcf)
      print("LLAMANDO A LAS VARIANTES")
      print(comando)
      system(comando)
      
      
      
      
      
    } else if (dir.exists(output_folder) &&
               !file.exists(output_file.vcf)) {
      comando.index <- paste("samtools index ", rmdup_file)
      
      print(comando.index)
      system(comando.index)
      
      comando <-
        paste("freebayes -f",
              fasta_file,
              rmdup_file,
              ">",
              output_file.vcf)
      print("LLAMANDO A LAS VARIANTES")
      print(comando)
      system(comando)
      
      
      
      
    } else if (dir.exists(output_folder) &&
               file.exists(output_file.vcf)) {
      message("Ya se han llamado a las variantes")
    }
    
    
    
    
  }




bcfnorm <-
  function(input_directory = output_directory,
           reference_genome = folder_fasta) {
    extension = c(unlist(lapply(list.files(reference_genome, pattern = "fa$"), function(x)
      file_ext(x))),
      unlist(lapply(list.files(reference_genome, pattern = "fasta$"), function(x)
        file_ext(x))))
    print(extension)
    if (extension == "fasta" | extension == "fa") {
      extension <- extension
    } else{
      errorCondition("No existe archivo fasta")
    }
    fasta_file <-
      list.files(reference_genome,
                 pattern = paste0(".", extension, "$"),
                 full.names = T)
    
    input_file <-
      file.path(input_directory, "calling_variants", "freebayes.vcf")
    output_file <-
      file.path(input_directory, "calling_variants", "freebayes_norm.vcf")
    
    if (!file.exists(output_file)) {
      command <-
        paste(
          "bcftools norm --fasta",
          fasta_file,
          " --check-ref -w --multiallelics '-both' --site-win 1000 --output-type v ",
          input_file,
          " >",
          output_file
        )
      
      print(command)
      system(command)
    } else{
      message("Ya se ha normalizado")
    }
    
    
  }



snpeff <-
  function(input_directory = output_directory,
           dir_snpeff = dir_snpeff_,
           reference_genome = folder_fasta) {
    output_folder <- file.path(input_directory, "anotacion")
    
    extension = c(unlist(lapply(list.files(reference_genome, pattern = "fa$"), function(x)
      file_ext(x))),
      unlist(lapply(list.files(reference_genome, pattern = "fasta$"), function(x)
        file_ext(x))))
    print(extension)
    if (extension == "fasta" | extension == "fa") {
      extension <- extension
    } else{
      errorCondition("No existe archivo fasta")
    }
    fasta_file <-
      list.files(reference_genome,
                 pattern = paste0(".", extension, "$"),
                 full.names = T)
    fasta_ref <- unlist(strsplit(fasta_file, "/"))
    fasta_ref <- fasta_ref[grep("hg", ignore.case = T, fasta_ref)]
    fasta_ref <- unlist(strsplit(fasta_ref, "[A-Za-z]"))
    fasta_ref <- ifelse(fasta_ref == "", NA, fasta_ref)
    referencia_version <-
      paste0("hg", fasta_ref[complete.cases(fasta_ref)])
    file_input <-
      file.path(input_directory, "calling_variants", "freebayes_norm.vcf")
    
    snpeff_program <-
      list.files(dir_snpeff, "snpEff.jar", full.names = T)
    snpsift_program <-
      list.files(dir_snpeff, "SnpSift.jar", full.names = T)
    if (referencia_version == "hg19") {
      output_folder <-
        file.path(input_directory, paste0("anotacion", referencia_version))
      
      file_output1.hg19 <-
        file.path(output_folder, "annotated_hg19.vcf")
      file_output2.hg19 <-
        file.path(output_folder, "annotated_hg19_pred.vcf")
      file_output3.hg19 <-
        file.path(output_folder, "annotated__hg19_clin.vcf")
      
      campos_data_hg19 <-
        list.files(file.path(dir_snpeff, "data", "dbNSFP", "hg19"), pattern = "gz$",full.names = T)
      clinvar_hg19 <-
        list.files(file.path(dir_snpeff, "data", "clinvar", "GRCh37"), "gz$",full.names = T)
      if (!dir.exists(output_folder)) {
        dir.create(output_folder)
        if (!file.exists(file_output1.hg19) |
            !file.exists(file_output2.hg19) |
            !file.exists(file_output3.hg19)) {
          
          if(!file.exists(file_output1.hg19)){
            comando1_hg19 <-
              paste(
                "java -Xmx4g -jar ",
                snpeff_program,
                referencia_version,
                "-v",
                file_input,
                ">",
                file_output1.hg19
              )
            
            print(comando1_hg19)
            
            system(comando1_hg19)
            
            
          }else{
            message("Ya se ha andotado 1 version hg19")
          }

          if(!file.exists(file_output2.hg19)){
            comando2_hg19 <-
              paste(
                "java -Xmx4g -jar ",
                snpsift_program,
                "DbNsfp -db",
                campos_data_hg19,
                "-f Uniprot_acc,Interpro_domain,1000Gp3_AC,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,SIFT_pred,SIFT_score,Polyphen2_HDIV_pred,Polyphen2_HDIV_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score,MutationTaster_pred,MutationTaster_score,phastCons100way_vertebrate,phastCons100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP30way_mammalian_rankscore,phyloP17way_primate_rankscore,integrated_fitCons_score -v",
                file_output1.hg19,
                ">",
                file_output2.hg19
                
                
              )
            
            print(comando2_hg19)
            system(comando2_hg19)
            
            
          }else{
            message("Ya se ha anotado 2 version hg19")
          }
          if(!file.exists(file_output3.hg19)){
            comando3_clinvar.hg19 <-
              paste("java -Xmx4g -jar",
                    snpsift_program,
                    "annotate",
                    file_output2.hg19,
                    clinvar_hg19,
                    "-v >",
                    file_output3.hg19)
            print(comando3_clinvar.hg19)
            
            system(comando3_clinvar.hg19)
            
            
            
          }else{
            message("Ya se ha anotado 3 version hg19")
          }

          
        } else{
          message("Ya se ha anotado el exoma")
        }
        
      } else if (dir.exists(output_folder)) {
        if (!file.exists(file_output1.hg19) |
            !file.exists(file_output2.hg19) |
            !file.exists(file_output3.hg19)) {
          
          if(!file.exists(file_output1.hg19)){
            
            comando1_hg19 <-
              paste(
                "java -Xmx4g -jar ",
                snpeff_program,
                referencia_version,
                "-v",
                file_input,
                ">",
                file_output1.hg19
              )
            
            print(comando1_hg19)
            
            system(comando1_hg19)
            
          }else{
            message("Ya se ha anotado y el directorio existe 1 version hg19")
          }
          if(!file.exists(file_output2.hg19)){
            
            comando2_hg19 <-
              paste(
                "java -Xmx4g -jar ",
                snpsift_program,
                "DbNsfp -db",
                campos_data_hg19,
                "-f Uniprot_acc,Interpro_domain,1000Gp3_AC,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,SIFT_pred,SIFT_score,Polyphen2_HDIV_pred,Polyphen2_HDIV_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score,MutationTaster_pred,MutationTaster_score,phastCons100way_vertebrate,phastCons100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP30way_mammalian_rankscore,phyloP17way_primate_rankscore,integrated_fitCons_score -v",
                file_output1.hg19,
                ">",
                file_output2.hg19
                
                
              )
            
            print(comando2_hg19)
            system(comando2_hg19)
            
          }else{
            message("Ya se ha anotado y el directorio existe 2 version hg19")
          }
          if(!file.exists(file_output3.hg19)){
            
            
            comando3_clinvar.hg19 <-
              paste("java -Xmx4g -jar",
                    snpsift_program,
                    "annotate",
                    file_output2.hg19,
                    clinvar_hg19,
                    "-v >",
                    file_output3.hg19)
            print(comando3_clinvar.hg19)
            
            system(comando3_clinvar.hg19)
            
            
            
          }else{
            message("Ya se ha anotado todo el exoma 3 version hg19")
          }



          
          
          
        } else{
          message("Ya sea ha anotado el exoma
")
        }
      }
      
    } else if (referencia_version == "hg38") {
      output_folder <-
        file.path(input_directory, paste0("anotacion", referencia_version))
      
      file_output1.hg38 <-
        file.path(output_folder, "annotated_hg38.vcf")
      file_output2.hg38 <-
        file.path(output_folder, "annotated_hg38_pred.vcf")
      file_output3.hg38 <-
        file.path(output_folder, "annotated__hg38_clin.vcf")
      
      
      campos_data_hg38 <-
        list.files(file.path(dir_snpeff, "data", "dbNSFP", "hg38"), pattern = "gz$",full.names = T)
      
      
      clinvar_hg38 <-
        list.files(file.path(dir_snpeff, "data", "clinvar", "GRCh37"), "gz$",full.names = T)
      
      if (!dir.exists(output_folder)) {
        dir.create(output_folder)
        if (!file.exists(file_output1.hg38) |
            !file.exists(file_output2.hg38) |
            !file.exists(file_output3.hg38)) {
          
          
          if(!file.exists(file_output1.hg38)){
            
            comando1_hg38 <-
              paste(
                "java -Xmx4g -jar ",
                snpeff_program,
                referencia_version,
                "-v",
                file_input,
                " >",
                file_output1.hg38
              )
            print(comando1_hg38)
            system(comando1_hg38)
            
          }else{
            message("Ya se ha anotado el exoma 1 version hg38")
          }
          if(!file.exists(file_output2.hg38)){
            
            
            comando2_hg38 <-
              paste(
                "java -Xmx4g -jar ",
                snpsift_program,
                "DbNsfp -db",
                campos_data_hg38,
                "-f Uniprot_acc,Interpro_domain,1000Gp3_AC,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,SIFT_pred,SIFT_score,Polyphen2_HDIV_pred,Polyphen2_HDIV_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score,MutationTaster_pred,MutationTaster_score,phastCons100way_vertebrate,phastCons100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP30way_mammalian_rankscore,phyloP17way_primate_rankscore,integrated_fitCons_score -v",
                file_output1.hg38,
                ">",
                file_output2.hg38
              )
            print(comando2_hg38)
            system(comando2_hg38)
          }else{
            message("Ya se ha anotado el exoma 2 version hg38")
          }
          if(!file.exists(file_output3.hg38)){
            
            
            comando3_clinvar.hg38 <-
              paste("java -Xmx4g -jar",
                    snpsift_program,
                    "annotate",
                    clinvar_hg38,
                    file_output2.hg38,
                    "-v >",file_output3.hg38)
            
            print(comando3_clinvar.hg38)
            system(comando3_clinvar.hg38)
            
          }else{
            message("Ya se ha anotado todo el exoma 3 version hg38")
          }

      
        } else if (dir.exists(output_folder)) {
          if (!file.exists(file_output1.hg19) |
              !file.exists(file_output2.hg19) |
              !file.exists(file_output3.hg19)) {
            
            if(!file.exists(file_output1.hg38)){
              
              comando1_hg38 <-
                paste(
                  "java -Xmx4g -jar ",
                  snpeff_program,
                  referencia_version,
                  "-v",
                  file_input,
                  ">",
                  file_output1.hg38
                )
              print(comando1_hg38)
              system(comando1_hg38)
            }else{
              message("Ya se ha anotado y el dorectorio esta creado 1 version hg38")
            }
            
            if(!file.exists(file_output2.hg38)){
              
              comando2_hg38 <-
                paste(
                  "java -Xmx4g -jar ",
                  snpsift_program,
                  "DbNsfp -db",
                  campos_data_hg38,
                  "-f Uniprot_acc,Interpro_domain,1000Gp3_AC,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,SIFT_pred,SIFT_score,Polyphen2_HDIV_pred,Polyphen2_HDIV_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_score,MutationTaster_pred,MutationTaster_score,phastCons100way_vertebrate,phastCons100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP30way_mammalian_rankscore,phyloP17way_primate_rankscore,integrated_fitCons_score -v",
                  file_output1.hg38,
                  ">",
                  file_output2.hg38
                )
              print(comando2_hg38)
              system(comando2_hg38)
            }else{
              message("Ya se ha anotado y el directorio esta creado 2 version hg38")
            }
            if(!file.exists(file_output3.hg38)){
              
              
              comando3_clinvar.hg38 <-
                paste("java -Xmx4g -jar",
                      snpsift_program,
                      "annotate",
                      clinvar_hg38,
                      file_output2.hg38,
                      ">",
                      file_output3.hg38)
              
              print(comando3_clinvar.hg38)
              system(comando3_clinvar.hg38)
              
              
              
            }else{
              
              message("Ya se ha anotado todo el exoma  3 el directorio creado version hg38")
            }

            


            
          } else{
            message("Ya sea ha anotado el exoma
")
          }
        }
        
      }
    } else{
      errorCondition("NO se tomo bien la referencia")
    }
    
    
    
    
    
    
    
    
  }




pipeline <- "../pipeline"
muestra <- "muestraToy"
folder_fq <- "fastq_files"
referencia <- "19" #19 o 38
cromosoma_genoma <- "cromosoma8"
folder_fasta <-
  file.path("../datos", paste0("genomaHg", referencia), cromosoma_genoma)
dir_snpeff_ <- "~/tools/exomas_tools/snpEff/"
files_folder <- file.path(pipeline, muestra, folder_fq)
output_directory <- file.path(pipeline, muestra, paste0("output_dir_",referencia))

fastqc_R()
#
bwa_mem2()
#
filter1()
#
RmDup()
#
freebayes()
#
bcfnorm()

snpeff()