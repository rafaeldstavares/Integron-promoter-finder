int.PromoterFINDER <- function(file.path, integron.class){
  cat("Reading file: ", file.path, "\n")
  cat("Identifying promoter sequences of", integron.class, "integrons", "\n")
  # Automatically set the working directory
  setwd(dirname(rstudioapi::getSourceEditorContext()[[2]]))
  # Install packages if not available
  if (!require("tidyverse", quietly = TRUE))
    install.packages("tidyverse")
  if (!require("seqinr", quietly = TRUE))
    install.packages("seqinr")
  if (!require("spgs", quietly = TRUE))
    install.packages("spgs")
  # Loading required packages
  library(seqinr)
  library(tidyverse)
  library(spgs)
  cat(" Dependencies loaded sucessfully", "\n")
  # Importing the data from a fasta file
  test.seq <- seqinr::read.fasta(file = file.path,
                                 as.string = TRUE)
  cat(" Fasta file read sucessfully", "\n")
  name.seq <- names(test.seq)
  sequences <- seqinr::getSequence(test.seq, as.string = TRUE)
  # Importing the table with promoter sequences
  promoter.database <- read_tsv(file = 'promoter.database.tsv')
  promoID <- (promoter.database %>% dplyr::filter(class == integron.class))$promoterID
  cat(" promoter database imported sucessfully", "\n")
  # FORWARD SEQUENCES
  # Creating an empty table
  tab.results.fwd <- tibble()
  # Loop to identify each promoter type in the input sequences
  cat(" Initializaing promoter sequence detection | Forward sequences", "\n")
  for(i in promoID){
    # Get the position of the hits
    cat(" Detecting the promoter ", i, "\n")
    promoter.localization <- str_locate_all(string = sequences, 
                                            pattern = promoter.database$sequence_pattern[promoter.database == i])
    # Attribute the position of the hit to a table with all information
    tab <- tibble(seqID = name.seq,
                  promoterID = i,
                  length = length.seq,
                  start = unlist(lapply(promoter.localization, function(x) x[1]))-nchar(word(str_replace_all(promoter.database$sequence_pattern[promoter.database == i], '\\+', ' '), 1)),
                  end = unlist(lapply(promoter.localization, function(x) x[2]))-nchar(word(str_replace_all(promoter.database$sequence_pattern[promoter.database == i], '\\+', ' '), -1))) %>% 
      dplyr::filter(!is.na(start)) %>% 
      dplyr::select(-c(length, start.x, end.x))
    cat("  Detection step sucessfull", i, "\n")
    # Bind the results together
    tab.results.fwd <- rbind(tab.results.fwd,
                             tab) %>% 
      group_by(seqID, start) %>% 
      dplyr::slice_max(nchar(promoterID)) %>% 
      mutate(promoterCLASS = word(str_replace(promoterID, '\\.', ' '),1)) %>% 
      select(seqID, promoterCLASS, promoterID, start, end) %>% 
      mutate(strand = 1)
    cat("  Hits added to the output table", i, "\n")
  }
  # REVERSE SEQUENCES
  # Creating an empty table
  tab.results.rev <- tibble()
  # Loop to identify each promoter type in the input sequences
  cat(" Initializaing promoter sequence detection | Reverse sequences", "\n")
  for(i in promoID){
    # Get the position of the hits
    cat(" Detecting the promoter ", i, "\n")
    promoter.localization <- str_locate_all(string = spgs::reverseComplement(sequences), 
                                            pattern = promoter.database$sequence_pattern[promoter.database == i])
    # Attribute the position of the hit to a table with all information
    tab <- tibble(seqID = name.seq,
                  promoterID = i,
                  length = length.seq,
                  start = unlist(lapply(promoter.localization, function(x) x[1]))-nchar(word(str_replace_all(promoter.database$sequence_pattern[promoter.database == i], '\\+', ' '), 1)),
                  end = unlist(lapply(promoter.localization, function(x) x[2]))-nchar(word(str_replace_all(promoter.database$sequence_pattern[promoter.database == i], '\\+', ' '), -1))) %>% 
      dplyr::filter(!is.na(start))
    cat("  Detection step sucessfull", i, "\n")
    # Bind the results together
    tab.results.rev <- rbind(tab.results.rev,
                             tab) %>% 
      group_by(seqID, start) %>% 
      dplyr::slice_max(nchar(promoterID)) %>% 
      mutate(promoterCLASS = word(str_replace(promoterID, '\\.', ' '),1)) %>% 
      select(seqID, promoterCLASS, promoterID, start, end) %>% 
      mutate(strand = -1)
    cat("  Hits added to the output table", i, "\n")
  }
  # Return an output table
  tab.results <- rbind(tab.results.fwd, tab.results.rev) %>% 
    arrange(seqID,start, end)
  return(tab.results)
}
