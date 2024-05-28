

Greetings, integrons enthusiasts!!

In this repository you find a custom R function to identify integron's endogenous promoters in fasta formatted sequence files. Briefly, this function searchs for exact matchs of the correctly spaced -10 and -35 box of integron promoters described to date (reviewed by [Fonseca and Vicente, 2022](https://doi.org/10.3390/microorganisms10020224)).


# Integron-promoter-finder

## Arguments required for usage

- **file.path:** path to the fasta file (mandatory).

- **integron.class:** string specifying the promoter search for integron class (mandatory). Valid strings are: *"class 1"*, *"class 2"* and *"class 3"*. Only individual searches for each class is currently supported.


## Required packages

- [tidyverse]()
- [seqinr]()
- [spgs]()

library()
  library(tidyverse)
  library(spgs)











