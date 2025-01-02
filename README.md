# msc_thesis_model
This is the code for Hanif Kawousi's msc thesis model, depicting the evolution of genes in digital songbirds and ecological experiments where they are subjected to increased field of predation. 



! These are the main code files for Hanif Kawousi's Msc-project.
! To run it, make sure to have two_gene_model.f90 (the model code),
! the Makefile and main.f90 in the same directory.
! Then run "make" in the terminal for running in evolutionary mode, or
! run "make gen_(your evolutionary generation number) to run in ecological mode.
! For the ecological mode, make sure to have the csv file produced by the evolutionary mode
! for your wanted generation in the same directory.
! ex: For the thesis, I have run 300 generations, producing 300 csv files.
! For ecological experimentations, in the terminal, I ran the following commands: 
!  "make gen_300" - this compiles the code with the generation number 300.
! If changes are made in the code and you do not want to recompile with all the
! outputs in the terminal, you can run "make model". 
! If you want to clear all produced data from previous compilation, run 
! "make clean". 
