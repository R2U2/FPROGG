# MLTL to SMT Translation
This directory contains the program called by FPROGG to ingest an MLTL formula and return a formula that can be used to query a SAT/SMT solver.

By defualt this is the Fast MLTL to Boolean translator from 
Hariharan et al (2023). https://doi.org/10.1007/978-3-031-42626-1_6.

If you would like to change the underlying engine, simply replace the binary file in this directory.