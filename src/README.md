# FPROGG 
## Abstract
Mission-time Linear Temporal Logic (MLTL) is a popular finite variant of Linear Temporal Logic with closed interval integer bounds 
    which system designers and engineers use to write specifications that can be monitored with runtime verification (RV). 
    Rigorously testing these runtime monitors is a critical necessity, but designing adequate tests for them is challenging. 
    Therefore, domain experts need an automated method for generating tests at scale. 
    We introduce the open source Formula Progression Generator (FPROGG) tool to construct benchmarks for MLTL specifications through formula progression. 
    FPROGG takes a formula and a desired trace length as input and returns an MLTL benchmark instance. 
    An MLTL benchmark is a three-tuple consisting of a formula; a trace of variable assignments; and an execution sequence of verdict-timestamp tuples. 
    We also contribute a set of pre-computed MLTL benchmarks and their corresponding experimental evaluations. The tool is available at https://temporallogic.org/research/FPROGG
    

## Dependencies
- Flex and Bison
- Z3 must be installed to usr/bin
- Make
- The artifact from [Hariharan et al](10.1007/978-3-031-42626-1_6)
To compile, type "make run" in the terminal from the "Benchmark Generator" directory

## Reproducing Results
- Run `FPROGG_Batch_Testing.sh -S 5` to regenerate the benchmarks. 
fig 4 can then be reproduced by parsing `scripts/total_time.log`


## Supported patterns
There are currently 4 main patterns available to generate benchmarks. Almost-SAT, Almost-UNSAT, Median-SAT, and Semi-Random-SAT. 
	
###	Almost Satisfiable Pattern
Generates a trace that maximizes satisfiability through time via formula progression.
###	Almost Unsatisfiable Pattern
Generates a trace that minimizes satisfiability through time via satisfying the negation of the given formula.
###	Median Satisfiable Pattern
Generates a trace that satisfies roughly half of all timesteps by switching from Almost-Satisfiable to Almost-Unsatisfiable at the halfway mark or vice-versa. Sensitive to formula becoming arbitrarily true or false before the midpoint and being unable to complete trace generation.

###	Semi-Random Satisfible Pattern
Randomly picks time intervals to switch between Almost-Satisfible and Almost Unsatisfiable patterns. Very sensitive to becoming arbitrarily true or false before complete trace generation.
 	 


## Instructions for running
1. Compile the tool by running the 'make run' command in the 'FPROGG' directory.

2. Place your formulas into *.mltl* files. 
	- We support a functional subset of the  standard propostional logic operators  (!, &, |),  and of the temporal operators. ( G[a,b], U[a,b], R[a,b] )
	- All temporal operators have closed interval integer bounds immediately following the operator. 
	- Comment lines in *.mltl* files are denoted with  a # character at the start of a line.
	- Midline comments are not supported.
	- Multiple formulas in a file will generate multiple traces. To combine them, replace all newline characters with the '&' operator. 
	

3. To run FPROGG, execute:
	 
	`./FPROGG [option] [path to mltl file] [length] [optional random seed for random patterns]`
	
	Options: `--help`, `-S`, `-U`, `-M1`, `-M2`, `-SR`
		`-S`: Almost-Satisfiable Pattern
		`-U`: Almost-Unsatisfiable Pattern
		`-M1`: Median-Satisfiable Pattern with Almost-Satisfiable for the first half and Almost-Unsatisfiable for the second half of the generated trace
		`-M2`: Median-Satisfiable Pattern with Almost-Unsatisfiable for the first half and Almost-Satisfiable for the second half of the generated trace
		`-SR`: Semi-Random-Satisfiable Pattern
		`--help`: displays useage help message
		
	Length: desired integer length of generated trace, must be greater than 1.
		
4. Output
	The trace is sent to a *.csv* file in the same directory as the provided *.mltl* file. The name will contain the line number of the formula for *.mltl* files with more than formula.  
	The produced Oracle is sent to a *.txt* file in the same directoty as the generated *.csv* file.
	Both are named [mltl file name][formula index within file]. 
	
	The first formula in a file will have a trace file named *[mltl file name]_0.csv*

## Examples
Consider the formula G[0,1] a0 in 'ex.mltl'

Running the tool under the almost-satisfiable pattern with a length of 4 will yield the trace and Oracle:
	`./FPROGG -S ./ex.mltl 4`	
	
	1	T
	1	T
	1	T
	1
	
Similarly, the same formula under the Almost-Unsatisfiable pattern yields:
	`./FPROGG -U ./ex.mltl 4`

	0	F
	0	F
	0	F
	0
	
The Median-Satisfiabiltiy pattern with the Almost-Satisfiable pattern first yields:
	`./FPROGG -M1 ./ex.mltl 4`

	1	T	
	1	T
	1	F
	0
	
And the Median-Satisfiability pattern with the Almost-Unsatisfiable pattern first: 
	`./FPROGG -M2 ./ex.mltl 4`

	0	F
	0	F
	1	T
	1
	

	
The Semi-Random patten will produce traces similar to the Median patterns.
	
## Misc
FPROGG is partially built off the artifact from [Li et al]()
We use the artifact from [Hariharan et al](10.1007/978-3-031-42626-1_6) for the fast MLTL to propostional logic translation.
We use the Z3 theorem prover for SAT solving.