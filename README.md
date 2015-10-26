# pbwt__shapeit

This is an edited version of PBWT (Positional Burrows-Wheeler Transform) code.
Comparing to original PBWT code, this version changed the pbwtMain.c (Main function), pbwt.h (add some function definition) 
and add one pbwtShapeIt.c (my implementation code).

This version adds the function to shape the genotype.
The input is .pbwt file and output is haplotype info. after shapeit.

Brief usage instructions
------------------------
Typing

    pbwt

by itself gives a list of commands with brief descriptions.

A quick synopsis for usage is:

shapeIt with fixed heterozygous number (3) in each block, output to the file

	./pbwt -read XXX.pbwt -shapeIt1 <file>   

shapeIt with extensible heterozyogous in each block, maxGeno is the max num of heterozyogous in each block,output to the file

	./pbwt -read XXX.pbwt -shapeIt2 maxGeno  <file> 

shapeIt with heterozyogous only probability and original probability, percent 0 represents only original probability, percent 100 represents only heterozyogous only probability. And output to the file
	
	./pbwt -read XXX.pbwt -shapeIt3 percent  <file>    
