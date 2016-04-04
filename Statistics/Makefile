##########################################
# Makefile for compiling driver.m        #
#  - 10/30/91  - Started                 #
##########################################


CFLAGS		= -O -Wall -g


MOD		= 
MODC		= 



MAKEFILE	= makefile.
MANPAGES	= 
DOC		= 
MISC		= 
HFILES		= 
MFILES		= driver.m DesStat.m
OFILES		= driver.o DesStat.o
SCRIPTS		=
PROGRAM		= driver



$(PROGRAM): $(OFILES)
	cc $(OFILES) -lm -o $(PROGRAM)
$(OFILES): $(HFILES) 
