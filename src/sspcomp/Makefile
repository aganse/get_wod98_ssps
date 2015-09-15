# Makefile to compile sspcm2 sndspd function and little ssp program

CC=gcc

sspcomp: sspcomp.o sspcm2.o Makefile
	gcc -O -pedantic -o sspcomp sspcomp.o sspcm2.o -lm

clean:
	\rm *.o sspcomp
