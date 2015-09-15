# Makefile to compile oclfilt, with or without GMT database-bathy functionality

CC = gcc
CFLAGS = -O -pedantic -ansi -Wall
LIBS = -lm

oclfilt: oclfilt.c getOCLStationData.c ocl.h
	${CC} ${CFLAGS} -o oclfilt oclfilt.c getOCLStationData.c ${LIBS}

outputAllLatsLons: outputAllLatsLons.c getOCLStationData.c ocl.h
	${CC} ${CFLAGS} -o outputAllLatsLons outputAllLatsLons.c \
	getOCLStationData.c ${LIBS}

clean:
	# deleting object files and temp files
	\rm -f *.o *.~ oclfilt

