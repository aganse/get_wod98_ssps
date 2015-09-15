# Top-level makefile to compile oclfilt and sspcomp, for use with get.wod98.ssps

all:
	cd src/oclfilt; make; cp oclfilt ../..; cd ../..
	cd src/sspcomp; make; cp sspcomp ../..; cd ../..

clean:
	cd src/oclfilt; make clean; cd ../..
	cd src/sspcomp; make clean; cd ../..
	\rm oclfilt sspcomp
