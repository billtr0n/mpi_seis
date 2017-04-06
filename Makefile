DBG = -Wall

all:
	make param
	make awpdata
	make xapiir
	make peak

param:
	mpicc -c $(DBG) fd3dparam.c

awpdata:
	mpicc -c $(DBG) awp_data.c

peak:
	mpicc $(DBG) -o pgv2 pgv2.c fd3dparam.o awp_data.o xapiir.o

xapiir:
	gfortran -c -O1 xapiir.f

extract:
	mpicc $(DBG) -o extrts extrts.c awp_data.o

clean:
	rm *.o pgv2
