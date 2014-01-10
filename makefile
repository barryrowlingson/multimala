RANLIB = /home/sut3/AEGISS/LGCP/RANLIB
RANOBJECTS = $(RANLIB)/ranlib.o $(RANLIB)/com.o $(RANLIB)/linpack.o 
LDFLAGS = -L$(RANLIB)/  -lm 
CFLAGS = -O -I$(RANLIB)/ #this means if can't find things, find this directory
CC = gcc 
#FC = gf

multimala: multimala.o
	$(CC) -o multimala multimala.o -lm -O

lilleg:
	R CMD SHLIB -o lilleg2.so lilleg2.c

countvar: countvar.o
	$(CC) -o countvar countvar.c -lm -O
#	R CMD SHLIB -o countvar.so countvar.c -lm

poisson.sim: poisson.sim.c
	$(CC) $(CFLAGS) -o poisson.sim.o  -c poisson.sim.c
	R CMD SHLIB -o r.poisson.sim.so poisson.sim.o $(RANOBJECTS)

