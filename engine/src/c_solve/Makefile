CC = gcc
all: liblp_solve.a scs_solve.o
liblp_solve.a: scs_solve.o
	ar rcs liblp_solve.a scs_solve.o 
scs_solve.o: scs_solve.c
	$(CC) -c scs_solve.c -fPIC -O3 -llpack -lscsdir -lblas -o scs_solve.o
