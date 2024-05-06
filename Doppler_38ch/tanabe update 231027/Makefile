CC = gcc -g

all:	gs

gs	:../common/flux.o ../common/gscore_ST.o ../common/gs.o ../common/celi.o ../common/spline.o
	$(CC) -o gs ../common/flux.o ../common/gscore_ST.o ../common/gs.o ../common/celi.o ../common/spline.o -lm

flux.o	:../common/flux.c
	$(CC) -c ../common/flux.c -o ../common/flux.o

gscore_ST.o	:../common/gscore_ST.c
	$(CC) -c ../common/gscore_ST.c -o ../common/gscore_ST.o

gs.o	:../common/gs.c
	$(CC) -c ../common/gs.c -o ../common/gs.o

celi.o :../common/celi.c
	$(CC) -c ../common/celi.c -o ../common/celi.o

spline.o :../common/spline.c
	$(CC) -c ../common/spline.c -o ../common/spline.o

clean:
	(rm -f gs.exe ../common/*.o)
