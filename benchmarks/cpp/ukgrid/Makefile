CC=g++
CFLAGS= -std=c++11 -Wall -g -pedantic -O3 
LDFLAGS= -std=c++11 -Wall -pedantic -O3 -static 
EXEC=executable
VPATH=../classes

all: $(EXEC)

executable: Case.o Array.o SparseMatrix.o main.o
	$(CC) -o executable Case.o Array.o SparseMatrix.o main.o $(LDFLAGS)

SparseMatrix.o: SparseMatrix.cpp Case.h Array.h
	$(CC) -o SparseMatrix.o -c $< $(CFLAGS)

Array.o: Array.cpp
	$(CC) -o Array.o -c $< $(CFLAGS)

Case.o: Case.cpp
	$(CC) -o Case.o -c $< $(CFLAGS)

main.o: main.cpp Case.h SparseMatrix.h Array.h 
	$(CC) -o main.o -c main.cpp $(CFLAGS)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
