GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../../eigen/ -I../../../femtool/

tp2-1: tp2-1.o
	$(GCC) $(FLAGS) tp2-1.o -o tp2-1

tp2-1.o: tp2-1.cpp
	$(GCC) $(FLAGS) $(INC) -c tp2-1.cpp -o tp2-1.o

clean:
	rm *.o