GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../eigen/ -I../../femtool/

tp1-1: tp1-1.o
	$(GCC) $(FLAGS) tp1-1.o -o tp1-1

tp1-1.o: tp1-1.cpp
	$(GCC) $(FLAGS) $(INC) -c tp1-1.cpp -o tp1-1.o

clean:
	rm *.o