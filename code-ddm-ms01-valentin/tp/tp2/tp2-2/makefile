GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../../eigen/ -I../../../femtool/

tp2-2: tp2-2.o
	$(GCC) $(FLAGS) tp2-2.o -o tp2-2

tp2-2.o: tp2-2.cpp
	$(GCC) $(FLAGS) $(INC) -c tp2-2.cpp -o tp2-2.o

clean:
	rm *.o