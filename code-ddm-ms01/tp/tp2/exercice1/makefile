GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../../eigen/ -I../../../femtool/


tp2-exo1: tp2-exo1.o
	$(GCC) $(FLAGS) tp2-exo1.o -o tp2-exo1

tp2-exo1.o: tp2-exo1.cpp
	$(GCC) $(FLAGS) $(INC) -c tp2-exo1.cpp -o tp2-exo1.o

clean:
	rm *.o
