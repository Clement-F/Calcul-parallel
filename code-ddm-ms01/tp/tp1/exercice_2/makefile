GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../../eigen/ -I../../../femtool/


tp1-exo2: tp1-exo2.o
	$(GCC) $(FLAGS) tp1-exo2.o -o tp1-exo2

tp1-exo2.o: tp1-exo2.cpp
	$(GCC) $(FLAGS) $(INC) -c tp1-exo2.cpp -o tp1-exo2.o

clean:
	rm *.o
