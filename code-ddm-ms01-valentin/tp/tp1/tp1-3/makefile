GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../../eigen/ -I../../../femtool/

tp1-3: tp1-3.o
	$(GCC) $(FLAGS) tp1-3.o -o tp1-3

tp1-3.o: tp1-3.cpp
	$(GCC) $(FLAGS) $(INC) -c tp1-3.cpp -o tp1-3.o

clean:
	rm *.o

cleantxt:
	rm *.txt