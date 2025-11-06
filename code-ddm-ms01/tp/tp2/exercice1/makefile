GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../../../eigen/ -I../../../femtool/


ex: ex.o
	$(GCC) $(FLAGS) ex.o -o ex

ex.o: ex.cpp
	$(GCC) $(FLAGS) $(INC) -c ex.cpp -o ex.o

clean:
	rm *.o
