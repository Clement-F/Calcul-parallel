GCC=g++ -std=c++23 -Wall
FLAGS=-O3
INC=-I../eigen/ -I../femtool/

tp1-exo1: tp1-exo1.o
	$(GCC) $(FLAGS) tp1-exo1.o -o tp1-exo1

tp1-exo1.o: tp1-exo1.cpp
	$(GCC) $(FLAGS) $(INC) -c tp1-exo1.cpp -o tp1-exo1.o



tp1-exo2: tp1-exo2.o
	$(GCC) $(FLAGS) tp1-exo2.o -o tp1-exo2

tp1-exo2.o: tp1-exo2.cpp
	$(GCC) $(FLAGS) $(INC) -c tp1-exo2.cpp -o tp1-exo2.o



tp2-exo1: tp2-exo1.o
	$(GCC) $(FLAGS) tp2-exo1.o -o tp2-exo1

tp2-exo1.o: tp2-exo1.cpp
	$(GCC) $(FLAGS) $(INC) -c tp2-exo1.cpp -o tp2-exo1.o



tp2-exo2: tp2-exo2.o
	$(GCC) $(FLAGS) tp2-exo2.o -o tp2-exo2

tp2-exo2.o: tp2-exo2.cpp
	$(GCC) $(FLAGS) $(INC) -c tp2-exo2.cpp -o tp2-exo2.o



tp3-exo1: tp3-exo1.o
	$(GCC) $(FLAGS) tp3-exo1.o -o tp3-exo1

tp3-exo1.o: tp3-exo1.cpp
	$(GCC) $(FLAGS) $(INC) -c tp3-exo1.cpp -o tp3-exo1.o



cleantxt:
	rm *.txt

cleanexe:
	rm *.exe

cleanmesh:
	rm *.mesh

cleansol:
	rm *.sol

clean:
	rm *.o