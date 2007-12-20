all: main
main: main.cpp
	g++ main.cpp -Wall -pedantic -O2 -I ./include -I ${HOME}/individual_project/local_mtl/include -o main
clean:
	rm -f *.o main
