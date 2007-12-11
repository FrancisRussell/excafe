all: main
main: main.cpp
	g++ main.cpp -Wall -pedantic -I ./include -o main
clean:
	rm -f *.o main
