run: main.cpp piping.o configuration.o segment.o construct.o
	g++ -std=c++11 -I./matrix_stru3 -o run main.cpp piping.o configuration.o segment.o construct.o

construct.o: construct.cpp construct.h
	g++ -std=c++11 -I./matrix_stru3 -c construct.cpp

piping.o: piping.cpp piping.h
	g++ -c piping.cpp

configuration.o: configuration.cpp configuration.h
	g++ -c configuration.cpp

segment.o: segment.cpp segment.h
	g++ -std=c++11  -I./matrix_stru3 -c segment.cpp

clear:
	rm *.o construct
