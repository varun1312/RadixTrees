all: clean comp 

debug: clean comp_debug 

opt: clean comp_opt

clean:
	rm -rf a.out

comp:
	g++ nbrt.cpp -pthread -std=c++11 -o nbrt.o

run:
	./nbrt.o 32 70 20 10 50 1000 6543

comp_opt:
	g++ -O3 -pthread -std=c++11 nbrt.cpp -o nbrt.o

run_opt:
	./nbrt.o 32 70 20 10 50 1000 6543 

comp_debug:
	g++ -g -pthread -std=c++11 nbrt.cpp -o nbrt.o

run_debug:
	gdb ./nbrt.o 32 70 20 10 50 1000 6543 
