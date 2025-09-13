util.o: util.cpp
	g++ -c util.cpp

bnsch.o: bnsch.cpp
	g++ -c bnsch.cpp

bnsch.out: util.o bnsch.o
	g++ -o bnsch.out bnsch.o util.o

file:
	mkdir -p data1
	mkdir -p data2
	mkdir -p data3

clean:
	rm *.o

all: file bnsch.out clean 

delete:
	rm -rf data1 data2 data3
	rm *.out
	find . -maxdepth 1 -name '*.m' ! -name 'show.m' ! -name 'draw_area.m' ! -name 'draw_mass.m' -delete
	rm *.eps
