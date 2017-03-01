CC=g++
INCLUDES = -I /sw/include -I/sw/opt/boost-1_58/include/ -I/home/vskokov/lib/include
LIBS= -L /sw/lib/ -L/home/vskokov/lib/lib 
CFLAGS= -std=c++11 -O2 -lgsl -lgslcblas -lfftw3


all: 
	g++ -O3 -o jimwlk.x JIMWLK_TMD_beta_LATEST_dipole_probability.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3
Qs:
	g++ -O3 -o jimwlk_qs.x JIMWLK_Qs.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3
clean:
	rm -rf jimwlk.x
