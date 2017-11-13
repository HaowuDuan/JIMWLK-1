CC=g++
INCLUDES = -I /sw/include -I/sw/opt/boost-1_58/include/ -I/home/vskokov/lib/include
LIBS= -L /sw/lib/ -L/home/vskokov/lib/lib 
CFLAGS= -std=c++11 -O2 -lgsl -lgslcblas -lfftw3


all: 
	g++ -O3 -o jimwlk.x JIMWLK_TMD_beta_LATEST_dipole_probability.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3
MV: 
	g++ -O3 -o mv.x JIMWLK_TMD_beta_LATEST_dipole_probability_mv.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3
Qs:
	g++ -O3 -o jimwlk_qs.x JIMWLK_Qs.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3

SIP:
	g++ -O3 -o sip.x JIMWLK_SIP.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3

SIPq:
	g++  -I /home/vskokov/lib//include  -I /sw/include -I/sw/opt/boost-1_58/include/ -L /sw/lib/ -L  /home/vskokov/MC_DiJet/interp2d -L /home/vskokov/lib/lib  -std=c++11  -O3 -o sip_q.x JIMWLK_SIP_quark.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3


SIPrc:
	g++  -I/usr/local/include -L/usr/local/lib  -I /home/vskokov/lib//include  -I /sw/include -I/sw/opt/boost-1_58/include/ -L /sw/lib/ -L  /home/vskokov/MC_DiJet/interp2d -L /home/vskokov/lib/lib  -std=c++11  -O3 -o sip_rc.x JIMWLK_rc_SIP.cpp zheevc3.c zheevv3.c  -lgsl -lgslcblas -lfftw3

clean:
	rm -rf jimwlk.x
