.DEFAULT_GOAL := compile

CC := gcc
CC := ${CC}
CFLAGS := -Wall -O3 -g
INCLUDES := -I/usr/local/include
LIBS := /usr/local/lib/libgnuastro.a /usr/local/lib/libchealpix.so \
		-lgit2 -ltiff -llzma -ljpeg -L/usr/local/lib -lwcs -lcfitsio \
		/usr/lib/x86_64-linux-gnu/libcurl-gnutls.so -lz \
		/usr/local/lib/libgsl.so /usr/local/lib/libgslcblas.so -lm -lc -lpthread -pthread

PROGNAME := testprog
SHELL := /usr/bin/env bash
REGION := 


compile: ${PROGNAME}.c
	if [[ -a healpix-test.fits ]];  \
	then				\
	    rm healpix-test.fits;	\
	fi;				\
	${CC} ${CFLAGS} $< ${INCLUDES} -o ${PROGNAME} ${LIBS} && ./${PROGNAME}

plot:
	python3 scripts/plot.py

ds9: ${REGION}
	ds9 test-pv.fits -zscale -zoom to fit -regions $< 
