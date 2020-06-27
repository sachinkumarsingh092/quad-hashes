.DEFAULT_GOAL := compile

CC := gcc
CC := ${CC}
CFLAGS := -Wall -O3 -g
INCLUDES := -I/usr/local/include
LIBS := /usr/local/lib/libgnuastro.so.10  /usr/local/lib/libwcs-7.1.a \
	/usr/local/lib/libcfitsio.so.8 /usr/local/lib/libchealpix.so -lm
PROGNAME := testprog
SHELL := /usr/bin/env bash



compile: ${PROGNAME}.c
	if [[ -a healpix-test.fits ]];  \
	then				\
	    rm healpix-test.fits;	\
	fi;				\
	${CC} ${CFLAGS} ${PROGNAME}.c ${INCLUDES} -o ${PROGNAME} -pthread ${LIBS} && ./${PROGNAME}

plot: 
	python3 scripts/plot.py
