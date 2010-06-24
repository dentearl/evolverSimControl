include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags = ${cflags_opt}
cflags += -I../../sonLib/inc
progs = ${binPath}/simUtil_paralogBlockMasker

CFLAGS=${cflags} -std=c99 -pedantic

all: ${progs}

${binPath}/simUtil_paralogBlockMasker: simUtil_paralogBlockMasker.c
	${CC} ${cflags} -I ${libPath} -I ${kentInc} -o $@ simUtil_paralogBlockMasker.c ${libPath}/sonLib.a

clean: 
	rm -f *.o
