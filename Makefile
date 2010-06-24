include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

cflags = ${cflags_opt}
cflags += -I../../sonLib/inc
progs = ${binPath}/simUtil_paralogBlockMasker
py_progs = $(notdir $(wildcard *.py))

CFLAGS=${cflags} -std=c99 -pedantic

all: ${progs} ${py_progs:%=${binPath}/%}

${binPath}/simUtil_paralogBlockMasker: simUtil_paralogBlockMasker.c
	${CC} ${cflags} -I ${libPath} -I ${kentInc} -o $@ simUtil_paralogBlockMasker.c ${libPath}/sonLib.a

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

clean: 
	rm -f *.o
	rm -f ${py_progs:%=${binPath}/%} ${binPath}/simUtil_*
