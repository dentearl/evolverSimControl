<<<<<<< HEAD:Makefile
SHELL:=/bin/bash -e
export SHELLOPTS=pipefail

binPath = bin
py_progs = $(notdir $(wildcard src/*.py))

all: ${py_progs:%=${binPath}/%}

${binPath}/%: src/%
=======
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
	${CC} ${cflags} -I ${libPath} -o $@ simUtil_paralogBlockMasker.c ${libPath}/sonLib.a

${binPath}/%: %
>>>>>>> utilities:Makefile
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

<<<<<<< HEAD:Makefile
clean:
	rm -f ${py_progs:%=${binPath}/%} ${binPath}/simCtrl_*
=======
clean: 
	rm -f *.o
	rm -rf ${py_progs:%=${binPath}/%} ${binPath}/simUtil_*
>>>>>>> utilities:Makefile
