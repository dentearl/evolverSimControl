SHELL:=/bin/bash -e
export SHELLOPTS=pipefail

binPath = bin
libPath = ../../sonTrace/sonLib/lib
progs = ${binPath}/simUtil_paralogBlockMasker
py_progs = $(notdir $(wildcard src/*.py))
cflags += -I${libPath}/inc -std=c99 -pedantic

all: ${py_progs:%=${binPath}/%} ${progs}

${binPath}/%: src/%
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

${binPath}/simUtil_paralogBlockMasker: src/simUtil_paralogBlockMasker.c
	${CC} ${cflags} -I ${libPath} -o $@.tmp $< ${libPath}/sonLib.a
	mv $@.tmp $@

clean:
	rm -f ${py_progs:%=${binPath}/%} ${binPath}/simCtrl_* ${binPath}/simUtil_* 
	rm -f *.o
