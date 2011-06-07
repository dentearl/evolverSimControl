SHELL:=/bin/bash -e
export SHELLOPTS=pipefail

.PHONY: all clean

binPath = bin
libPath = lib
py_progs = simCtrl_runSim.py simCtrl_postSimMafExtractor.py
libraries = libSimControl.py libSimControlClasses.py

all: ${py_progs:%=${binPath}/%} ${progs} $(foreach l,${libraries}, ${libPath}/$l)

${libPath}/%: src/%
	@mkdir -p $(dir $@)
	touch ${libPath}/__init__.py
	cp -f $< $@.tmp
	mv $@.tmp $@

${binPath}/%: src/%
	@mkdir -p $(dir $@)
	cp -f $< $@.tmp
	chmod 775 $@.tmp
	mv $@.tmp $@

clean:
	rm -f ${py_progs:%=${binPath}/%} ${binPath}/simCtrl_* ${binPath}/simUtil_* ${libPath}/
	rm -f *.o
