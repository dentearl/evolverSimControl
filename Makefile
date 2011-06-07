SHELL:=/bin/bash -e
export SHELLOPTS=pipefail

.PHONY: all clean

binPath = bin
libPath = lib
py_progs = simCtrl_runSim.py simCtrl_postSimMafExtractor.py simCtrl_postSimFastaExtractor.py 
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
	rm -rf ${binPath} ${libPath}/
