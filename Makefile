SHELL:=/bin/bash -e
export SHELLOPTS=pipefail

.PHONY: all clean test

binPath = bin
libPath = lib
extPath = external
py_progs = simCtrl_runSim.py simCtrl_checkSimStatus.py simCtrl_postSimAnnotDistExtractor.py simCtrl_postSimFastaExtractor.py simCtrl_postSimMafExtractor.py
externals= evolver_codon_report.pl evolver_drawrev evolver_gene_deactivate.sh evolver_gff_featurestats2.sh evolver_handle_mobiles.pl evolver_mobile_report.pl
libraries = libSimControl.py libSimControlClasses.py

all: ${py_progs:%=${binPath}/%} $(foreach l,${libraries}, ${libPath}/$l) $(foreach f,${externals},${binPath}/$f)

${libPath}/%: src/% __init__.py
	@mkdir -p $(dir $@)
	touch ${libPath}/__init__.py
	cp -f $< $@.tmp
	mv $@.tmp $@

${binPath}/%: src/%
	@mkdir -p $(dir $@)
	cp -f $< $@.tmp
	chmod 755 $@.tmp
	mv $@.tmp $@

${binPath}/%: external/%
	@mkdir -p $(dir $@)
	cp -f $< $@.tmp
	chmod 755 $@.tmp
	mv $@.tmp $@

__init__.py:
	touch $@

clean:
	rm -rf ${binPath} ${libPath}/ __init__.py*

test:
	python src/simCtrl_testDependencies.py
