SHELL:=/bin/bash -e
export SHELLOPTS=pipefail

binPath = bin
py_progs = $(notdir $(wildcard src/*.py))

all: ${py_progs:%=${binPath}/%}

${binPath}/%: src/%
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

clean:
	rm -f ${py_progs:%=${binPath}/%} ${binPath}/simCtrl_*
