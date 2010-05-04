include ../../../include.mk
binPath = ../../../bin
libPath = ../../../lib

py_progs = $(notdir $(wildcard *.py))

all: ${py_progs:%=${binPath}/%}

${binPath}/%: %
	@mkdir -p $(dir $@)
	cp -f $< $@
	chmod 775 $@

clean:
	rm -f ${py_progs:%=${binPath}/%} ${binPath}/simCtrl_*
