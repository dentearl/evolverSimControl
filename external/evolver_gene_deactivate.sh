#!/bin/bash
# Copyright (C) 2008-2011 by
# George Asimenos, Robert C. Edgar, Serafim Batzoglou and Arend Sidow.
# 
# All rights reserved. Reproduced and distributed here with permission.
# 
##############################

if (( $# != 4 ))
then
	echo Usage: $0 parent_gff intermediate_gff output_gff evo_bin
	exit 1
fi

parent_gff="$1"
med_gff="$2"
o_gff="$3"
evo="$4"

old=`cat "$parent_gff" | sed -n 's/^\([^\t]\+\).\+\(gene_index [0-9]\+\).\+/\1\t\2/p' | uniq | sort | uniq | wc -l`
new=`cat "$med_gff" | sed -n 's/^\([^\t]\+\).\+\(gene_index [0-9]\+\).\+/\1\t\2/p' | uniq | sort | uniq | wc -l`
let diff=new-old
echo "Changes in number of genes: $old -> $new"
if (( diff > 0 ))
then
	echo "Deleting $diff genes..."
	$evo -delrandgenes "$med_gff" -n -$diff -out "$o_gff"
else
	echo "No gene deactivation will take place."
	cp "$med_gff" "$o_gff"
fi
