#!/bin/bash
# Copyright (C) 2008-2011 by
# George Asimenos, Robert C. Edgar, Serafim Batzoglou and Arend Sidow.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
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
