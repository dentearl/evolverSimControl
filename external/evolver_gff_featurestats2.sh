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

set -e

usage(){
    
    echo Usage: $0 gff_featurestats2.py genome1.gff genome2.gff genome1_name genome2_name genome1.seq.rev genome2.seq.rev path_to_cvt >&2
	 echo >&2
	 echo This will run gff_featurestats2.py with the gff files and genome names, >&2
	 echo and with the correct genome sizes extracted from the seq.rev files using cvt. >&2
	 exit 2
}

if (( $# != 8 ))
then
    echo "Not enough variables!" >&2
    usage
fi

gff_featurestats="$1"
g1_gff="$2"
g2_gff="$3"
g1_name="$4"
g2_name="$5"
g1_seq="$6"
g2_seq="$7"
cvt="$8"

if [ ! -f $gff_featurestats ];then
    echo "ERROR, gff_featurestats does not exist at $gff_featurestats" >&2
    usage
fi

if [ ! -f $g1_gff ];then
    echo "ERROR, genome 1 gff does not exist at $g1_gff" >&2
    usage
fi

if [ ! -f $g2_gff ];then
    echo "ERROR, genome 2 gff does not exist at $g2_gff" >&2
    usage
fi

if [ ! -f $g1_seq ];then
    echo "ERROR, genome 1 seq does not exist at $g1_seq" >&2
    usage
fi

if [ ! -f $g2_seq ];then
    echo "ERROR, genome 2 seq does not exist at $g2_seq" >&2
    usage
fi

if [ ! -f $cvt ];then
    echo "ERROR, cvt does not exist at $cvt" >&2
    usage
fi

size1=$($cvt -showgenomesizes $g1_seq 2>/dev/null | awk '{print $1}')
size2=$($cvt -showgenomesizes $g2_seq 2>/dev/null | awk '{print $1}')

"$gff_featurestats" "$g1_gff" "$g2_gff" "$g1_name" "$g2_name" "$size1" "$size2"
