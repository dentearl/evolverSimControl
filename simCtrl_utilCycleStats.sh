#!/bin/bash
# 16 dec 2009
# dent earl, dearl (a) soe ucsc edu
# prints out useful information about
# a simulation and its reconstruction.
#
##############################

function usage(){
    echo "USAGE: $0 cycleDirectory/"
    exit 2
}
if [ "$#" -ne 1 ]; then
    usage
fi
if [ ! -d $1 ]; then
    usage
fi
DIR=$1
echo -e "Leaf\tSize (bp)\tPost-Mask (bp)\t%"
for a in $(ls $DIR/); do
    if [ -d "$DIR/$a/" ]; then
        if [ -f "$DIR/$a/seq.masked.fa" ]; then
            #MASK=$(cat $DIR/$a/repeatMask/faSize.rmsk.txt | perl -nle 'print $1 if(m/^%(\S+) masked total/); ')
            MASK=$(eval_maskedFaLength.py $DIR/$a/seq.masked.fa)
            SIZE=$(evolver_cvt -showgenomesizes $DIR/$a/seq.rev | awk '{print $1}')
            MASKED=$(echo "scale=8; $MASK/$SIZE" | bc -l)
            echo -e "$a\t$SIZE\t$MASK\t$MASKED"
        fi
    fi
done
RECON=$(eval_MAFLength.py $DIR/reconstruction.maf)
echo "Reconstruction size: $RECON (bp)"
echo -e "\nComparison\tTrue Size\tSensitivity (#)\tSpecificity (#)"
for a in $( ls $DIR*xml ); do
    NAME=$(echo "$a" | perl -na -F/ -e '$a = $F[$#F]; $a =~ s/.aln.clean.maf.compare.xml//; print $a')
    TRUE_SIZE=$( eval_MAFLength.py $DIR/$NAME.aln.maf )
    SENS=$(cat $a | perl -nle 'print $1 if(m/homology_tests fileA="\S+?" fileB="\S+?" totalTests="\S+?" totalTrue="\S+?" totalFalse="\S+?" average="(\S+?)"/); ' | head -1 /dev/stdin )
    SENS_CNT=$(cat $a | perl -nle 'print $1 if(m/homology_tests fileA="\S+?" fileB="\S+?" totalTests="\S+?" totalTrue="(\d+?)\.0+" totalFalse="\S+?" average="\S+?"/); ' | head -1 /dev/stdin )
    SPEC=$(cat $a | perl -nle 'print $1 if(m/homology_tests fileA="\S+?" fileB="\S+?" totalTests="\S+?" totalTrue="\S+?" totalFalse="\S+?" average="(\S+?)"/); ' | tail -1 /dev/stdin )
    SPEC_CNT=$(cat $a | perl -nle 'print $1 if(m/homology_tests fileA="\S+?" fileB="\S+?" totalTests="\S+?" totalTrue="(\d+?)\.0+" totalFalse="\S+?" average="\S+?"/); ' | tail -1 /dev/stdin )
    echo -e "$NAME\t$TRUE_SIZE\t$SENS ($SENS_CNT)\t$SPEC ($SPEC_CNT)"
done
