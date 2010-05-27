#!/bin/bash
# cycleStandAlone.sh
#  a generic shell wrapper for an evolver cycle
#  dent earl, dearl a soe ucsc edu
#
# $ ./cycle.sh root/ child params/ 0.002
# 
# the last value, the cycle step size (in branch length)
# is optional, the default value is 0.001.
#
# October 2009 - 0.3
########################################
usage() {
    echo "USAGE: $0 parentDir/ childDir[ must not exist ] globalParamsDir/"
    exit 2
}
if [ ! -n "$1" ];then
    usage
fi
if [ -n "$5" ];then
    usage
fi
PAR_DIR=$1                               # existing parent (aka ancestor) directory.
CHI_DIR=$2                               # child (descendant) dircetory, to be created.
GPARAMS_DIR=$3                           # global parameters directory.
STEP_SIZE='0.001'
if [ -n "$4" ];then
    STEP_SIZE=$4
fi
set -beEu -o pipefail                    # from mark d !
evo=/cluster/home/dearl/sonTrace/src/eval/evolver/bin/evo64
cvt=/cluster/home/dearl/sonTrace/src/eval/evolver/bin/cvt
transalign=/cluster/home/dearl/sonTrace/src/eval/evolver/bin/transalign
PAR_DIR=${PAR_DIR%/}
CHI_DIR=${CHI_DIR%/}
PAR=$(echo $PAR_DIR | tr '\/' '\n' | tail -1 )
CHI=$(echo $CHI_DIR | tr '\/' '\n' | tail -1 )

if [ ! -d "$PAR_DIR" ];then
    echo "Error, parent, $PAR_DIR is not a directory"
    usage
fi
if [ -d "$CHI_DIR" ];then
    echo "Error, child, $CHI_DIR is not a directory"
    usage
fi
if [ ! -d "$GPARAMS_DIR" ];then
    echo "Error, params, $GPARAMS_DIR is not a directory"
    usage
fi

########################################
# step one, create child directory and internal
# structure
mkdir $CHI_DIR
mkdir $CHI_DIR/logs
mkdir $CHI_DIR/stats
mkdir $CHI_DIR/chr
mkdir $CHI_DIR/inter

########################################
# Inter step 
########################################
$evo -interchr $PAR_DIR/seq.rev                       \
    -inannots $PAR_DIR/annots.gff                     \
    -aln $CHI_DIR/inter/inter.aln.rev                 \
    -outchrnames $CHI_DIR/inter/inter.chrnames.txt    \
    -outannots $CHI_DIR/inter/inter.outannots.gff     \
    -outseq $CHI_DIR/inter/inter.outseq.rev           \
    -outgenome $CHI.inter                             \
    -branchlength $STEP_SIZE                          \
    -statsfile $CHI_DIR/stats/inter.stats             \
    -model $GPARAMS_DIR/model.txt -seed 1 -logevents  \
    -log $CHI_DIR/logs/inter.log
if [ "$?" -gt 0 ];then
    echo 'Error in first inter step.'
    exit
fi
########################################
# Read out all chromosome names 
########################################
N=0
for chrName in ` cat $CHI_DIR/inter/inter.chrnames.txt `;do
    CHROMS[$N]=$chrName
    ((N++))
done
########################################
# Transalign 1
########################################
if [ -f $PAR_DIR/root.aln.rev ];then
    $transalign -in1 $CHI_DIR/inter/inter.aln.rev   \
        -in2 $PAR_DIR/root.aln.rev                  \
        -out $CHI_DIR/inter/$PAR.inter.aln.rev      \
        -log $CHI_DIR/logs/transalign.inter.log &
else 
    # this is the base case, when the parent *is* the root
    ln -s inter.aln.rev $CHI_DIR/inter/$PAR.inter.aln.rev
fi
if [ "$?" -gt 0 ];then
    echo 'Error in first transalign step.'
    exit
fi
########################################
# Intra steps PARALLEL 1
########################################
i=0
for chr in ${CHROMS[@]}
do
  $evo -inseq $CHI_DIR/inter/inter.outseq.rev       \
      -chrname $chr                                 \
      -seed 1                                       \
      -branchlength $STEP_SIZE                      \
      -mes $GPARAMS_DIR/mes.fa	                    \
      -inannots $CHI_DIR/inter/inter.outannots.gff  \
      -statsfile $CHI_DIR/stats/$chr.stats          \
      -codonsubs $CHI_DIR/chr/$chr.codonsubs	    \
      -outannots $CHI_DIR/chr/$chr.outannots.gff    \
      -outgenome $CHI.intra                         \
      -model $GPARAMS_DIR/model.txt                 \
      -aln $CHI_DIR/chr/$chr.aln.rev	            \
      -outseq $CHI_DIR/chr/$chr.outseq.rev          \
      -log $CHI_DIR/logs/evo.$chr.log & 
  JOBS[$i]=$!
done
while [ "$i" -gt 0 ];do
    wait ${JOBS[$i]}
    if [ "$?" -gt 0]; then
        echo 'Error in intra chr parallel step.'
        exit
    fi
    ((i--))
done
wait
########################################
# Merge steps PARALLEL 2
########################################
i=0
for chr in ${CHROMS[@]}
do # this is just making three commands that include all chromosomes.
  if [ "$i" -eq 0 ]; then
      CAT_COMMAND="$CHI_DIR/chr/$chr.outannots.gff"
      EVO_MERGE="$CHI_DIR/chr/$chr.aln.rev"
      CVT_MERGE="$CHI_DIR/chr/$chr.outseq.rev"
      ((i++)) # ignore first element, it's done above
  else
      CAT_COMMAND="$CAT_COMMAND $CHI_DIR/chr/$chr.outannots.gff"
      EVO_MERGE="$EVO_MERGE,$CHI_DIR/chr/$chr.aln.rev"
      CVT_MERGE="$CVT_MERGE,$CHI_DIR/chr/$chr.outseq.rev"
  fi
done
cat $CAT_COMMAND > $CHI_DIR/annots.gff & 
$evo -mergechrs $EVO_MERGE                 \
    -outgenome $CHI.intra                  \
    -out $CHI_DIR/chr/intra.aln.rev &
JOBS[0]=$!
$cvt -mergerevseqs $CVT_MERGE -out $CHI_DIR/seq.rev &
JOBS[1]=$!
wait ${JOBS[0]}
if [ "$?" -gt 0 ]; then
    echo 'Error in evo merge chr algns step.'
    exit
fi
wait ${JOBS[1]}
if [ "$?" -gt 0 ]; then
    echo 'Error in cvt merge chr seqs  step.'
    exit
fi
########################################
# Transalign 2, CDSalign PARALLEL 3
########################################
$transalign -in1 $CHI_DIR/inter/$PAR.inter.aln.rev \
    -in2 $CHI_DIR/chr/intra.aln.rev                \
    -out $CHI_DIR/root.aln.rev                     \
    -log $CHI_DIR/logs/transalign.log &
JOBS[0]=$!
$evo -cdsalns $CHI_DIR/chr/intra.aln.rev        \
    -alns $CHI_DIR/cdsalns.rev                  \
    -annots1 $CHI_DIR/inter/inter.outannots.gff \
    -annots2 $CHI_DIR/annots.gff                \
    -outgenome $CHI.final                       \
    -log $CHI_DIR/logs/cdsalns.log &
JOBS[1]=$!
wait ${JOBS[0]}
if [ "$?" -gt 0 ]; then
    echo 'Error in final transalign step.'
    exit
fi
wait ${JOBS[1]}
if [ "$?" -gt 0 ]; then
    echo 'Error in cdsalns step.'
    exit
fi

########################################
exit 0
