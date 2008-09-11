#! /bin/sh
make clean

rm ./aclocal.m4
rm -rf ./autom4te.cache
rm ./config.guess
rm ./config.sub
rm ./config.cache
rm ./config.log
rm ./configure
rm ./depcomp
rm ./INSTALL
rm ./install-sh
rm ./ltmain.sh
rm ./Makefile.in
rm ./missing
rm ./mkinstalldirs
rm ./README
rm ./bppSuite/Makefile.in
rm ./config.status
rm ./libtool
rm ./Makefile
rm ./bppSuite/Makefile
rm -rf ./bppSuite/.deps

rm Examples/*.profile
rm Examples/*.messages
rm Examples/*.mat
rm Examples/warnings
rm Examples/LSU.bionj.dnd
rm Examples/LSU.pars*
rm Examples/LSU.ML*
rm Examples/LSU.asr.csv
rm Examples/LSU.ancestors.*
rm Examples/boot.dnd
rm Examples/treeList2_rooted.dnd
rm Examples/*~
rm Examples/LSUSim.fasta

