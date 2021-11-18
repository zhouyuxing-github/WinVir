#!/bin/env bash
#################################
# author:	"ZhouYuXing"                                                           
# copyright:	"Copyright 2021, Southwest Minzu University"  
# version:	"2.6"                                                                           
# maintainer: "Zhouyuxing"                                                       
# email: "1037782920@qq.com"                                                
#####################################

######     get bin directory    ######
SOURCE="$0"
while [ -h "$SOURCE"  ]; do # resolve $SOURCE until the file is no longer a symlink
    bin="$( cd -P "$( dirname "$SOURCE"  )" && pwd  )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /*  ]] && SOURCE="$bin/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
bin="$( cd -P "$( dirname "$SOURCE"  )" && pwd  )"
#######         finished         #######


cat $bin/databases/*cut*  $bin/databases/*fa  $bin/databases/*fas  $bin/databases/*fasta> $bin/databases/db
time $bin/programs/blast-2.11.0+/bin/makeblastdb.exe  -dbtype nucl  -in $bin/databases/db -out $bin/databases/db && \
awk '/>/{print $0}' $bin/databases/db > $bin/databases/db_idName && \
flag=0 && \
echo "blastdb was generated at `date`" && \
echo "All done... Remove this script plz..."

if [ "$flag" != 0 ];
then
	echo  -e "\n\nwinVir installation failed. The path of the program must be in English with no blank space (use uderline but not space!)"
	echo  -e "such as \"/d/bio soft/winVir2.0/\" is incorrect; \"/d/bio_soft/winVir2.0/\" is correct"
	echo  -e "And make sure your fasta dataset is in /database folder."
	exit
fi
