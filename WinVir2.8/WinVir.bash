#!/bin/env bash
#################################
# author:	"ZhouYuXing" 
# copyright:	"Copyright 2021, Southwest Minzu University"
# version:	"2.6"
# maintainer: "Zhouyuxing" 
# email: "1037782920@qq.com"
#####################################


while getopts 'i:c:h' arg
do
	case $arg in
		i)
			inFile="$OPTARG";;
		c)
			cutOff="$OPTARG";;
		h)
			echo -e "\nwinVir:\n       Usage: `basename $0`  -i inFile.fastq.gz\n       winVir is a convenient virus annotation tool for metagenomic data on windows"
			exit
	esac
done
flag=0
if [ "$inFile" = "" ];
then
	flag=1
fi
if [ "$cutOff" = "" ];
then
	$cutOff=80
fi
if [ "$flag" = 1 ];
then
	echo -e "\nwinVir:\n       Usage: `basename $0`  -i inFile.fastq.gz"
	exit
fi

##############         get bin, base  and check       ################
bin="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" && \
dirname="$( cd "$( dirname "$inFile" )" && pwd )" && \
basename=`basename  ${inFile%.*}` && \

STR='GNU/Linux is an operating system'
SUB=' '
if [[ "$bin" == *"$SUB"* ]]; then
  echo -e "\n\nWinVir:\n       The path to store the program must be in English with no blank space and illegal characters!"
  echo  -e "       such as \"/d/bio soft/WinVir2.0/\" is incorrect; \"/d/bio_soft/WinVir2.0/\" is correct."
  exit
fi
################################################

#############    check the existence of the database     #############
if [ ! -f $bin/databases/db.ndb ];then
    echo "Start installing winVir..."
	sh $bin/install.bash
else
    echo ""
fi && \
################################################

if [ ! -f $bin/databases/db.ndb ];then #check the existence of the database again in case of  whistling into the wind
    exit
else
    echo "No problem with the database"
fi && \

if [ ! -d $dirname/${basename%.*}_res ];then
    mkdir $dirname/${basename%.*}_res && \
    echo -e "Directory $dirname/${basename%.*}_res was made!\n" && \
    echo "WinVir start at `date`" && \
    echo "Start converting fastq to fasta..." && \
    zcat  $inFile | awk '{getline seq;getline plus;getline qual;sub("@",">",$0); if (seq !~ /NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/) print $0 "\n"seq}' | tr -s ' ' '_'> $dirname/${basename%.*}_res/${basename%.*}.fasta && \
    echo -e "It takes about 1.2h to process a 2Gb fastq.gz on personal computers, please wait...\nStart reads_blast..." && \
    time $bin/programs/blast-2.11.0+/bin/blastn.exe  -num_threads 4 -outfmt 7 -db $bin/databases/db -query $dirname/${basename%.*}_res/${basename%.*}.fasta  -out $dirname/${basename%.*}_res/${basename%.*}_blastRes.xls
    echo "Start sorting by bit score..." && \
    sort -nrk 12 $dirname/${basename%.*}_res/${basename%.*}_blastRes.xls | awk '$4>80{print $0}' > $dirname/${basename%.*}_res/${basename%.*}_blastResSorted.xls && \
    awk '{if (!keys[$1]) print $0; keys[$1] = 1;}'  $dirname/${basename%.*}_res/${basename%.*}_blastResSorted.xls > $dirname/${basename%.*}_res/${basename%.*}_blastResSortedUniq.xls && \
    perl $bin/programs/reads_blast_Add_idName.pl   $bin/databases/db_idName   $dirname/${basename%.*}_res/${basename%.*}_blastResSortedUniq.xls > $dirname/${basename%.*}_res/${basename%.*}_blastResSortedUniqAddidName.xls && \
    rm $dirname/${basename%.*}_res/${basename%.*}_blastRes.xls $dirname/${basename%.*}_res/${basename%.*}_blastResSorted.xls && \
    awk '{sum[$2]+=1}END{for(i in sum)print i"\t"sum[i]}'   $dirname/${basename%.*}_res/${basename%.*}_blastResSortedUniq.xls | awk '/gi/||/gb/ {print}' |sort -nrk 2 >  $dirname/${basename%.*}_res/${basename%.*}_readsStat.xls && \
    perl $bin/programs/reads_blast_Add_idName.pl   $bin/databases/db_idName   $dirname/${basename%.*}_res/${basename%.*}_readsStat.xls > $dirname/${basename%.*}_res/${basename%.*}_readsStatAddidName.xls && \
    rm $dirname/${basename%.*}_res/${basename%.*}_readsStat.xls && \
    echo "Species table was generated at `date`!" && \

    awk 'NR==FNR{a[$1]=$0;next}{if ($2 in a) print  $0 "\t" a[$2]}'  $dirname/${basename%.*}_res/${basename%.*}_readsStatAddidName.xls  $dirname/${basename%.*}_res/${basename%.*}_blastResSortedUniq.xls |tr -s ',\ '  '_' | sort -rk 14   > $dirname/${basename%.*}_res/${basename%.*}_eachReadDetail.xls
    echo "Start extracting fasta from each_read_detail..." && \
    sed -i '1i\SWUN_ZYX'  $dirname/${basename%.*}_res/${basename%.*}_eachReadDetail.xls && \
    perl $bin/programs/fromContigsNameGetFasta.pl   $dirname/${basename%.*}_res/${basename%.*}_eachReadDetail.xls  $dirname/${basename%.*}_res/${basename%.*}.fasta > $dirname/${basename%.*}_res/${basename%.*}_matched.fa && \
    rm  $dirname/${basename%.*}_res/${basename%.*}.fasta && \

    echo "File $dirname/${basename%.*}_res/${basename%.*}.fasta was removed..." && \

    sed -i "s/\t/_/14" $dirname/${basename%.*}_res/${basename%.*}_eachReadDetail.xls && \
    mkdir $dirname/${basename%.*}_res/fasta_datasets && \
	echo "$dirname/${basename%.*}_res/fasta_datasets was generated."  && \
    for x in `awk '{print $14}'  $dirname/${basename%.*}_res/${basename%.*}_eachReadDetail.xls | uniq `
    do
    grep  $x $dirname/${basename%.*}_res/${basename%.*}_eachReadDetail.xls >> $dirname/${basename%.*}_res/fasta_datasets/${x}.xls
    done
    echo " $dirname/${basename%.*}_res/fasta_datasets/excels were generated, start making result table..." && \

###############     Std. Dev of Reads       #################
    for i in `ls $dirname/${basename%.*}_res/fasta_datasets/*xls`
    do
        species=`basename  $i`
        avr3=`cat $i|awk '{sum+=$3} END {print sum/NR}' `
        avr4=`cat $i|awk '{sum+=$4} END {print sum/NR}' `
        readsNum=`wc -l $i`
        sd=`awk '{print $9}'  $i | awk '{x[NR]=$0; s+=$0; n++} END{a=s/n; for (i in x){ss += (x[i]-a)^2} sd = sqrt(ss/n); print sd}' `
        printf "${species%.*}\t%f\t%f\t%d\t%f\n " "$avr3" "$avr4" "$readsNum" "$sd" >> $dirname/${basename%.*}_res/fasta_datasets/STDEVP.xls
    done
    sort -nrk 4 $dirname/${basename%.*}_res/fasta_datasets/STDEVP.xls > $dirname/${basename%.*}_res/${basename%.*}_result.xls && \
    sed -i '1i\Species\tMeanIdentity\tMeanQueryCover\tReadsNum\tSt_Dev_of_Reads_Location'  $dirname/${basename%.*}_res/${basename%.*}_result.xls && \
    echo "Standard deviation of reads location was generated!" && \
	
#############     extract reads from blastn_res      ##############
    for file in `ls $dirname/${basename%.*}_res/fasta_datasets/*.xls`
    do
            if [ ! -f $file ];then
                    echo "$file not exsist..."
            else
                    sed -i '1i\winVir_result'  $file
                    echo -e "Start extracting reads of `basename ${file%.*}`..."
                    perl $bin/programs/fromContigsNameGetFasta.pl $file $dirname/${basename%.*}_res/${basename%.*}_matched.fa >> ${file%.*}.fas
            fi
    done && \
	
    rm  $dirname/${basename%.*}_res/${basename%.*}_readsStatAddidName.xls  $dirname/${basename%.*}_res/${basename%.*}_blastResSortedUniq.xls $dirname/${basename%.*}_res/${basename%.*}_eachReadDetail.xls  $dirname/${basename%.*}_res/${basename%.*}_matched.fa $dirname/${basename%.*}_res/${basename%.*}_blastResSortedUniqAddidName.xls   $dirname/${basename%.*}_res/fasta_datasets/STDEVP.fas   $dirname/${basename%.*}_res/fasta_datasets/*xls && \
	#mv  $dirname/${basename%.*}_res/${basename%.*}_result.xls  $dirname/${basename%.*}_res/${basename%.*}_result.xls 
    echo -e "\n\nwinVir:\n       Conguratulation! Reads corresponding to each species were extracted, you can assemble them using SeqMan...\n       Low \"MeanQueryCover\" and \"St_Dev_of_Reads_Location\" suggests the result might be a false possitive...\n\n       Result files were generated in \"$dirname/${basename%.*}_res\"\n       All done at `date`, Bye Bye~"

else
    rm $dirname/${basename%.*}_res/${basename%.*}_matched.fa  $dirname/${basename%.*}_res/${basename%.*}_readsStatAddidName.xls  $dirname/${basename%.*}_res/${basename%.*}_eachReadDetail.xls  $dirname/${basename%.*}_res/${basename%.*}_blastResSortedUniq.xls  $dirname/${basename%.*}_res/${basename%.*}_blastResSortedUniqAddidName.xls 
	echo -e "\nwinVir:\n       File \"$dirname/${basename%.*}_res\" already exists...\n       Please remove it and try again..."
fi 
