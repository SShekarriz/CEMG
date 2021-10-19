#!/bin/bash
####################################################################################
DATE=`date '+%Y-%m-%d-%H:%M'`
exec > >(tee -i CEMG_PIPER_${DATE}.log)
exec 2>&1
####################################################################################

Don=SHCM15
proj=$PWD
org=$PWD/raw-reads
trim=$PWD/raw-reads_trimmed
bin=$PWD/bin
assem=$PWD/assembled
Res=$PWD/results
mkdir -p $trim $assem $Res



echo "#####################################"
echo "####### Req softwares         #######"
Trimomatic="/dataone/common/software/Trimmomatic-0.38/trimmomatic-0.38.jar"
echo $Trimmomatic
Samtools="/dataone/common/software/samtools-1.11/samtools"
echo $Samtools
InterLeave=$PWD/bin/interleave_fastq.py
echo $InterLeave
Spades=$PWD/bin/SPAdes-3.15.2-Linux/bin/spades.py
echo $Spades
LongExractor=$PWD/bin/LongContigExtractor.py
echo $LongExractor

echo "#####################################"


echo "#######################################"
echo "Trimming the fastq files; trimmomatic"
echo "#######################################"
Adapters=$bin/NEB_Costum.fa
mkdir -p $trim
cd $org
for R1 in *_R1.fastq.gz
do
        sample=${R1%_R1*}
	R2=${sample}_R2.fastq.gz
	P1=$trim/${sample}_1P.fastq.gz
	P2=$trim/${sample}_2P.fastq.gz
	U1=$trim/${sample}_1U.fastq.gz
	U2=$trim/${sample}_2U.fastq.gz
	if [ -f $trim/${sample}_1P.fastq.gz ]; then
        echo "The file '$trim/${sample}_1P.fastq.gz' exists."
        else
        echo "The file '$trim/${sample}_1P.fastq.gz' is not found."
        java -jar $Trimomatic PE -threads 30 \
        $R1 $R2 $P1 $U1 $P2 $U2 \
        ILLUMINACLIP:$Adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
        fi
done

echo "####################################################"
echo "#######building CEMG library, co-assem pp + stool ##"
echo "####################################################"

Don_lib=$assem/${Don}_CEMG
Don_temp=$assem/${Don}_temp
mkdir -p $Don_temp

if [ -f ${Don_lib}_inter.fastq.gz ]; then
	echo "The file '${Don_lib}_inter.fastq.gz' exists."
else
        echo "The file '${Don_lib}_inter.fastq.gz' is not found."

	cd $trim
	cp *P.fastq.gz $Don_temp/
	gunzip $Don_temp/*.fastq.gz
	
	for sample_R1 in *_1P.fastq.gz; do
                sample_R2=${sample_R1%_1P*}_2P.fastq.gz
		sample=${sample_R1%_1P*}
                echo "${Don} contains: $sample_R1, $sample_R2"
                python2.7 $InterLeave $Don_temp/${sample}_1P.fastq $Don_temp/${sample}_2P.fastq $Don_temp/${sample}_inter.fastq
	done
	echo "removing SHCM15ana10b_inter.fastq from CEMG pool, this is not part of CEMG"
	mv $Don_temp/SHCM15ana10b_inter.fastq $Don_temp/SHCM15ana10b_inter_OUTofCEMG.fastq
        echo "'$Don_temp/*_inter.fastq' into $Don_lib"
        cat $Don_temp/*_inter.fastq >> ${Don_lib}_inter.fastq
        echo "Cleaning up the tem directory"
        gzip $Don_temp/*_inter.fastq
	gzip ${Don_lib}_inter.fastq
        rm $Don_temp/*_1P*.fastq
	rm $Don_temp/*_2P*.fastq

fi

echo "#####################################"
echo "#### Assembling the donor libs ######"
echo "#####################################"

if [ -f ${Don_lib}_spade/contigs.fasta ]; then
	echo "The file '${Don_lib}_spade/contigs.fasta' exists."
else
	echo "The file '${Don_lib}_spade/contigs.fasta' is not found."
	echo "assembling '${Don_lib}_inter.fastq'"
        python3 $Spades --meta --pe1-12 ${Don_lib}_inter.fastq.gz -o ${Don_lib}_spade \
                -t 150 --memory 2000
fi

echo "########################################"
echo "####### single stool assembly ## #### ##"
echo "########################################"

Don_stool=$assem/${Don}_DMG
out=${Don_stool}_spade/contigs.fasta
if [ -f $out ]; then
	echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
	#python3 $Spades --meta --pe1-12 $Don_temp/${Don}Stool_inter.fastq -o ${Don_stool}_spade \
        #-t 50 --memory 1500
fi


echo "########################################"
echo "## Adding lib name to contig headers  ##"
echo "########################################"
assemQ=$assem/quality
mkdir -p $assemQ

out=$assemQ/${Don}_CEMG_contigs_edit_headers.txt
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."

        cd $assem
        for dir in *_spade; do
                echo $dir
                sample=${dir%_spade}
                echo $sample
                sed "/^>/s/^>/>${sample}__/g" $dir/contigs.fasta > $dir/contigs_edit.fasta
                grep -e "^>" $dir/contigs_edit.fasta > $assemQ/${sample}_contigs_edit_headers.txt
        done
fi

echo "########################################"
echo "## Comparing contigs cumulative plots ##"
echo "########################################"

out=$Res/assembly_cumulative.png
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
        cd $assem
        Rscript $bin/CumuContigs.R quality $Res/assembly_cumulative.png
fi


echo "################################"
       	echo "Donors: $Don"
        Don_temp=$assem/${Don}_temp
        contig_file=$assem/${Don}_CEMG_spade/contigs_edit.fasta
        echo $contig_file

        echo "################################"
        echo "#####selecting large contigs####"
        echo "################################"

        dRef=${contig_file%.fasta}_1k.fa
        if [ -f $dRef ]; then
        echo "The file '$dRef' exists."
        else
        echo "The file '$dRef' is not found."
        echo "Selecting contig >= 1kb"
        python2.7 $LongExractor $contig_file
        mv ${contig_file}_out $dRef
        fi


        echo "################################"
        echo "#indexing contigs Large contigs#"
        echo "################################"

        if [ -f ${dRef}.fai ]; then
        echo "The file '${dRef}.fai' exists."
        else
        echo "The file '${dRef}.fai' is not found."
        echo "Indexing the contig file now"
        bwa index $dRef
        samtools faidx $dRef
        fi





