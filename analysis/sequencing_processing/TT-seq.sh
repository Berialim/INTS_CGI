# TT-seq single end
# using star for RNAseq mapping
# using spike in for normalize 
# RNAseq singal end
echo "usage:"
echo "bash ????.sh species(hg38. ce11)  spike-in-species(dm6)"
if [ ! -n "$3" ];then
	echo "no enough parameters"
	exit 8
fi

species=$1 
spike_in=$2


source ~/.bashrc
source activate qc

# need all fastq file in fq files
cd fq 
for i in *.fq.gz
do
    ii=`echo $i | sed "s/.*_L0._//g"`
    # remove the sequencing mark (FXXXX_L01_sample.fq.gz)
    mv $i $ii
done
# mv ambiguous reads and undecoded reads into useless folder
mkdir useless
mv ambiguous.fq.gz useless
mv undecoded.fq.gz useless
x=`ls *.fq.gz | sed 's/.fq.gz//g'`
cd ..
##trimadapter
mkdir trimmed
for line in $x
do
    trim_galore -q 30 fq/${line}.fq.gz -o trimmed &
done
wait
# remove the mark of trimmed 
cd trimmed
for i in *.fq.gz
do
    mv ${i} ${i%%_*}.fq.gz
done
cd ..
# # remove the rRNA reads
# mkdir de_rRNA
touch rRNA_report
for line in $x
do
    echo 'rRNA_sample:'$line >> rRNA_report
    bowtie -p 16 -x ~/reference/${species}/ribosome/rDNA -q trimmed/${line}.fq.gz -S useless.sam 2>>rRNA_report
done
rm *.sam

# get the quality of fastq
mkdir fastqc
for line in $x
do
    fastqc trimmed/${line}.fq.gz -o fastqc &
done
wait
source ~/.bashrc
source activate star
mkdir aligned
for line in $x
do
    STAR --genomeDir ~/reference/${species}/star_index/ --runThreadN 20 \
    --readFilesIn trimmed/${line}.fq.gz --outFileNamePrefix aligned/${line} \
    --outSAMattributes Standard --readFilesCommand zcat \
    --outFilterMultimapNmax 1 
    mv aligned/${line}Aligned.out.sam aligned/${line}.sam
    mv aligned/${line}Log* fastqc
    samtools sort -@ 16 aligned/${line}.sam -o aligned/${line}_o1.bam
    rm aligned/${line}.sam
    samtools view -@ 16 -h aligned/${line}_o1.bam >> aligned/${line}.sam
    # remove low quality and duplicates
    samtools view -@ 16 -q 30 -F 1024 aligned/${line}_o1.bam >> aligned/${line}.sam
    samtools sort -@ 16 aligned/${line}.sam -o aligned/${line}.bam
    rm aligned/${line}.sam
    rm aligned/${line}_o1.bam
    samtools index -@ 16 aligned/${line}.bam
done

mkdir temp_spikein
for line in $x
do
    echo "sample: ${line}_spikein" >> align_report
    STAR --genomeDir ~/reference/${spike_in}/star_index/ --runThreadN 20 \
    --readFilesIn   trimmed/${line}.fq.gz --outFileNamePrefix temp_spikein/${line}_spikein \
    --outSAMattributes Standard --readFilesCommand zcat	\
    --outFilterMultimapNmax 1 
    mv temp_spikein/${line}_spikeinAligned.out.sam temp_spikein/${line}.sam
    mv temp_spikein/${line}_spikeinLog* fastqc
    samtools sort -@ 16 temp_spikein/${line}.sam -o temp_spikein/${line}_o1.bam
    rm temp_spikein/${line}.sam
    samtools view -@ 16 -h temp_spikein/${line}_o1.bam >> temp_spikein/${line}.sam
    # remove low quality and duplicates
    samtools view -@ 16 -q 30 -F 1024 temp_spikein/${line}_o1.bam >> temp_spikein/${line}.sam
    samtools sort -@ 16 temp_spikein/${line}.sam -o temp_spikein/${line}.bam
    rm temp_spikein/${line}.sam
    rm temp_spikein/${line}_o1.bam
    samtools index -@ 16 temp_spikein/${line}.bam
done

conda activate

cd temp_spikein

multiBamSummary bins --bam *.bam --scalingFactors scalingfactors.txt -p 16

cd ..
# generate bw with spikein scalingfactor
mkdir bw
cat temp_spikein/scalingfactors.txt | while read line
do
    arr=(${line/'/\t'/ })
    i=`echo ${arr[0]} | sed "s/.bam//g"`
    sf=${arr[1]}
    bamCoverage -b aligned/$i.bam -o bw/${i}_fw.bw --normalizeUsing None \
    --binSize 10 --filterRNAstrand forward -p 16 --scaleFactor $sf
    bamCoverage -b aligned/$i.bam -o bw/${i}_rev.bw --normalizeUsing None \
    --binSize 10 --filterRNAstrand reverse -p 16 --scaleFactor $sf
    bamCoverage -b aligned/$i.bam -o bw/${i}.bw --normalizeUsing None \
    --binSize 10 -p 16 --scaleFactor $sf
done
# generate negative bw file
cd bw
for i in *_rev.bw
do
    bash ~/shell/neg_bw.sh $i &
done
cd ..
multiqc fastqc