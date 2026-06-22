# using star for RNAseq mapping
# RNAseq singal end
echo "usage:"
echo "bash ????.sh quality_cutoff(20/30) species(hg38. ce11)"
if [ ! -n "$2" ];then
	echo "no enough parameters"
	exit 8
fi

species=$2  

source ~/.bashrc
source activate qc
if [[ $1 == 20 ]] || [[ $1 == 30 ]]
then
    echo $1
else
    echo "without quality_cutoff(20/30)"
    exit 8
fi
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
    trim_galore -q $1 fq/${line}.fq.gz -o trimmed &
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
    STAR --genomeDir ~/reference/${species}/star_index/ --runThreadN 82 \
    --readFilesIn trimmed/${line}.fq.gz --outFileNamePrefix aligned/${line} \
    --outSAMattributes Standard --readFilesCommand zcat \
    --outFilterMultimapNmax 1 
    mv aligned/${line}Aligned.out.sam aligned/${line}.sam
    mv aligned/${line}Log* fastqc
    samtools sort -@ 82 aligned/${line}.sam -o aligned/${line}_o1.bam
    rm aligned/${line}.sam
    samtools view -@ 82 -h aligned/${line}_o1.bam >> aligned/${line}.sam
    # remove low quality and duplicates
    samtools view -@ 82 -q 30 -F 1024 aligned/${line}_o1.bam >> aligned/${line}.sam
    samtools sort -@ 82 aligned/${line}.sam -o aligned/${line}.bam
    rm aligned/${line}.sam
    rm aligned/${line}_o1.bam
    samtools index -@ 82 aligned/${line}.bam
done

conda activate

cd aligned
cp ../*.csv ./
mkdir ../bw
mkdir ../bw_nosep
for i in *.bam
do
    bamCoverage -b $i -o ../bw/${i%.*}_fw.bw --normalizeUsing RPKM --binSize 10 --filterRNAstrand forward -p 82
    bamCoverage -b $i -o ../bw/${i%.*}_rev.bw --normalizeUsing RPKM --binSize 10 --filterRNAstrand reverse -p 82
    bamCoverage -b $i -o ../bw_nosep/${i%.*}.bw --normalizeUsing RPKM --binSize 10 -p 82
done
cd ../bw
for i in *_rev.bw
do
    bash ~/shell/neg_bw.sh $i &
done
cd ..
multiqc fastqc