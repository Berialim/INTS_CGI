mkdir trimmed
mkdir fastqc1

# move all fastq in fq dictionary
echo "usage:"
echo "bash chipseq_pe.sh species(e.g. hg38) fragment_length(from QSep) effectiveGenomeSize(for deeptools bamCoverage)"
if [ ! -n "$3" ];then
        echo "not enough parameter"
        exit 8
fi
species=$1
fragment_length=$2
genomesize=$3

source /home/fanlab/.bashrc
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
x=`ls *_1.fq.gz | sed 's/_1.fq.gz//g'`
cd ..
##trimadapter
mkdir trimmed
for line in $x
do
    trim_galore --paired -q 30 fq/${line}_1.fq.gz  fq/${line}_2.fq.gz  -o trimmed &
done
wait
# remove the mark of trimmed 
cd trimmed
for i in *.fq.gz
do
    mv ${i} ${i%%_val*}.fq.gz
done
cd ..

for i in trimmed/*.fq.gz
do 
    fastqc ${i} -o fastqc1& 
done
wait
multiqc fastqc1 -o qc
mkdir aligned
cd trimmed
for i in *_1.fq.gz
do
    x=`echo $i | sed "s/_1.fq.gz//g"`
    echo "sample: ${x}" >> ../qc/align_report
    # align to genome
    bowtie2 -p 82 -x /reference/$species/bowtie2/genome -1 ${i} -2 ${x}_2.fq.gz  -S ../aligned/${x}.sam 2>> ../qc/align_report
    # generate bam file from sam file 
    samtools sort -@ 82 -o ../aligned/${x}.bam  ../aligned/${x}.sam
    # remove duplicates
    picard MarkDuplicates I=../aligned/${x}.bam  O=../aligned/${x}.markdup.bam  M=../aligned/${x}.markdup.txt
    # remove temp files
    rm ../aligned/${x}.sam
    rm ../aligned/${x}.bam
    mv ../aligned/${x}.markdup.bam ../aligned/${x}.bam
    # generate bam index
    samtools index -@ 82 ../aligned/${x}.bam
done
cd ../qc
python /home/fanlab/code/alignment_se.py

# generatre bigwig file for read coverage in genome(visualize in igv)
mkdir ../bw
cd ../aligned
for i in *.bam
do
    bamCoverage --bam ${i} -o ../bw/${i%.*}.bw --binSize 10 --normalizeUsing RPGC \
     --effectiveGenomeSize $genomesize -e $fragment_length -p 82
done

rm -r trimmed