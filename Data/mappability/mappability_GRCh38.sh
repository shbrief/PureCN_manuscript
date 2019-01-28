# To generate a mappability file with the GEM library
# Calculate mappability, set kmer size to length of mapped reads
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#wget https://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2/download
#tar xvf download

THREADS=24
KMER=76

export PATH=$PATH:/home/sehyun/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin
PREF="/home/sehyun/reference/GRCh38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set"
REFERENCE="${PREF}.fna"

gem-indexer -T ${THREADS} -i ${REFERENCE} -o ${PREF}_index
gem-mappability -T ${THREADS} -I ${PREF}_index.gem -l ${KMER} -o ${PREF}_${KMER} -m 2 -e 2
gem-2-wig -I ${PREF}_index.gem -i ${PREF}_${KMER}.mappability -o ${PREF}_${KMER}

# Convert to bigWig format, for example using the UCSC wigToBigWig tool
cut -f1,2 ${REFERENCE}.fai > ${PREF}.sizes 

# I found the unexpected letter, "AC" in my .wig file --> remove it
cp ${PREF}_${KMER}.wig GRCh38_no_alt_${KMER}.wig
sed -e s/AC//g -i GRCh38_no_alt_${KMER}.wig
/home/sehyun/tools/wigToBigWig GRCh38_no_alt_${KMER}.wig ${PREF}.sizes ${PREF}_${KMER}.bw

