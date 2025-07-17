### Activate Conda environment
#conda activate trycycler_env
conda list -n phame_env > trycycler_env_packages.txt
conda env export > trycycler_env.yaml

### Avoid filling the disk with huge temporary files
export TMPDIR=/data/djs217/tmp

### Download FASTQ files of genomic sequence reads from Sequence Read Archive
for i in SRR3254744 SRR33638550 SRR33638551 SRR34108134 SRR3254743 SRR33638434 SRR34103947; do
    echo $i
    fasterq-dump $i
    gzip $i*.fastq
done

### Rename the ONT FASTQ files
mv SRR34103947.fastq.gz Cala2.SRR34103947.fastq.gz 
mv SRR34108134.fastq.gz Noks1.SRR34108134.fastq.gz 

### Apply QC on the ONT FASTQ files
filtlong --min_length 1000 --keep_percent 95 Cala2.SRR34103947.fastq.gz | gzip > Cala2.SRR34103947.filtlong.fastq.gz 
filtlong --min_length 1000 --keep_percent 95 Noks1.SRR34108134.fastq.gz | gzip > Noks1.SRR34108134.filtlong.fastq.gz

### To save space, delete the original FASTQ files
rm Cala2.SRR34103947.fastq.gz
rm Noks1.SRR34108134.fastq.gz 

### Generate sub-samples of ONT sequence reads
trycycler subsample --reads Cala2.SRR34103947.filtlong.fastq.gz --out_dir Cala2.read_subsets --count 12 --threads 12 --genome_size 80m
trycycler subsample --reads Noks1.SRR34108134.filtlong.fastq.gz --out_dir Noks1.read_subsets --count 12 --threads 12 --genome_size 80m

### gzip the subsample files that will be used for long-read assemblies
#for i in 01 02 03 04 05 06 07 08 09 10 11 12 ; do
    echo $i
    gzip Cala2.read_subsets/sample_$i.fastq
    gzip Noks1.read_subsets/$i
done

### Perform long-read de-novo assemblies
for strain in Noks1 Cala2; do
    echo $strain

    ### Generate assemblies for $strain
    mkdir $strain.assemblies
    
    ### Flye assembly
    for i in 01 02  03 04 05 06 07 08 09 10 11 12; do
	if [ ! -f "$strain.assemblies/assembly_$i.gfa" ]; then
	    echo Do Flye assembly $i for $strain
	    flye --nano-hq $strain.read_subsets/sample_$i.fastq.gz --threads 12 --out-dir assembly_$i --genome-size 80m
	    cp assembly_$i/assembly.fasta $strain.assemblies/assembly_$i.fasta
	    cp assembly_$i/assembly_graph.gfa $strain.assemblies/assembly_$i.gfa
	    rm -r assembly_$i
	else
	    echo "$strain.assemblies/assembly_$i.fasta already exists, so skip this Flye assembly"
	fi
    done
done


### Set up BUSCO
conda create -n busco_env -c bioconda -c conda-forge busco
conda activate busco_env
conda list -n busco_env > busco_env_packages.txt
conda env export > busco_env.yaml
https://busco-data.ezlab.org/v5/data/lineages/stramenopiles_odb12.2025-07-01.tar.gz
tar xvf stramenopiles_odb12.2025-07-01.tar.gz


### QC on Nokes1 assemblies, using BUSCO
for i in 01 02  03 04 05 06 07 08 09 10 11 12; do
    echo running BUSCO on Noks1 assembly $i
    busco -i Noks1.assemblies/assembly_"$i".fasta -l stramenopiles_odb10 -o Noks1."$i".busco -m genome --cpu 8 -f
done

### Make a summary plot of BUSCO results
mkdir Noks1.all.busco
cd  Noks1.all.busco
ln -s ../Noks1.*.busco/*.json .
cd -
busco --plot Noks1.all.busco



### Set up mash and seaborn
conda create -n mash_env
conda activate mash_env
conda install -c bioconda mash
pip install pandas seaborn matplotlib


### Run Mash pipeline

# Set working directory
WORKDIR="Noks1.assemblies"
MASH_OUT="Noks1.mash_dist.tsv"
MASH_MATRIX="Noks1.mash_distance_matrix.csv"
MASH_HEATMAP="Noks1.mash_heatmap.png"

echo "[Step 1] Creating Mash sketch..."
mash sketch -o mash_sketches $WORKDIR/*.fa*

echo "[Step 2] Calculating pairwise distances..."
mash dist mash_sketches.msh mash_sketches.msh > $MASH_OUT

echo "[Step 3] Creating distance matrix and heatmap..."
python3 mash_to_matrix.py $MASH_OUT $MASH_MATRIX $MASH_HEATMAP

echo "✅ Done. Outputs:"
echo "- $MASH_OUT"
echo "- $MASH_MATRIX"
echo "- $MASH_HEATMAP"




conda activate trycycler_env
medaka_consensus -i Noks1.SRR34108134.filtlong.fastq.gz -d Noks1.assemblies/assembly_06.fasta -o Noks1.06.medaka_output  -t 8 -m r1041_e82_400bps_sup_v4.3.0








### Set up Kraken
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz
tar -xvzf k2_standard_08gb_20250402.tar.gz
kdir kraken
mv k2_standard_08gb_20250402.tar.gz kraken/
mv database* kraken/
mv hash.k2d inspect.txt opts.k2d names.dmp nodes.dmp kraken/
mv taxo.k2d kraken/
mv  library_report.tsv seqid2taxid.map unmapped_accessions.txt kraken/
mv ktaxonomy.tsv kraken/

conda create -n kraken_env
conda activate kraken_env
conda install -c bioconda kraken2

### Run Kraken on Noks1 assembly
kraken2 --db ./kraken --threads 8 --output Noks1.kraken2_output.txt --report Noks1.06.kraken2_report.txt Noks1.assemblies/assembly_06.fasta


###  Get all taxids descending from Bacteria (taxid 2)
awk -F '\t\|\t|\t' '$3 == 2 {print $1}' kraken/nodes.dmp > bacterial_taxids.txt

###  Get all taxids descending from Fungi (taxid 4751)    
awk -F '\t\|\t|\t' '$3 == 4751 {print $1}' kraken/nodes.dmp > fungal_taxids.txt

### Get a list of the Noks1 contigs that match bacterial taxa
awk 'NR==FNR {bact[$1]; next} $1 == "C" && ($3 in bact) {print $2, $3}' bacterial_taxids.txt Noks1.06.kraken2_output.txt > Noks1.06.bacterial_contigs_list.txt

### Get a list of the Noks1 contigs that match fungal taxa 
awk 'NR==FNR {fungi[$1]; next} $1 == "C" && ($3 in fungi) {print $2, $3}' fungal_taxids.txt Noks1.06.kraken2_output.txt > Noks1.06.fungal_contigs_list.txt

