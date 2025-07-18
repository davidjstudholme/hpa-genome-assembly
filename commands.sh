### Activate Conda environment
#conda activate trycycler_env
conda list -n trycycler_env > trycycler_env_packages.txt
conda env export > trycycler_env.yaml

### Avoid filling the disk with huge temporary files
export TMPDIR=/data/djs217/tmp

### Download FASTQ files of genomic sequence reads from Sequence Read Archive
for i in SRR3254744 SRR33638550 SRR33638551 SRR34108134 SRR3254743 SRR33638434 SRR34103947; do
    echo $i
    fasterq-dump $i
    gzip $i*.fastq
done

### Rename the Illumina FASTQ files
mv SRR3254743_1.fastq.gz Cala2.SRR3254743_1.fastq.gz
mv SRR3254743_2.fastq.gz Cala2.SRR3254743_2.fastq.gz
mv SRR3254744_1.fastq.gz Noks1.SRR3254744_1.fastq.gz
mv SRR3254744_2.fastq.gz Noks1.SRR3254744_2.fastq.gz
mv SRR33638434_1.fastq.gz Cala2.SRR33638434_1.fastq.gz
mv SRR33638434_2.fastq.gz Cala2.SRR33638434_2.fastq.gz
mv SRR33638550_1.fastq.gz Noks1.SRR33638550_1.fastq.gz
mv SRR33638550_2.fastq.gz Noks1.SRR33638550_2.fastq.gz
mv SRR33638551_1.fastq.gz Noks1.SRR33638551_1.fastq.gz
mv SRR33638551_2.fastq.gz Noks1.SRR33638551_2.fastq.gz

### Perform QC on Illumina FASTQ files
trim_galore -q 30 --paired Noks1.SRR33638551_1.fastq.gz Noks1.SRR33638551_2.fastq.gz
trim_galore -q 30 --paired Cala2.SRR33638434_1.fastq.gz Cala2.SRR33638434_2.fastq.gz

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
for i in 01 02 03 04 05 06 07 08 09 10 11 12 ; do
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
wget https://busco-data.ezlab.org/v5/data/lineages/stramenopiles_odb12.2025-07-01.tar.gz
tar zxvf stramenopiles_odb12.2025-07-01.tar.gz
conda create -n busco_env -c bioconda -c conda-forge busco
conda activate busco_env
conda list -n busco_env > busco_env_packages.txt
conda env export > busco_env.yaml

### QC on Noks1 assemblies, using BUSCO
for i in 01 02 03 04 05 06 07 08 09 10 11 12; do
    echo running BUSCO on Noks1 assembly $i
    busco -i Noks1.assemblies/assembly_"$i".fasta -l stramenopiles_odb10 -o Noks1."$i".busco -m genome --cpu 8 -f
done

### QC on Cala2 assemblies, using BUSCO
for i in 01 02 03 04 05 06 07 08 09 10 11 12; do
    echo running BUSCO on Cala2 assembly $i
    busco -i Cala2.assemblies/assembly_"$i".fasta -l stramenopiles_odb10 -o Cala2."$i".busco -m genome --cpu 8 -f
done

### Make a summary plot of BUSCO results
mkdir Noks1.all.busco
cd  Noks1.all.busco
ln -s ../Noks1.*.busco/*.json .
cd -
busco --plot Noks1.all.busco

### Make a summary plot of BUSCO results
mkdir Cala2.all.busco
cd  Cala2.all.busco
ln -s ../Cala2.*.busco/*.json .
cd -
busco --plot Cala2.all.busco

### Set up mash and seaborn
conda create -n mash_env
conda activate mash_env
conda install -c bioconda mash
pip install pandas seaborn matplotlib
conda list -n mash_env > mash_env_packages.txt
conda env export > mash_env.yaml

### Run Mash pipeline
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

### Tidy up
mkdir Noks1.busco
mv *Noks1*.busco Noks1.busco/

### Run Mash pipeline
WORKDIR="Cala2.assemblies"
MASH_OUT="Cala2.mash_dist.tsv"
MASH_MATRIX="Cala2.mash_distance_matrix.csv"
MASH_HEATMAP="Cala2.mash_heatmap.png"
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

### Tidy up
mkdir Cala2.busco
mv *Cala2*.busco Cala2.busco/

### Polishing assembly with long reads - Medaka
conda activate trycycler_env
medaka_consensus -i Noks1.SRR34108134.filtlong.fastq.gz -d Noks1.assemblies/assembly_06.fasta \
		 -o Noks1.06.medaka_output -t 8 -m r1041_e82_400bps_sup_v4.3.0

### Polish with short reads - PyPolca
pypolca run --force --careful -a Noks1.06.medaka_output/consensus.fasta \
	-1 Noks1.SRR33638551_1_val_1.fq.gz \
	-2 Noks1.SRR33638551_2_val_2.fq.gz  \
	-t 16 -o Noks1.06.pypolca
cp Noks1.06.pypolca/pypolca_corrected.fasta Noks1.06.medaka.pypolca.fasta
gzip Noks1.06.medaka.pypolca.fasta

### Polishing assembly with long reads - Medaka
conda activate trycycler_env
medaka_consensus -i Cala2.SRR34103947.filtlong.fastq.gz -d Cala2.assemblies/assembly_07.fasta \
		 -o Cala2.07.medaka_output -t 8 -m r1041_e82_400bps_sup_v4.3.0

### Polish with short reads - PyPolca
pypolca run --force --careful -a Noks1.06.medaka_output/consensus.fasta \
	-1 Cala2.SRR33638434_1_val_1.fq.gz \
	-2 Cala2.SRR33638434_2_val_2.fq.gz \
	-t 16 -o Cala2.06.pypolca
cp Cala2.07.pypolca/pypolca_corrected.fasta Cala2.07.medaka.pypolca.fasta
gzip Cala2.07.medaka.pypolca.fasta
