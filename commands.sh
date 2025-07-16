### Activate Conda environment
#conda activate trycycler_env
#conda install minimap2=2.24
#conda list -n phame_env > trycycler_env_packages.txt
#conda env export > trycycler_env.yaml

### Avoid filling the disk with huge temporary files
export TMPDIR=/data/djs217/tmp

### Download FASTQ files of genomic sequence reads from Sequence Read Archive
#for i in SRR3254744 SRR33638550 SRR33638551 SRR34108134 SRR3254743 SRR33638434 SRR34103947; do
#    echo $i
#    fasterq-dump $i
#    gzip $i*.fastq
#done

### Rename the ONT FASTQ files
#mv SRR34103947.fastq.gz Cala2.SRR34103947.fastq.gz 
#mv SRR34108134.fastq.gz Noks1.SRR34108134.fastq.gz 

### Apply QC on the ONT FASTQ files
#filtlong --min_length 1000 --keep_percent 95 Cala2.SRR34103947.fastq.gz | gzip > Cala2.SRR34103947.filtlong.fastq.gz 
#filtlong --min_length 1000 --keep_percent 95 Noks1.SRR34108134.fastq.gz | gzip > Noks1.SRR34108134.filtlong.fastq.gz

### To save space, delete the original FASTQ files
#rm Cala2.SRR34103947.fastq.gz
#rm Noks1.SRR34108134.fastq.gz 

### Generate sub-samples of ONT sequence reads
#trycycler subsample --reads Cala2.SRR34103947.filtlong.fastq.gz --out_dir Cala2.read_subsets --count 12 --threads 12 --genome_size 80m
#trycycler subsample --reads Noks1.SRR34108134.filtlong.fastq.gz --out_dir Noks1.read_subsets --count 12 --threads 12 --genome_size 80m

### gzip the subsample files that will be used for long-read assemblies
#for i in 01 02 03 04 05 06 07 08 09 10 11 12 ; do
#    echo $i
#    gzip Cala2.read_subsets/sample_$i.fastq
#    gzip Noks1.read_subsets/$i
#done

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

    ### Raven assembly - takes too much disk space for remporary files - generates relatively poor assemblies
    #for i in 02 04 06 10 ; do
    #	if [ ! -f "$strain.assemblies/assembly_$i.fasta" ]; then
    #	    echo Do Raven assembly $i for $strain
    #	    raven --threads 12 --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_$i.gfa $strain.read_subsets/sample_$i.fastq.gz > $strain.assemblies/assembly_$i.fasta
    #	else
    #	    echo "$strain.assemblies/assembly_$i.fasta already exists, so skip this Raven assembly"
    #	fi
    #   done
    
    ### Canu assembly - takes too much disk space for temporary files and is very slow
    #  for i in 08 11 12 05; do
    #	if [ ! -f "$strain.assemblies/assembly_"$i".fasta" ]; then
    #	    echo Do Canu assembly $i for $strain
    #	    canu -p canu -d $i.canu_temp -fast genomeSize=80m useGrid=false maxThreads=12 -nanopore $strain.read_subsets/sample_"$i".fastq.gz
    #	    canu_trim.py $i.canu_temp/canu.contigs.fasta > $strain.assemblies/assembly_"$i".fasta
    #	    rm -rf $i.canu_temp
    #	else
    #	    echo "$strain.assemblies/assembly_$i.fasta already exists, so skip this Canu assembly"
    #	fi
    #   done
done
