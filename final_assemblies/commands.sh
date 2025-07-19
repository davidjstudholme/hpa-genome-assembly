### Set up BUSCO
wget https://busco-data.ezlab.org/v5/data/lineages/stramenopiles_odb12.2025-07-01.tar.gz
tar zxvf stramenopiles_odb12.2025-07-01.tar.gz
#conda create -n busco_env -c bioconda -c conda-forge busco
conda activate busco_env
conda list -n busco_env > busco_env_packages.txt
conda env export > busco_env.yaml

### QC on Noks1 assemblies, using BUSCO
for i in *.fna *.fasta; do
    echo running BUSCO on $i
    busco -i $i -l stramenopiles_odb10 -o $i.busco -m genome --cpu 8 -f
done

### Make a summary plot of BUSCO results
mkdir busco.all
cd  busco.all
ln -s ../*.busco/*.json .
cd -
busco --plot busco.all

### QUAST
conda activate quast_env
conda list -n quast_env > quast_env_packages.txt
conda env export > quast_env.yaml

quast.py *.fna *.fasta
