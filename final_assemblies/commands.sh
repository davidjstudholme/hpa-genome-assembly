### Set up BUSCO
wget https://busco-data.ezlab.org/v5/data/lineages/stramenopiles_odb12.2025-07-01.tar.gz
tar zxvf stramenopiles_odb12.2025-07-01.tar.gz
conda create -n busco_env -c bioconda -c conda-forge busco
conda activate busco_env
conda list -n busco_env > busco_env_packages.txt
conda env export > busco_env.yaml

### QC on Noks1 assemblies, using BUSCO
for i in *.fna *.fasta; do
    echo running BUSCO on $i
    busco -i $i -l stramenopiles_odb10 -o busco.$i -m genome --cpu 8 -f
done

### Make a summary plot of BUSCO results
mkdir busco.all
cd  busco.all
ln -s ../busco.*/*.json .
cd -
busco --plot busco.all

### Set up mash and seaborn
conda create -n mash_env
conda activate mash_env
conda install -c bioconda mash
pip install pandas seaborn matplotlib
conda list -n mash_env > mash_env_packages.txt
conda env export > mash_env.yaml

### Run Mash pipeline
WORKDIR="./"
MASH_OUT="mash_dist.tsv"
MASH_MATRIX="mash_distance_matrix.csv"
MASH_HEATMAP="mash_heatmap.png"
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
