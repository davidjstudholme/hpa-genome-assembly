conda create -n seqkit_env
conda activate seqkit_env
conda install seqkit

conda list -n seqkit_env > seqkit_env_packages.txt
conda env export > seqkit_env.yaml

seqkit grep -v -f contigs_for_removal.txt ../Noks1.06.pypolca/pypolca_corrected.fasta > Noks1.assembly.cleaned.fasta

gzip Noks1.assembly.cleaned.fasta

