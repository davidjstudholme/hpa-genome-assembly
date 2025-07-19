conda create -n seqkit_env
conda activate seqkit_env
conda install seqkit

conda list -n seqkit_env > seqkit_env_packages.txt
conda env export > seqkit_env.yaml

seqkit grep -v -f contigs_for_removal.txt ../Cala2.07.pypolca/pypolca_corrected.fasta > Cala2.assembly.cleaned.fasta

gzip Cala2.assembly.cleaned.fasta

