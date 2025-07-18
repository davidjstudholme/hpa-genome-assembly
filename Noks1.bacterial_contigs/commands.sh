conda activate circlator_env
conda list -n circlator_env > circlator_env_packages.txt
conda env export > circlator_env.yaml

for i in 387 609 721 799 508; do
    echo contig $i
    circlator fixstart bacterial.contig_"$i".fasta  contig_"$i".circlator
done

gzip *.fasta *.fa 
