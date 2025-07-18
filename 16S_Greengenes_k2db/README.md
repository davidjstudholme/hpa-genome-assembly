# 16S Greengenes Kraken2/Bracken Databsae
The following sections describes the contents of the database, 
programs/command lines used to generate the files, and 
how users can utilize these files to analyze sequenced reads 

Please note that this is a nucleotide 16S database. 

## Example Usage
Given paired 150bp reads, the following commands will estimate 
**Genus** abundance with the provided files. 
For clarity, each command line is split into different lines 
and annotated to explain each parameter. 

```
    kraken2 --db 16S_Greengenes_k2db 
        --threads 4                 #number of threads
        --report SAMPLE.kreport2    #kraken-style report (REQUIRED FOR BRACKEN)
        --paired SAMPLE_1.fq SAMPLE_2.fq > SAMPLE.kraken2
    bracken -d 16S_Greengenes_k2db    
        -i SAMPLE.kreport2          #REPORT output from kraken2 
        -o SAMPLE.bracken           #tab-delimited text file with read counts pre/post abundance calc
        -w SAMPLE_bracken.kreport2  #kraken-style report with bracken read counts
        -r 150                      #select read length (paired reads - use length of one mate) 
        -l G                        #level at which to calculate abundance (G = genus, S = species, etc) 
```

Users should not directly copy the above command lines. 
Select options best suitable for your computer resources/purposes.
Remove line breaks between parameters and remove all annotations (#).

## Files Included
Kraken 2 Database files
* `hash.k2d` 
* `opts.k2d`
* `taxo.k2d` 
Bracken Files
* `database50mers.kmer_distrib`
* `database75mers.kmer_distrib`
* `database100mers.kmer_distrib`
* `database150mers.kmer_distrib`
* `database200mers.kmer_distrib`
* `database250mers.kmer_distrib`
* `seqid2taxid.map`

## Programs/Versions Used
Kraken 2 version downloaded 03/05/2020 (https://github.com/DerrickWood/kraken2/commits/master/)
Bracken v2.5.2 (https://github.com/jenniferlu717/Bracken/releases) 

## Command Lines used to Generate Files
The Kraken 2 files (`*.k2d`) were generated using 
    ```
    kraken2-build --special greengenes --db 16S_Greengenes_k2db --threads 16 
    ``` 

Bracken files (`*.kmer_distrib`) were generated using
    ```
    bracken-build -k 35 -l 50 -d 16S_Greengenes_k2db -t 35 
    bracken-build -k 35 -l 75 -d 16S_Greengenes_k2db -t 35 
    bracken-build -k 35 -l 100 -d 16S_Greengenes_k2db -t 35 
    bracken-build -k 35 -l 150 -d 16S_Greengenes_k2db -t 35 
    bracken-build -k 35 -l 200 -d 16S_Greengenes_k2db -t 35 
    bracken-build -k 35 -l 250 -d 16S_Greengenes_k2db -t 35 
    ```

# Author Information
Jennifer Lu 
jennifer.lu717@gmail.com
[http://ccb.jhu.edu/people/jennifer.lu]

Last Updated: 03/25/2020 
