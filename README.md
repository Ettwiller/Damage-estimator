### Damage-estimator
For more information please read our pre-print in Biorxiv : http://biorxiv.org/content/early/2016/08/23/070334


##REQUIREMENT : 
Prior to running the Damage estimator workflow please download and install the following programs :

**SAMTOOLS** (http://samtools.sourceforge.net/)
**R** (https://www.r-project.org/)
**GGPLOT2** (http://ggplot2.org/ in R : install.packages("ggplot2", dependencies=TRUE))

##OVERVIEW :

This set of programs are designed to estimate the DNA damage when the DNA is sequenced using Illumina plateform on paired-end mode. 
There is 3 steps (starting from an aligned bam file) :

- Split the paired end reads into R1 and R2 using ```split_mapped_reads.pl``` (universal)
- Estimate the damage / damage across reads / damage across reads and context using ```estimate_damage.pl```, ```estimate_damage_location.pl``` and ```estimate_damage_location_context.pl``` respectively
- Plot the result using R.


##CONSIDERATIONS :
Damage estimation is based on the systematics mutation rate difference between the first in pair and the second in pair reads. Therefore it is essential that the sequencing is done using Illumina in a paired end mode. BWA mapping is recommended and mapping in paired-end mode is required. The estimation of damage is a global estimation based in a imbalance between R1 and R2 variant frequency (see publication). Thus a minimum of 2 million mapped read is the limit for detection. 


##TYPICAL WORKFLOW :

###Primary analysis : 
- Adaptor trimming (using your favorite algorithm).
- Mapping reads to the genome : BWA mem paired-end mode. Make a bam file from the resulting sam file (samtools view -bS sam_file | samtools sort - bam_file) and index the resulting bam file (samtools index bam_file). Other mapping algorithm can be used but no guaranteed. Local realignment using GATK can be optionally performed.  

###Split the reads to estimate the imbalance :

- Create first in pair mapped reads file (bam1) and second in pair mapped reads (bam2) and derived respective mpileup files (mpileup1 and mpileup2) using ```split_mapped_reads.pl```

###No context no position on the reads (just the basic degree of damage):

- Calculate damage using ```estimate_damage.pl``` 
- Concatenate all the damages into one single file for plotting.
- Plot result using ```plot_damage.R```


![alt tag](https://github.com/Ettwiller/Damage-estimator/blob/master/figures/Output_of_plot_damage_R.png)

###No context but random sampling :

- Calculate damage using ```random_sampling_and_estimate_damage.pl```
- Concatenate all the damage files into one single file for plotting.
- Plot result using ```plot_random_sampling_damage.R ```

![alt tag](https://github.com/Ettwiller/Damage-estimator/blob/master/figures/Output_of_plot_random_sampling_damage_R.png)



###No context but damage relative to read position: 

- Calculate damage using ```estimate_damage_location.pl```
- concatenate all the damages into one single file for plotting. 
- plot result using ```plot_damage_location.R```

![alt tag](https://github.com/Ettwiller/Damage-estimator/blob/master/figures/Output_of_plot_damage_location_R.png)

###Context and read positions : 

- Calculate damage using ```estimate_damage_location_context.pl```
- concatenate all the damages into one single file for plotting. 
- plot result using ```plot_damage_location_context.R```

![alt tag](https://github.com/Ettwiller/Damage-estimator/blob/master/figures/Output_of_plot_damage_location_context_R.png)


## EXAMPLES

### Getting the GIV scores (Similar to Figure 2) :
-```perl split_mapped_reads.pl -bam bam_file.bam -genome genome.fasta -mpileup1 file1.mpileup -mpileup2 file2.mpileup -Q 20 (DEFAULT 0) -q 20 (DEFAULT 10)```

-```perl estimate_damage.pl --mpileup1 file1.mpileup --mpileup2 file2.mpileup --id TE_shear --qualityscore 25 (DEFAULT 30) > TE_shear.damage```

-```Rscript --vanilla plot_damage.R TE_shear.damage TE_shear.png``` OR ```Rscript --vanilla plot_damage.R TE_shear.damage TE_shear.pdf```

### Getting the damage relative to read position (Similar to Figure 1C) :

-```perl split_mapped_reads.pl -bam bam_file.bam -genome genome.fasta -mpileup1 file1.mpileup -mpileup2 file2.mpileup -Q 20 (DEFAULT 0) -q 20 (DEFAULT 10)```

-```perl estimate_damage_location.pl --mpileup1 file1.mpileup --mpileup2 file2.mpileup --id TE_shear --qualityscore 25 (DEFAULT 30) --out TE_shear.damage```

-```Rscript --vanilla  plot_damage_location.R TE_shear.damage TE_shear.pdf

### Getting the damage relative to read position (Similar to Supplementary Figure 1) :

-```perl split_mapped_reads.pl -bam bam_file.bam -genome genome.fasta -mpileup1 file1.mpileup -mpileup2 file2.mpileup -Q 20 (DEFAULT 0) -q 20 (DEFAULT 10)```

-```perl estimate_damage_location_context.pl --mpileup1 file1.mpileup --mpileup2 file2.mpileup --id TE_shear --qualityscore 25 (DEFAULT 30) --out TE_shear.damage --context 2```

-```Rscript --vanilla  plot_damage_location_context.R TE_shear.damage TE_shear.pdf



##DETAILS OF THE ANALYSIS :

### 1. split_mapped_reads.pl :

EXAMPLE :
```perl split_mapped_reads.pl -bam bam_file.bam -genome genome.fasta -mpileup1 file1.mpileup -mpileup2 file2.mpileup -Q 20 (DEFAULT 0) -q 20 (DEFAULT 10)```

DESCRIPTION :
```split_mapped_reads.pl``` is the program to split the mapped first in paired reads into one temporary bam file and the mapped second in paired reads into another temporary bam file. From these files are derived mpileup file 1 and mpileup file 2 resperctively. 

OPTIONS :

```-bam``` (REQUIRED): Paired-end aligned reads in bam format.

```-genome``` (REQUIRED): fasta file containing the genome used for mapping. 

```-Q``` (OPTIONAL): Phred score quality threshold (Sanger encoding). Only keep the bases with a Q score above a given threshold (default 0).

```-q``` (OPTIONAL): mapping quality. Only keep the reads that passes a given threshold (default 10). 

OUTPUT :
The outputs of ```split_mapped_reads.pl``` is 2 mpileup files generated by samtools containing all the positions in the genome with at least one read. The file in -mpileup1 correspond to the first in paired reads and the file in -mpileup2 correspond to the second in paired reads. Both files are intermediate to be used by the ```estimate_damage_location.pl``` program.


### 2. estimate_damage.pl :

EXAMPLE :
```perl estimate_damage.pl --mpileup1 file1.mpileup --mpileup2 file2.mpileup --id water_shearing --qualityscore 25 (DEFAULT 30) > out_file ```

DESCRIPTION :

OPTIONS :

```--mpileup1``` REQUIRED mpileup file 1 (output of split_mapped_reads.pl)

```--mpileup2``` REQUIRED mpileup file 2 (output of split_mapped_reads.pl)

```--id``` REQUIRED id of the experiment. Can contain any letters or number. No space (ex: idx12)

```--qualityscore MIN``` OPTIONAL, default=30 : Discard the match or mismatch if the base on a read has less than MIN base quality. Important parameters. The lower this limit is, the less the damage is apparent.

```--max_coverage_limit MAX (DEFAULT 100)``` : If a position has equal or more than MAX reads (R1 or R2), the position is not used to calculate the damage. This option is to avoid high coverage regions of the genome being the main driver for the damage estimation program. Increase it to very a high number if your coverage is high. 

```--min_coverage_limit MIN (DEFAULT 1)``` : If a position has equal or less than MIN reads (R1 or R2), the position is not used to calculate the damage. This option is to calculate damage only in on-target regions (in cases of enrichment protocol such as exome ....).



OUTPUT : 
The output of the estimate_damage.pl program in out_file consist of 6 columns : [1] raw count of variant type [2] variant type (ex. G_T, G to T) [3] id (from the --id option) [4] frequency of variant [5] family (the variant type and reverse complement) [6] GIV-score  . If you have followthe standard protocol for acoustic shearing during library preparation you should obtain a GIV score for G_T around 2. You can directly feed this file to plot_damage.R or concatenates with other files and plot the concatenated file with plot_damage.R  : 

Rscript --vanilla plot_damage.R out_file (from the estimate_damage.pl program) name_of_figure.png``` OR ```Rscript --vanilla plot_damage.R out_file (from the estimate_damage.pl program) name_of_figure.pdf


### 3. random_sampling_and_estimate_damage.pl :

EXAMPLE :

DESCRIPTION :

OPTIONS :




### 4. estimate_damage_location.pl :

EXAMPLE :
```perl estimate_damage_location.pl  --mpileup1 file1.mpileup --mpileup2 file2.mpileup ---out file.damage --id idx12 --qualityscore 25 (DEFAULT 30)  --max_coverage_limit 200 --min_coverage_limit 10 ```

DESCRIPTION :

```estimate_damage_location.pl``` Is the main program to calculate the mutation rate for all the substitution including indels. The output of the program can be visualized using plot_damage_location.R. 

OPTIONS :

```--mpileup1``` REQUIRED mpileup file 1 (output of split_mapped_reads.pl)

```--mpileup2``` REQUIRED mpileup file 2 (output of split_mapped_reads.pl)

```--out``` REQUIRED Name of the output file

```--id``` REQUIRED id of the experiment. Can contain any letters or number. No space (ex: idx12) 

```--max_coverage_limit MAX (DEFAULT 100)``` : If a position has equal or more than MAX reads (R1 or R2), the position is not used to calculate the damage. This option is put in place in order to avoid high coverage regions of the genome being the main driver for the damage estimation program.

```--min_coverage_limit MIN (DEFAULT 1)``` : If a position has equal or less than MIN reads (R1 or R2), the position is not used to calculate the damage. This option is put in place in order to calculate damage only in on-target regions (in cases of enrichment protocol such as exome ....)

```--qualityscore MIN``` : Discard the match or mismatch if the base on a read has less than MIN base quality. Important parameters. The lower this limit is, the less the damage is apparent. 

OUTPUT :
The output of ```estimate_damage_location.pl``` is a table delimited file that can be directly used by ```plot_damage_location.R``` to visualized the damage function of the read positions. The columns are the following :
"id (from the --id option","Variant type","read (R1 or R2)","count (freq)","abs (absolute counts)","position on the read"

### 5. estimate_damage_location_context.pl :

EXAMPLE :
```perl estimate_damage_location_context.pl  --mpileup1 file1.mpileup --mpileup2 file2.mpileup ---out file.damage --id idx12 --qualityscore 25 (DEFAULT 30)  --max_coverage_limit 200 --min_coverage_limit 10 --context 1

OPTIONS :

Options are the same as for the estimate_damage_location.pl excpet for ```--context```

```--context``` correspond to the context of the damage. There are 3 possibilities : 
	[1] ```--context 1 ``` The damage is analysed function of the 5' nucleotide (C_[base], G_[base], T_[base] and A_[base])
    
    [2] ```--context 2 ``` The damage is analysed function of the 3' nucleotide ([base]_C, [base]_G, [base]_T and [base]_A)
    
    [3] ```--context 3 ``` The damage is analysed function of the 5' and 3' nucleotides (C_[base]_T, C_[base]_C, C_[base]_G, C_[base]_A, G_[base]_T ....)
    
With [base] being the variant analyzed.

OUTPUT :
The output of ```estimate_damage_location_context.pl``` is a table delimited file that can be directly used by ```plot_damage_location_context.R``` to visualized the damage function of the read positions within nucleotide context. The columns are the following :

"id (from the --id option)", "Variant type","read (R1 or R2)","counts (freq)","position on the read", "context", "absolute counts"
