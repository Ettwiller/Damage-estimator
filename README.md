### Damage-estimator

##REQUIREMENT : 
Prior to running the Damage estimator workflow please download and install the following programs :

**SAMTOOLS** (http://samtools.sourceforge.net/)


##OVERVIEW
This set of programs are designed to estimate the DNA damage when the DNA is sequenced using Illumina. 
The repository contains 1 basic programs :

 1. ```estimate_damage_location.pl``` 
 
 ##CONSIDERATIONS :
 Damage estimation is based on the systematics mutation rate difference between the first in pair and the second in pair reads. Therefore it is essential that the sequencing is done using Illumina in a paired end mode. BWA mapping is recommended and mapping in paired-end mode is required.  

