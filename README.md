This script takes MaxQuant evidence.txt file with ADPr sites and the matching fasta file to extract localization scores and filter them to only retain high-quality sites.
It also annotates the site within the Uniprot ID, and transfers MBR intensities for modified peptides within groups if they were detected by MS/MS. 

File names with replicates should have the suffix "_01" (underscore and two digets), but otherwise identical names. This is used for form experiment groups. 
