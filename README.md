# SNPVariantCalling
### Identify SNP Variant in Genomes of an organsim


The python script will be run from the command line according to the format below:

**python <script_name.py> [-v VCFFILE] [-g GFFFILE] [-f FASTAFILE]**

-v, -g and –f are used to specify the input files. The files can be provided in any order as long as the parameters are specified before each file. For the script to run, the three files (VCF file, Gff file and Fasta file) must be provided on the command line.

#### DESCRIPTION OF THE INPUT FILES REQUIRED AND THEIR FORMATS

The required files are the VCF file, GFF file and the Genome Fasta file.

-	The VCF file contains information on the variants chromosome, position of the SNP, reference allele, alternate allele, the quality value and others. The input VCF file must be in .gz zipped format, otherwise an error will be thrown and system will exit.

-	The GFF file contains annotation of the genome fasta file and include information such as sequence ID (similar to variant chromosome in VCF file), source, strand, start and end coordinates, score, feature type and other attributes. This file must be provided in .gff format.

-	The Fasta file is the file containing the nucleotide sequence of variant chromosomes. Each record in the genome fasta file contain the chromosome ID (same as variant chromosome ID), description and nucleotide sequence. This file must be provided in fasta format on the command line.


#### DESCRIPTION OF THE OUTPUT FILES AND AN EXPLANATION OF THEIR CONTENTS

The output of the script are three files namely output.tsv, variants_barplot.png and SNPvariantlog.txt. The output files are described below:

1. output.tsv is a tab-separated table with one row for each variant with QUAL > 20. There are nine columns in the table which are CHROM, POS, REF, ALT, TYPE, TRANSCRIPT, PROTEIN LOCATION, REF AA and ALT AA. The description of the columns are as follow:

•	CHROM is the chromosome ID of the variant.
•	POS is the position of variation in nucleotide sequence.
•	REF is the reference nucleotide at the position of the variant.
•	ALT is the alternate nucleotide in the variant
•	TYPE indicate whether a variant is non-coding (if the variant is not feature in CDS), synonymous (if it doesn’t result in amino acid change) or non-synonymous (if it results in change of amino acids).
•	TRANSCRIPT is the transcript id of the variant if it is found in coding region. If the variant is not in a coding region, NA is inserted as the value.
•	PROTEIN LOCATION is the coordinate of the amino acid from start of the protein.
•	REF AA is the reference amino acid at the location of the variant.
•	ALT AA is the alternate amino acid at the location of the variant.

2. A barplot named Variant_Barplot.png showing the proportion of the variants based on coding type (Non-coding, Synonymous and Non-synonymous).

3. A log file named SNPVariantlog.txt containing the filenames given at the command line, the count of the variants where quality is less than 20. Also if the script runs correctly, the location of the output files will be logged. Errors are captured should in case an invalid file format is given at the command line or any other error is encountered in the file and are reported in the log file.
