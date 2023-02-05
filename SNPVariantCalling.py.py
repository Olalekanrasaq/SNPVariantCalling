# import required libraries
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import gffutils
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gzip, math, logging, os
from Bio import SeqIO


# Parse command line arguments
# options -v, -g and -f are specified to tag the files as seperate files when provided in the CLI
# the input files are made required, the program will not run if any of the files is not provided 
parser = ArgumentParser(description='determine effect of SNP variant', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-v", "--vcfFile", required=True, help="Input VCF file")
parser.add_argument("-g", "--gffFile", required=True, help="Input gff file")
parser.add_argument("-f", "--fastaFile", required=True, help="Input genome fasta file")
args = vars(parser.parse_args())

# set up logging
# to log info about the scripts and error into a text file
logger = logging.getLogger()
logger.setLevel(logging.INFO)
script_log = logging.FileHandler('SNPvariantlog.txt')
script_log.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(script_log)

# setting variables for files provided in the command line
vcf = args['vcfFile']
gff = args['gffFile']
fasta = args['fastaFile']

# logging confirmation of the filenames given at the command line
logger.info('The input files are {}, {}, and {}\n\n'.format(vcf, gff, fasta))

# defining a function that takes in a vcf files (gz format) and parse it into a list of dictionaries
# each row of the vcf file is a dictionary
# if a vcf file is not a valid gzip file, an error is raised and log into the log file
def vcf_parser(vcf_file):
    try:
        with gzip.open(vcf_file, mode='rt') as VCF:
            vcf_list = []
            for line in VCF:
                ###Skip metadata lines
                if line[0] != '#':
                    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,unknown = line.rstrip().split('\t')
                    vcf_dict = {'CHROM': CHROM, 'POS': POS, 'ID': ID, 'REF': REF, 'ALT': ALT,
                               'QUAL': QUAL, 'FILTER': FILTER, 'INFO': INFO, 'FORMAT': FORMAT, 'unknown': unknown}
                    vcf_dict['INFO'] = {item.split('=')[0]:item.split('=')[1] for item in vcf_dict['INFO'].split(';')}
                    vcf_list.append(vcf_dict)
    except gzip.BadGzipFile:
        logger.error("File {} is not a gzipped file (b'##')\n\n".format(vcf_file))
        raise SystemExit(1)
    return vcf_list

# defining a function to get the non-coding variants in the vcf files
# non-coding variants are those variants whose the snp position are not found in the cds region 
def NonCodingVariants(vcf_file, cds_df):
    #the function takes a vcf file and a dataframe of variants found in cds region
    #create an empty list to store dictionaries for variants not found in the coding variants dataframe 
    ncd_list = []
    #get the chrom and pos of variants in cds region and store them in a list
    snps = [snp for snp in zip(cds_df['CHROM'], cds_df['POS'])]

    for variant in vcf_parser(vcf_file):
        if float(variant['QUAL']) > 20: # get only variants with quality greater than 20
            var = (variant['CHROM'], int(variant['POS'])) #get the chrom and pos of variants
            if var not in snps: # if the variant chrom and pos not found among variants in cds region
                # such variant is not in cds region, therefore a dictionary is created
                ncd_dic = {'CHROM': variant['CHROM'], 'POS': int(variant['POS']), 'REF': variant['REF'],
                            'ALT': variant['ALT'], 'Type': 'Non-coding', 'Transcript': 'NA',
                            'Protein Location': 'NA', 'Ref AA': 'NA', 'Alt AA': 'NA'}
                ncd_list.append(ncd_dic) # append the dictionary to the list of variants not found in cds region
    # turn the lists to a dataframe
    ncd_df = pd.DataFrame(ncd_list, columns=['CHROM', 'POS', 'REF', 'ALT', 'Type', 'Transcript',
                                         'Protein Location', 'Ref AA', 'Alt AA'])
    return ncd_df

# defining a function to concatenate variants found in cds region and those not found in cds region
# the function takes two arguments which are the dataframe of the two set of variants
def OutputTSV(ncd_df, cds_df):
    tsv_df = pd.concat([cds_df, ncd_df], axis=0).sort_values('POS', ascending=True).reset_index(drop=True)
    return tsv_df  #the output is a dataframe which is a combination of the two set of variants

# create a database for the gff file using a library gffutils
if not os.path.isfile('PlasmoDB-54_Pfalciparum3D7.db'):
    db = gffutils.create_db(gff, dbfn='PlasmoDB-54_Pfalciparum3D7.db', force=True, keep_order=True,
                        merge_strategy='merge', sort_attribute_values=True)
else:
    db = gffutils.FeatureDB('PlasmoDB-54_Pfalciparum3D7.db', keep_order=True)

# parse the genome fasta file using seqIO
# records is a list of each fasta sequence found in the genome faster file
records = [record for record in SeqIO.parse(fasta, "fasta")]


cds_list = []  #create an empty list to store variants found in the cds region
count = 0  #generate the count of variants with quality less than 20
for entry in vcf_parser(vcf):
    if float(entry['QUAL']) <= 20: #if variant quality is less than 20, increment the count by 1
        count += 1
        continue  #skip to the next variant in the file
    else:
        # find the CDS feature containing the position of the variant from the file
        for feature in db.region(seqid=entry['CHROM'], featuretype='CDS'):
            if not (feature.start <= int(entry['POS']) <= feature.end):
                continue  #skip to the next CDS feature if the position is not in the current CDS feature
            else:
                 # for CDS feature that contain variant snp position, get the transcript that is the parent.
                for transcript in db.parents(feature.id, featuretype='mRNA'):
                    region = 0  #get the region of the snp pos from the start of the chrom CDS
                    sequence = '' #get the CDS sequence of the region
                    alt_sequence = ''  #get the CDS sequence of the region after replacing the ref nuc with alt nuc
                    if transcript.strand == '+':
                        # if the transcript is on the positive strand
                        # get all the CDS children ordered by the start coordinate
                        # one of these will be the CDS we obtained in line 112, containing our region
                        for cds in db.children(transcript.id, featuretype='CDS', order_by='start'):
                            # if the CDS we are looking at now is the same as the one containing our region
                            if cds == feature:
                                # find how far into the CDS our region is
                                # add these to the counter values we initialised above
                                region = region + (int(entry['POS']) - cds.start + 1)
                                
                            else:
                                # if the CDS is not the one containing our region, get the length
                                # add it to the counter values we initialised above
                                region = region + (cds.end - cds.start + 1)
                            #get the sequence and alternate sequence of the region 
                            for record in records:
                                if record.id == entry['CHROM']:
                                    sequence += record.seq[cds.start-1:cds.end]
                                    alt_seq = record.seq[:int(entry['POS'])-1] + entry['ALT'] + record.seq[int(entry['POS']):]
                                    alt_sequence += alt_seq[cds.start-1:cds.end]
            
                    # if the transcript is on the negative strand
                    # get all the CDS children reverse ordered by the start coordinate
                    # one of these will be the CDS we obtained in line 112, containing our region
                    elif feature.strand == '-':
                        for cds in db.children(transcript.id, featuretype='CDS', order_by='start', reverse=True):
                            if cds == feature:
                                # find how far into the CDS our region is
                                # add these to the counter values we initialised above
                                # we are  measuring from the 3' end of the CDS
                                region = region + (cds.end - int(entry['POS']) + 1)
                                
                            else:
                                # if the CDS is not the one containing our region, get the length
                                # add it to the counter values we initialised above
                                region = region + (cds.end - cds.start + 1)
                            for record in records:
                                if record.id == entry['CHROM']:
                                    sequence += record.seq[cds.start-1:cds.end]
                                    alt_seq = record.seq[:int(entry['POS'])-1] + entry['ALT'] + record.seq[int(entry['POS']):]
                                    alt_sequence += alt_seq[cds.start-1:cds.end]

                        # reverse sequence for negative strand
                        sequence = sequence.reverse_complement()
                        alt_sequence = alt_sequence.reverse_complement()
                    
                    transcript = transcript.id
                    # translate the sequence to get the protein sequence 
                    protein_seq = sequence.translate()
                    alt_protein_seq = alt_sequence.translate() 
                    # transform coordinates to protein coordinates by dividing by 3, then rounding up
                    protein_loc = math.ceil(region/3)
                    # get the ref and alt amino acid using the protein coordinates
                    ref_AA = protein_seq[protein_loc - 1]
                    alt_AA = alt_protein_seq[protein_loc - 1]
                    #if ref amino acid is the same as alt amino acid, the snp type is synonymous
                    #otherwise is non-synonymous
                    if ref_AA == alt_AA:
                        type_snp = 'Synonymous'
                        alt_AA = 'NA'
                    else:
                        type_snp = 'Non-synonymous'
            
                    #create a dictionary for such variant found in the coding region
                    cds_dic = {'CHROM': entry['CHROM'], 'POS': int(entry['POS']), 'REF': entry['REF'],
                                'ALT': entry['ALT'], 'Type': type_snp, 'Transcript': transcript,
                                'Protein Location': protein_loc, 'Ref AA': ref_AA, 'Alt AA': alt_AA}  
                    #append the dictionary to the list
                    cds_list.append(cds_dic)
            break
# convert the list of dictionaries to a dataframe    
cds_df = pd.DataFrame(cds_list, columns=['CHROM', 'POS', 'REF', 'ALT', 'Type', 'Transcript',
                                         'Protein Location', 'Ref AA', 'Alt AA'])

# log the count of variants with quality less than 20 into the log file
logger.info('There are {} variants with quality less than 20 in the vcf file\n\n'.format(count))

# get the non-coding variants from the vcf file
# the function has been detailed earlier in line 59-78
ncd_df = NonCodingVariants(vcf, cds_df)

# generate the final dataframe containing both set of variants (coding and non-coding)
# the function is detailed earlier in line 82-84
tsv_df = OutputTSV(ncd_df, cds_df)

# generate the barplot for the non-coding, synonymous and non-synonymous variants  
sns.countplot(x='Type', data=tsv_df, color='blue')
plt.title('Proportion of the variants based on features')
plt.savefig("Variants_Barplot.png") # save the image as PNG

# output the tsv dataframe to a tsv file
tsv_df.to_csv('output.tsv', sep='\t', index=False)

# get the path to the current working directory
path = os.getcwd()

# all output files are saved in the current working directory
logger.info('Output files are saved to {}'.format(path))