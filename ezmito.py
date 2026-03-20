#!/usr/bin/env python

import argparse
import time
import sys
import os
from pathlib import Path
import pandas as pd
import shutil
import subprocess
import matplotlib.pyplot as plt
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Nexus import Nexus
import warnings
import numpy as np
import statistics
import pyfiglet

# Remove warnings
warnings.filterwarnings("ignore")


# Start tracking time at the beginning of the script
start_time = time.time()

# Genetic code descriptions
code_help = (
	"Specify the genetic code to use:\n"
	"  2  - Vertebrate Mitochondrial Code\n"
	"  3  - Yeast Mitochondrial Code\n"
	"  4  - Mold, Protozoan, and Coelenterate Mitochondrial Code\n"
	"  5  - Invertebrate Mitochondrial Code\n"
	"  9  - Echinoderm and Flatworm Mitochondrial Code\n"
	"  13 - Ascidian Mitochondrial Code\n"
	"  14 - Alternative Flatworm Mitochondrial Code\n"
	"  16 - Chlorophycean Mitochondrial Code\n"
	"  21 - Trematode Mitochondrial Code\n"
	"  22 - Scenedesmus obliquus Mitochondrial Code\n"
	"  23 - Thraustochytrium Mitochondrial Code\n"
	"  33 - Cephalodiscidae Mitochondrial UAA-Tyr Code"
)



names={	'ATP6': ['ATP6', 'A6','MT-ATP6', 'ATPASE6'] ,
		'ATP8': ['ATP8','A8', 'MT-ATP8', 'ATPASE8'] ,
		'COX1': ['COX1','MT-CO1', 'CO1',  'COXI', 'MTCO1', 'MTCOX1','MTCOXI'] ,	   
		'COX2': ['COX2','MT-CO2', 'CO2',  'COXII', 'MTCO2', 'MTCOX2','MTCOXII'] ,
		'COX3': ['COX3','MT-CO3', 'CO3',  'COXIII', 'MTCO3', 'MTCOX3','MTCOXIII'] ,
		'COB': ['COB','MT-CYB', 'CYTB', 'MTCYB'] ,
		'ND1': ['NAD1', 'MT-ND1', 'MTND1', 'NADH1', 'ND1'] ,
		'ND2': ['NAD2', 'MT-ND2', 'MTND2', 'NADH2', 'ND2' ] ,
		'ND3': ['NAD3', 'MT-ND3', 'MTND3', 'NADH3', 'ND3'] ,
		'ND4': ['NAD4', 'MT-ND4', 'MTND4', 'NADH4', 'ND4'] ,
		'NDL': ['NADL','MT-ND4L', 'MTND4L', 'NADH4L', 'ND4L', 'NDL', 'NAD4L', 'NADHL'] ,
		'ND5': ['NAD5', 'MT-ND5', 'MTND5', 'NADH5', 'ND5'] ,
		'ND6': ['NAD6', 'MT-ND6', 'MTND6', 'NADH6', 'ND6']
		}
		
		
ALIAS_TO_STANDARD = {
	alias.upper(): standard
	for standard, aliases in names.items()
	for alias in aliases
}

## General functions

# Concatenation
def concatenate(path, strand):

	path = Path(path)

	files = list(path.glob("*.nexus"))

	nexi = [(f.name, Nexus.Nexus(str(f))) for f in files]

	combined = Nexus.combine(nexi)

	output_file = path / f"{strand}.nexus"

	with open(output_file, "w") as f:
		combined.write_nexus_data(filename=f)

	return str(output_file)
	
	
# Translation
def translate(filename, genetic_code, outdir):
		processed_filename, filename_extension = os.path.splitext(filename)
		processed_filename = os.path.join(outdir,os.path.basename(processed_filename)+'_translated.fasta')
		
		with open(processed_filename, 'w') as outfile:
			for record in SeqIO.parse(filename, 'fasta'):
				record.seq = record.seq.translate(table=genetic_code, to_stop=False)
				SeqIO.write(record, outfile, 'fasta')
				
		return processed_filename
	
	
# NEXUS to FASTA
def nexus2fasta(nexus):
	output_fasta = nexus.replace('.nexus', '.fasta')
	
	# Create a list to store all the records
	records = []
	
	with open(nexus, "r") as nexus_file:
		# Parse the fasta file and annotate each record with molecule type 'DNA'
		for record in SeqIO.parse(nexus_file, "nexus"):
			record.annotations['molecule_type'] = 'DNA'
			record.seq = record.seq.replace('-','')
			records.append(record)
	
	# Write all the records to the Nexus file
	with open(output_fasta, "w") as fasta_file:
		SeqIO.write(records, fasta_file, "fasta")
#	os.remove(nexus)
	return output_fasta

# FASTA to NEXUS
def fasta2nexus(fasta, moltype):
	output_nexus = fasta.replace('.fasta', '.nexus')
	
	# Create a list to store all the records
	records = []
	
	with open(fasta, "r") as fasta_file:
		# Parse the fasta file and annotate each record with molecule type 'DNA'
		for record in SeqIO.parse(fasta_file, "fasta"):
			record.annotations['molecule_type'] = moltype
			records.append(record)
	
	# Write all the records to the Nexus file
	with open(output_nexus, "w") as nexus_file:
		SeqIO.write(records, nexus_file, "nexus")
	return output_nexus


# Pad the sequences
def pad_sequence(fasta):
	# Read all sequences from the given file
	sequences = list(SeqIO.parse(fasta, 'fasta'))  # Read all sequences at once
	
	# Determine the maximum sequence length
	max_length = max(len(seq.seq) for seq in sequences)
	
	padded_sequences = []
	for sequence in sequences:
		seq_length = len(sequence.seq)
		if seq_length < max_length:
			# Pad the sequence with gaps ('-') to match the maximum length
			sequence.seq = Seq(str(sequence.seq).ljust(max_length, '-'))
		padded_sequences.append(sequence)
	
	# Write all padded sequences to a new file
	processed_filename = fasta[:fasta.rfind('_')]
	processed_filename = processed_filename+'_aligned.fasta'
	with open(processed_filename, 'w') as outfile:
		SeqIO.write(padded_sequences, outfile, 'fasta')
	return processed_filename

# Remove stop codons
def isStopCodon(filename, geneticCode):

	processed_filename, filename_extension = os.path.splitext(filename)
	processed_filename = processed_filename[:processed_filename.rfind('_')]
	processed_filename = processed_filename+'_checked.fasta'
	
	with open(processed_filename, 'w') as outfile:
		for sequence in SeqIO.parse(filename, 'fasta'):
			translated=sequence.seq.translate(table=geneticCode, to_stop=False) #translating the nucleotide sequence to aa

			if len(sequence.seq) % 3 == 0: #check if the sequence is made of codon or there are more or less nucleotide
				if str(translated).count('*') == 1 and str(translated).find('*')+1 == len(translated) : #check if stop codon is at the end
					sequence.seq=sequence.seq[:-3]					
					print('Warning: Stop codon found at the end of the following sequence: ' + str(filename[filename.rfind('/')+1:])+' of '+ str(sequence.id) +'. It will be automatically removed')
				elif str(translated).count('*') > 1:
					print('***CODE ERROR: Multiple stop codons found in the following sequence: ' +str(sequence.id)+' of ' + str(filename)+ '. Check again the data set and re-submit it to the Web Server***\n')
					sys.exit()

				else:
					sequence.seq=sequence.seq


			elif len(sequence.seq) % 3 == 2:
				if str(translated).count('*') >= 1:
					print('***CODE ERROR: Stop codons found in the following truncated sequence: ' +str(sequence.id)+' of ' + str(filename)+ '. Check again the data set and re-submit it to the Web Server***')
					sys.exit()

				else:
					print('Warning: Found a truncated sequence: ' +str(sequence.id)+' of ' + str(filename)+ '. It will be automatically adjusted and the analysis is continuing...')
					sequence.seq=sequence.seq[:-2]

			elif len(sequence.seq) % 3 == 1:
				if str(translated).count('*') >= 1:
					print('***CODE ERROR: Stop codons found in the following truncated sequence: ' +str(sequence.id)+' of ' + str(filename)+ '. Check again the data set and re-submit it to the Web Server***\n')
					sys.exit()

				else:
					print('Warning: Found a truncated sequence: ' +str(sequence.id)+' of ' + str(filename)+ '. It will be automatically adjusted and the analysis is continuing...')
					sequence.seq=sequence.seq[:-1]
		
					
			SeqIO.write(sequence, outfile, "fasta")
	
		return processed_filename


# Replacing _verify alphabet function of Bio.Alphabet
def _verify_alphabet(sequence, alphabet):
	alphabet = set(alphabet) 
	return all(letter in alphabet for letter in sequence)

# Check IUPAC ambiguities
def isIUPAC(filename):
	for sequence in SeqIO.parse(filename, 'fasta'):
		seq = Seq(str(sequence.seq).upper()) #forcing to read the sequence using only the ambiguos DNA alphabet
		test = _verify_alphabet(seq, 'ATCGNRYBDKMHVSW') #testing it
		if test == False: #crashes if it founds a non recognized nucleotide
			print('***CODE ERROR: Found a not recognized nucleotide in the following file ' + str(filename) + ' of the following id' + str(sequence.id)+'. Check your matrix and re-submit the compressed folder to the web server.***\n')
			sys.exit()

# Check the lengths
def check_length(filename):
	seq_length=[]
	for seqrecord in SeqIO.parse(filename ,'fasta'):
		seq_length.append(len(seqrecord.seq))
	avg=statistics.mean(seq_length)
	if len(seq_length) > 1:
		stdev=statistics.stdev(seq_length)
		treshold = avg+(stdev*2)
		for seqrecord in SeqIO.parse(filename ,'fasta'):
			if len(seqrecord.seq) > treshold or len(seqrecord.seq) < treshold:
				print('Warning: the length of ' + str(seqrecord.id)+ ' in ' + 
				os.path.basename(filename).replace('_noempty.fasta','') + ' differs more than three standard deviation')

# Check the FASTA file
def is_fasta(fasta):
	with open(fasta, "r") as handle:
		fasta_records = list(SeqIO.parse(handle, "fasta"))
		if not fasta_records:
			print("***CODE ERROR. Invalid or empty fasta file detected. Please re-submit a valid file.***")
			sys.exit("Error: No valid FASTA records found.")
		else:
			return fasta_records

# Check the GFF3 file
def is_gff3 (gff_file):
	from BCBio.GFF import GFFExaminer
	examiner = GFFExaminer()
	in_handle = open(gff_file)
	first = in_handle.readline()
	if '##gff-version 3' in first:
		return True
	else:
		print("***CODE ERROR. Invalid GFF3 file detected. Please re-submit a valid file. Check the GFF3 format at: http://www.ensembl.org/info/website/upload/gff3.html***")
		sys.exit("Error: No valid GFF3 file.")
	in_handle.close()

# Check if there are any repeated IDs in the fasta file
def is_duplicated(filename):
	seen_ids = set()

	for record in SeqIO.parse(filename, 'fasta'):
		if record.id in seen_ids:
			print(f"***CODE ERROR. Duplicated ID: {record.id} found in {filename}. Check your matrix and re-submit.***")
			sys.exit(f"Error: Duplicated ID '{record.id}' found.")

		seen_ids.add(record.id)

#	print("Quality check | No duplicated IDs found. Continuing the analysis...")

# Check if the input file is a valid fasta and process the sequences
def check_fasta(filename, outdir):

	fasta_counter = 0
	
	processed_filename, filename_extension = os.path.splitext(filename)
	processed_filename = os.path.join(outdir,os.path.basename(processed_filename)+'_noempty.fasta')

	fasta_records = is_fasta(filename)
	if fasta_records:
		# Writing the cleaned file
		with open(processed_filename, 'w') as cleaned_file:
			for record in fasta_records:
				record.description = record.description.replace(' ', '_')
				SeqIO.write(record, cleaned_file, 'fasta')
				fasta_counter += 1

#	if fasta_counter:
#		print("Quality check | Fasta file processed without errors. Continuing...")
	#look for duplicated ids
	is_duplicated(processed_filename)
	return processed_filename

# Remove gaps ('-') from sequences
def remove_gaps(filename):
	processed_filename, filename_extension = os.path.splitext(filename)
	processed_filename = processed_filename[:processed_filename.rfind('_')]
	processed_filename = processed_filename+'_degapped.fasta'

	with open(processed_filename, 'w') as outfile:
		for record in SeqIO.parse(filename, 'fasta'):
			# Check for gaps and remove them
			if '-' in str(record.seq):
				print(f"Warning: Found gap (-) in {record.id}. Removing gap...")
				record.seq = Seq(str(record.seq).replace("-", ""))
			# Write each record to the file
			SeqIO.write(record, outfile, 'fasta')

#	print("Quality check | Gaps removed from sequences. File is ready.")
	
	return processed_filename
	
def replace_dir(directory):
	# If the directory exists, remove it
	if os.path.exists(directory):
		shutil.rmtree(directory)
	
	# Create the new directory
	os.makedirs(directory, exist_ok=True)

# Placeholder function for each subcommand
def ez_circular_subcommand(args):
	from Bio.SeqRecord import SeqRecord
	import pybedtools
	import re

	fasta_file = args.input
	bed_file = args.bed
	output_fasta = os.path.join(args.outdir, 'output.fasta')
	output_bed = os.path.join(args.outdir, 'output.bed')
	gene_name = args.start
	linear = '' if args.feature == 'circular' else 'linear'
	
	
	print(pyfiglet.figlet_format("EZcircular")) 
	
	print(f"Running EZcircular with the following parameters:\ninput: {fasta_file}\nbed: {bed_file}\noutdir: {args.outdir}\nfeature: {linear}\nstarting gene: {gene_name}\n\n")

	# Remove parentheses from gene name if present
	if '(' in gene_name or ')' in gene_name:
		gene_name = re.sub(r'\([^)]*\)', '', gene_name)
	
	def make_circular(fasta_file, bed_file):
		# Modify the FASTA file
		with open(fasta_file + '_linear.fa', 'w') as circularized:
			for record in SeqIO.parse(fasta_file, 'fasta'):
				seq_len = len(record.seq)
				record.seq = record.seq + 'N' * 100
				SeqIO.write(record, circularized, 'fasta')
				break

		# Modify the BED file
		columns = ["chrom", "start", "end", "name", "score", "strand"]
		bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=columns)
		bed_df.loc[len(bed_df.index)] = [bed_df['chrom'][0], seq_len, seq_len + 100, 
			'### Artifact of circularization of EZcircularize ###', '', '']
		bed_df.to_csv(bed_file + '_linear.bed', index=False, header=None, sep='\t')
		return fasta_file + '_linear.fa', bed_file + '_linear.bed'
	
	def parse_bed(bed_file, gene_name):
		try:
			bedtool = pybedtools.BedTool(bed_file)
			interval = bedtool.filter(lambda x: gene_name.lower() in x.name.lower())[0]
			return int(interval.start), int(interval.end)
		except Exception as e:
			print(f"Error parsing BED file: {e}")
			return None

	def reorder_and_write_fasta(gene_coordinates, fasta_file, output_fasta, gene_name):
		# Ensure output directory exists
		replace_dir(os.path.dirname(output_fasta))
		
		if gene_coordinates is None:
			print(f"Error: Unable to extract coordinates for gene {gene_name} from BED.")
		else:
			gene_start, gene_end = gene_coordinates
			gene_record = SeqRecord(Seq(""), id="", description="")
			for record in SeqIO.parse(fasta_file, "fasta"):
				gene_record.seq += record.seq[gene_start:]  # Adjust coordinates to 0-based index
				gene_record.seq += record.seq[:gene_start]
				gene_record.id = record.id + '_arranged_from_' + gene_name
			with open(output_fasta, "w") as output_handle:
				SeqIO.write(gene_record, output_handle, "fasta")
		return gene_record.id


	
	def write_bed(bed_file, output_bed, gene_coordinates, gene_record_id):
		columns = ["chrom", "start", "end", "name", "score", "strand"]
		bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=columns)
		length_seq = len(next(SeqIO.parse(fasta_file, "fasta")).seq)
		bed_df_first = bed_df[bed_df["start"] >= gene_coordinates[0]].copy()
		bed_df_first["start"] -= gene_coordinates[0]
		bed_df_first["end"] -= gene_coordinates[0]
		bed_df_second = bed_df[bed_df["start"] < gene_coordinates[0]].copy()
		to_add = length_seq - int(bed_df['end'].iloc[-1]) + bed_df_first['end'].iloc[-1]
		bed_df_second["start"] += to_add
		bed_df_second["end"] += to_add
		bed_df_ordered = pd.concat([bed_df_first, bed_df_second])
		bed_df_ordered["chrom"] = gene_record_id
		bed_df_ordered.to_csv(output_bed, sep='\t', index=False, header=None)

	# Apply circularization if needed
	if linear == 'linear':
		fasta_file, bed_file = make_circular(fasta_file, bed_file)

	# Parse the BED file and reorder sequences
	gene_coordinates = parse_bed(bed_file, gene_name)
	gene_record_id = reorder_and_write_fasta(gene_coordinates, fasta_file, output_fasta, gene_name)
	write_bed(bed_file, output_bed, gene_coordinates, gene_record_id)

	# Cleanup temporary files if created
	if linear == 'linear':
		os.remove(fasta_file)
		os.remove(bed_file)

	warnings.resetwarnings()

	

def ez_codon_subcommand(args):

	from cai2.CAI import RSCU
	from matplotlib.backends.backend_pdf import PdfPages
	import matplotlib.patches as patches


	def ezcodon_main(path, genetic_code, strand, outdir):
		for filename in os.listdir(path):
			filename = os.path.join(path,filename)
			processed_file = check_fasta(filename, outdir)   # Check if the fasta file is valid and free of duplicates
			check_length(processed_file)
			degapped_file = remove_gaps(processed_file)  # Check and remove gaps
			os.remove(processed_file)
			isIUPAC(degapped_file)
			clean_file = isStopCodon(degapped_file, genetic_code)
			os.remove(degapped_file)
			padded_file = pad_sequence(clean_file)
			os.remove(clean_file)
			nexus_file = fasta2nexus(padded_file, 'DNA')
		concatenated = concatenate(outdir, strand)
		return concatenated
		
	def CodonUsage(fasta, folder, geneticCode, strand):
		df_codon=pd.DataFrame(columns=['Species', 'Codon', 'RSCU']) #create and empty df
		df_aa=pd.DataFrame(columns=['Species', 'AA', 'Freq'])
	
		for sequence in SeqIO.parse(fasta, 'fasta'): #for each entry in the fasta combined sequences
			#RSCU
			whole=str('')
			d={}
			codons = [sequence.seq[i:i+3] for i in range(0, len(sequence.seq), 3)] #divide the sequence in triplets
			for codon in codons: #loop to parse single codons
				seq= Seq(str(codon).upper()).translate(table=geneticCode, to_stop=False)
				test = _verify_alphabet(seq, 'ACDEFGHIKLMNPQRSTVWY')
				if test == True:
					string=str(codon)
					whole += string
					d[str(codon).upper()]=sequence.id
				else:
					pass
					print(f'Warning: the codon ({codon}) in {sequence.id} produces an ambiguos aminoacid and it will be excluded from the analysis')
			seqlist=[]
			seqlist.append(whole)
			rscu=RSCU(seqlist, genetic_code=int(geneticCode))
			rscu=pd.DataFrame(rscu.items(), columns=['Codon', 'RSCU']) #calculates the RSCU for each seq
			df= pd.DataFrame(d.items(),columns=['Codon','Species']) #dictionary converted to a dataframe
			df= df.merge(rscu, how='left', on='Codon')
			df_codon = pd.concat([df_codon, df], ignore_index=True, sort=False) #add the current df to the empty one
			#AAfreq
			d_aa={'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0 ,'V':0, 'W':0, 'Y':0}
			translated=sequence.seq.translate(table=geneticCode, to_stop=False)
			for aa in translated:
				for k,v in d_aa.items():
					if aa == k:
						d_aa[k]=v+1
			for k,v in d_aa.items():
				d_aa[k]=v/len(translated)*100
			df= pd.DataFrame(d_aa.items(),columns=['AA','Freq'])
			df['Species']=str(sequence.id)
			df_aa = pd.concat([df_aa, df], ignore_index=True, sort=False)

		#RSCU

		df_codon=df_codon.reset_index() #reset index
		df_codon=df_codon[['Species','Codon','RSCU']] #discard the 1st column (it was the previous index)
		

		df_codon['AA']='' #add a new empty AA column
		for i in range(0,len(df_codon)): #for each triplet
			triplet=str(Seq(df_codon.iloc[i][1]).translate(table=geneticCode, to_stop=False)) #translate it to AA
			df_codon.iloc[i,3]=triplet #Add this result to the new column o
		
		output_RSCU = os.path.join(folder,strand+'_RSCU.csv')
		df_codon.to_csv(output_RSCU)


		#AAfreq
		df_aa=df_aa.reset_index()
		df_aa=df_aa[['Species','AA','Freq']]
		output_AAfreq =  os.path.join(folder,strand+'_AAfreq.csv')
		df_aa.to_csv(output_AAfreq)
		
		return output_RSCU, output_AAfreq
		
	def AA_freq_plot(table, strand, outdir):
	# Load data
		if os.path.exists(table):
			df = pd.read_csv(table).drop(columns=['Unnamed: 0'])
			nsp = df['Species'].nunique()

			# For amino acid frequency plot
			output_fig = os.path.join(outdir, f'{strand}_Aminoacid_frequency.pdf')
			
			if nsp <= 20:
				# Line plot for species with <= 20 unique species
				plt.figure(figsize=(20, 15))
				for species in df['Species'].unique():
					species_data = df[df['Species'] == species]
					plt.plot(species_data['AA'], species_data['Freq'], marker='o', label=species)
					
				plt.title(f'{strand} strand Amino Acid Frequency', fontsize=15)
				plt.xlabel('Amino Acids', fontsize=12)
				plt.ylabel('Frequency', fontsize=12)
				plt.legend(title="Species")
				plt.tight_layout()
				
				# Enable the grid
				plt.grid(True, which='both', axis='both', linestyle='-', linewidth=0.5, color='#cccccc')
				plt.gca().set_axisbelow(True)
				plt.savefig(output_fig, dpi=300)
				plt.close()  # Close the plot to avoid display

			else:
				# Box plot for species with > 20 unique species
				plt.figure(figsize=(20, 15))
				df.boxplot(column='Freq', by='AA', grid=False, showfliers=False)
				plt.title(f'{strand} strand Amino Acid Frequency', fontsize=15)
				plt.xlabel('Amino Acids', fontsize=12)
				plt.ylabel('Frequency (%)', fontsize=12)
				plt.suptitle('')  # Suppress the automatic 'by' title
				plt.tight_layout()
				
				# Enable the grid
				plt.grid(True, which='both', axis='both', linestyle='-', linewidth=0.5, color='#cccccc')
				plt.gca().set_axisbelow(True)
				plt.savefig(output_fig, dpi=300)
				plt.close()  # Close the plot to avoid display
				
	def RSCU_plot(table, strand, outdir):
		if os.path.exists(table):
			codon = pd.read_csv(table).drop(columns=['Unnamed: 0'])
			species_groups = codon.groupby('Species')

			output_fig = os.path.join(outdir, f'{strand}_strand_RSCU.pdf')

			# Create a PDF file to store all plots
			with PdfPages(output_fig) as pdf:

				for name, df in species_groups:
					df['Col'] = df.groupby('AA').cumcount()

					# Sort the dataframe by AA alphabetically
					df = df.sort_values(by='AA')

					# Create a figure with two subplots: one for bar plot and one for tile plot
					fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [8, 2]}, figsize=(15, 10))

					# Get unique AA and assign positions after sorting
					unique_AA = df['AA'].unique()
					x_positions = np.arange(len(unique_AA))

					# Create a dictionary to map AA to x-axis positions
					aa_to_position = {aa: i for i, aa in enumerate(unique_AA)}

					# Initialize bottom to 0 for stacking bars
					bottoms = {aa: 0 for aa in unique_AA}
					
					# Define the colormap with a consistent number of unique colors
					colors = plt.get_cmap('Dark2', len(df['Col'].unique()))

					# Plot the stacked bars for codons on ax1
					for i, row in df.iterrows():
						x = aa_to_position[row['AA']]  # Get x position from the AA
						ax1.bar(x, row['RSCU'], bottom=bottoms[row['AA']], 
								color=colors(row['Col']), edgecolor='black', width=0.8)
						
						# Update the bottom position for stacking
						bottoms[row['AA']] += row['RSCU']

					# Set x-axis labels as amino acids
					ax1.set_xticks(x_positions)
					ax1.set_xticklabels(unique_AA, position=(1,-.05))  # Assign amino acids (AA) as labels on the x-axis

					# Add titles and labels
					ax1.set_title(f"{name} {strand} strand RSCU")
					ax1.set_ylabel('RSCU')

					# Enable the grid on ax1
					ax1.grid(True, which='both', axis='y', linestyle='-', linewidth=0.5, color='#cccccc')
					ax1.grid(True, which='both', axis='x', linestyle='-', linewidth=0.5, color='#cccccc')

					# Ensure the grid is in the background
					ax1.set_axisbelow(True)

					# Adjust layout
					ax1.set_xlim(-0.5, len(df['AA'].unique()) - 0.5)

					# Tile plot (Codons) on ax2
					amino_acids = df['AA'].unique()

					# Reverse the order for y-axis and iteration for proper stacking
					col_values_reversed = df['Col'].max() - df['Col']

					for idx, row in df.iterrows():
						aa = row['AA']
						codon = row['Codon']
						col_value = col_values_reversed[idx]  # Use the reversed value
						aa_idx = list(amino_acids).index(aa)
						codon_idx = col_value  # Use the reversed index
						
						# Draw the rectangle tile for the codon
						ax2.add_patch(patches.Rectangle((aa_idx, codon_idx), 1, 1, 
														fc=colors(row['Col']), ec='black'))

						# Add codon label
						ax2.text(aa_idx + 0.5, codon_idx + 0.5, codon, ha='center', va='center', 
								color='white', fontweight='bold', size=7)

					# Set the axis limits and labels
					ax2.axes.get_yaxis().set_visible(False)
					ax2.axes.get_xaxis().set_visible(False)

					# Set axis limits
					ax2.set_xlim(0, len(amino_acids))
					ax2.set_ylim(0, df['Col'].max() + 1)

					# Hide axis spines and ticks for ax2
					ax2.spines[:].set_visible(False)
					ax2.tick_params(left=False, bottom=False)

					# Remove grid for ax2
					ax2.grid(False)

					# Save the figure to the PDF
					pdf.savefig(fig)
					plt.close()


	genetic_code = args.code
	outdir = args.outdir
	
	replace_dir(outdir)
	
	tmp = os.path.join(outdir, 'tmp')
	replace_dir(tmp)
	
	plots = os.path.join(outdir, 'plots')
	replace_dir(plots)
	
	tables = os.path.join(outdir, 'tables')
	replace_dir(tables)
	
	print(pyfiglet.figlet_format("EZcodon")) 
	
	if args.heavy and not args.light:
		J_path = args.heavy
		J_tmp = os.path.join(tmp, 'J')
		replace_dir(J_tmp)
		print(f"Running EZcodon with the following parameters:\nheavy chain fasta files: {J_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	if args.light and not args.heavy:
		N_path = args.light
		N_tmp = os.path.join(tmp, 'N')
		replace_dir(N_tmp)
		print(f"Running EZcodon with the following parameters:\nlight chain fasta files: {N_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	if args.light and args.heavy:
		J_path = args.heavy
		N_path = args.light
		J_tmp = os.path.join(tmp, 'J')
		replace_dir(J_tmp)
		N_tmp = os.path.join(tmp, 'N')
		replace_dir(N_tmp)
		JN_tmp = os.path.join(tmp, 'JN')
		replace_dir(JN_tmp)
		print(f"Running EZcodon with the following parameters:\nheavy chain fasta files: {J_path}\nlight chain fasta files: {N_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	# Ensure at least one of --heavy or --light is provided
	if not args.heavy and not args.light:
		parser.error("At least one of --heavy (-J) or --light (-N) must be provided.")
	
	
	
	
	if args.heavy and not args.light:
		concatenated = ezcodon_main(J_path, genetic_code, 'J', J_tmp)
		final_file = nexus2fasta(concatenated)
		output_RSCU, output_AAfreq = CodonUsage(final_file, tables, genetic_code, 'J')
		AA_freq_plot(output_AAfreq, 'J', plots)
		RSCU_plot(output_RSCU, 'J', plots)
		
	if args.light and not args.heavy:
		concatenated = ezcodon_main(N_path, genetic_code, 'N', N_tmp)
		final_file = nexus2fasta(concatenated)
		output_RSCU, output_AAfreq = CodonUsage(final_file, tables, genetic_code, 'N')
		AA_freq_plot(output_AAfreq, 'N', plots)
		RSCU_plot(output_RSCU, 'N', plots)
		
	if args.light and args.heavy:
	
		concatenatedJ = ezcodon_main(J_path, genetic_code, 'J', J_tmp)
		concatenatedN = ezcodon_main(N_path, genetic_code, 'N', N_tmp)
		
		# copy the files to the JN folder
		src_files = os.listdir(J_tmp)
		for file_name in src_files:
			full_file_name = os.path.join(J_tmp, file_name)
			if os.path.isfile(full_file_name) and full_file_name.endswith('.nexus') and file_name != 'J.nexus':
				shutil.copy(full_file_name, JN_tmp)
		src_files = os.listdir(N_tmp)
		for file_name in src_files:
			full_file_name = os.path.join(N_tmp, file_name)
			if os.path.isfile(full_file_name) and full_file_name.endswith('.nexus') and file_name != 'N.nexus':
				shutil.copy(full_file_name, JN_tmp)
		concatenated = concatenate(JN_tmp, 'JN')
		final_file = nexus2fasta(concatenated)
		output_RSCU, output_AAfreq = CodonUsage(final_file, tables, genetic_code, 'JN')
		AA_freq_plot(output_AAfreq, 'JN', plots)
		RSCU_plot(output_RSCU, 'JN', plots)
		
	shutil.rmtree(tmp)

def ez_map_subcommand(args):
	from pycirclize.parser import Gff
	
	gff_file = args.gff
	colorJ = args.colorJ
	colorN = args.colorN
	feature = args.feature
	outdir = args.outdir
	
	replace_dir(outdir)
	
	print(pyfiglet.figlet_format("EZmap")) 
	
	print(f"Running EZmap with the following parameters:\nmt_feature: {feature}\nGFF file: {gff_file}\ncolors: {colorJ}, {colorN}\noutdir: {outdir}\n\n")
	
	
	LABEL_QUALIFIERS = ["product", "Product", "name", "Name", "gene", "Gene"]
	
	def get_label(feat):
		"""Get label from feature qualifiers using priority list"""
		return next((feat.qualifiers.get(q, [""])[0] for q in LABEL_QUALIFIERS if feat.qualifiers.get(q, [""])[0]), "")

	def extract_feats_ci(gff, types, strand=None):
		"""Case-insensitive feature extraction using gff.records"""
		types_lower = [t.lower() for t in types]
		feats = [
			rec.to_seq_feature()
			for rec in gff.records
			if rec.type.lower() in types_lower
			and (strand is None or rec.strand == strand)
		]
		return feats

	def is_control_region(rec):
		"""Check if a GffRecord is a control/AT-rich region"""
		t = rec.type.upper().replace(' ', '').replace('-', '').replace('_', '')
		return any(x in t for x in ['A+T', 'ATREGION', 'CONTROLREGION', 'CR', 'DLOOP', 'NCR', 'NONCODING', 'SEQUENCEFEATURE'])

	def get_at_feats(gff):
		"""Return SeqFeature list for AT/control region records"""
		return [rec.to_seq_feature() for rec in gff.records if is_control_region(rec)]
	
	def plot_linear(gff_file, colorJ, colorN, outdir):
				gff = Gff(gff_file)

				f_cds_feats = extract_feats_ci(gff, ["CDS", "tRNA", "rRNA", "gene"], strand=1)
				r_cds_feats = extract_feats_ci(gff, ["CDS", "tRNA", "rRNA", "gene"], strand=-1)
				AT_feat = get_at_feats(gff)

				plt.figure(figsize=(30, 5))
				plt.ylim(-3, 6)
				plt.xlim(0, gff.range_size + 1)

				plt.arrow(0, 0, gff.range_size + 1, 0, head_width=0, head_length=0, width=0.01, fc='black', ec='black', lw=1)

				cnt = 1
				for feat in f_cds_feats:
					start, end = int(str(feat.location.start)), int(str(feat.location.end))
					label = get_label(feat)
					plt.arrow(start, 0, end - start, 0, head_width=0.4, head_length=(end - start) * 0.25, width=0.25, length_includes_head=True, fc=colorJ, ec='black', lw=1)
					if len(label) > 10:
						plt.text((start + end) / 2, .5, label, ha='center', va='baseline', fontsize=10, rotation=90)
					else:
						if cnt % 2 == 0:
							plt.text((start + end) / 2, .5, label, ha='center', va='baseline', fontsize=10, rotation=90)
						else:
							plt.text((start + end) / 2, -2, label, ha='center', fontsize=10, rotation=90)
						cnt += 1

				for feat in r_cds_feats:
					end, start = int(str(feat.location.start)), int(str(feat.location.end))
					label = get_label(feat)
					plt.arrow(start, 0, end - start, 0, head_width=0.4, head_length=(start - end) * 0.25, width=0.25, length_includes_head=True, fc=colorN, ec='black', lw=1)
					if len(label) > 10:
						plt.text((start + end) / 2, .5, label, ha='center', va='baseline', fontsize=10, rotation=90)
					else:
						if cnt % 2 == 0:
							plt.text((start + end) / 2, .5, label, ha='center', va='baseline', fontsize=10, rotation=90)
						else:
							plt.text((start + end) / 2, -2, label, ha='center', fontsize=10, rotation=90)
						cnt += 1

				for feat in AT_feat:
					start, end = int(str(feat.location.start)), int(str(feat.location.end))
					label = 'A+T rich region / Control region'
					plt.arrow(start, 0, end - start, 0, head_width=0.4, head_length=(end - start) * 0.25, width=0.25, length_includes_head=True, fc=colorJ, ec='black', lw=1)
					if len(label) > 10:
						plt.text((start + end) / 2, .5, label, ha='center', va='baseline', fontsize=10, rotation=90)
					else:
						if cnt % 2 == 0:
							plt.text((start + end) / 2, .5, label, ha='center', va='baseline', fontsize=10, rotation=90)
						else:
							plt.text((start + end) / 2, -2, label, ha='center', fontsize=10, rotation=90)
						cnt += 1

				plt.gca().spines['top'].set_visible(False)
				plt.gca().spines['right'].set_visible(False)
				plt.gca().spines['left'].set_visible(False)
				plt.gca().yaxis.set_visible(False)
				plt.xticks(np.arange(0, gff.range_size + 1, 500), rotation=90)
				plt.xlabel('Genomic position (bp)')
				plt.savefig(os.path.join(outdir, 'mt_linear_output.pdf'), bbox_inches='tight')
				plt.close()

	def plot_circular(gff_file, colorJ, colorN, outdir):
		from pycirclize import Circos
		from pycirclize.parser import Gff

		gff = Gff(gff_file)
		circos = Circos(sectors={gff.name: gff.range_size})
		circos.text(gff.target_seqid, size=15)

		sector = circos.sectors[0]
		cds_track = sector.add_track((80, 100))
		cds_track.axis(fc="white", ec="none")

		region = extract_feats_ci(gff, ["region"])
		f_cds_feats = extract_feats_ci(gff, ["CDS", "tRNA", "rRNA", "gene"], strand=1)
		r_cds_feats = extract_feats_ci(gff, ["CDS", "tRNA", "rRNA", "gene"], strand=-1)
		AT_feat = get_at_feats(gff)

		cds_track.genomic_features(region, plotstyle="box", r_lim=(90, 91), fc="black", lw=.5)
		cds_track.genomic_features(f_cds_feats, plotstyle="arrow", r_lim=(80, 100), fc=colorJ, lw=1.0)
		cds_track.genomic_features(r_cds_feats, plotstyle="arrow", r_lim=(80, 100), fc=colorN, lw=1.0)
		cds_track.genomic_features(AT_feat, plotstyle="arrow", r_lim=(80, 100), fc=colorJ, lw=1.0)

		pos_list, labels = [], []
		all_feats = extract_feats_ci(gff, ["CDS", "tRNA", "rRNA", "gene"]) + AT_feat
		for feat in all_feats:
			start, end = int(str(feat.location.end)), int(str(feat.location.start))
			pos = (start + end) / 2
			t = feat.type.upper().replace(' ', '').replace('-', '').replace('_', '')
			if any(x in t for x in ['A+T', 'ATREGION', 'CONTROLREGION', 'CR', 'DLOOP', 'NCR', 'NONCODING', 'SEQUENCEFEATURE']):
				label = 'A+T rich region / Control region'
			else:
				label = get_label(feat)
			if label == "" or label.startswith("hypothetical"):
				continue
			pos_list.append(pos)
			labels.append(label)

		cds_track.xticks(
			pos_list, labels,
			label_orientation="vertical",
			show_bottom_line=True,
			label_size=6,
			line_kws=dict(ec="white"),
		)
		cds_track.xticks_by_interval(
			interval=500,
			outer=False,
			show_bottom_line=True,
			label_formatter=lambda v: f"{int(v)} bp",
			label_orientation="vertical",
			line_kws=dict(ec="black"),
		)

		fig = circos.plotfig()
		fig.savefig(f'{outdir}/circular_plot.pdf', bbox_inches='tight')
		plt.close(fig)
		
	is_gff3(gff_file)
	
	if feature == 'linear':
		plot_linear(gff_file, colorJ, colorN, outdir)
	
	else:
		plot_circular(gff_file, colorJ, colorN, outdir)
	

def ez_mix_subcommand(args):
	import matplotlib.cm as cm
	
	fasta_file = args.input
	length = args.length
	identity = args.identity*100
	outdir = args.outdir
	blastn = os.path.join(args.blastn, 'blastn')
	
	replace_dir(outdir)
	
	print(pyfiglet.figlet_format("EZmix")) 
	
	print(f"Running EZmix with the following parameters:\ninput: {fasta_file}\noutdir: {outdir}\nidentity: {identity}\nlength: {length}\npath to blastn: {blastn}\n\n")
	
	def run_blast(query, subject, outdir):
			output_file = os.path.join(outdir, "blout")
			command = f"{blastn} -task blastn -query {query} -subject {subject} -outfmt '6 qseqid qstart qend sseqid sstart send pident length evalue bitscore slen qlen' -evalue 0.5 > {output_file}"
			subprocess.run(command, shell=True)
			return output_file if os.path.exists(output_file) and os.path.getsize(output_file) > 0 else None

	# Function to parse BLAST output
	def parse_blast_output(blast_file, min_length, min_similarity):
		# Define column names based on the BLAST output format
		columns = ['qseqid', 'qstart', 'qend', 'sseqid', 'sstart', 'send', 'pc', 'length', 'evalue', 'bitscore', 'slen', 'qlen']
		
		# Read the BLAST output into a pandas DataFrame
		df = pd.read_csv(blast_file, sep="\t", header=None, names=columns)
		
		# Filter rows based on minimum similarity (pc) and length
		filtered_df = df[(df['pc'] >= min_similarity) & (df['length'] >= min_length)]
		
		# Select only the relevant columns (you can modify this based on what you need)
		hits_df = filtered_df[['qseqid','qstart', 'qend', 'sseqid' , 'sstart', 'send', 'pc']]

		# Convert the filtered DataFrame to a list of dictionaries
		hits = hits_df.to_dict(orient='records')
		
		return hits
		
	# Function to create a PDF plot
	def create_plot(names, df, max_length, output_pdf, min_similarity, min_length):
		plt.figure(figsize=(10,7))
		plt.title('EZmix output (min length: '+ str(min_length) +  'bp ; min percent identity: ' +  str(min_similarity) + '%)')
		plt.xlabel('Assembly, bp')
		plt.xlim(-100, max_length+100)
		plt.ylim(-1, len(names))
		
		n = len(df.qseqid.unique())
		cmap = cm.get_cmap('Accent')
		colors = [cmap(i) for i in range(n)]
		my_colors = {key:value for key, value in zip(df.qseqid.unique(), colors)}

		legend_handles = []
		for index, row in df.iterrows():
			plt.plot([0,row['length']],  [row['qseqid'],row['qseqid']], color='black')
			if row['qstart'] != '':
				plt.plot([row['qstart'],row['qend']],  [row['qseqid'],row['qseqid']], color='black', lw=5)
				plt.plot([row['sstart'],row['send']],  [row['sseqid'],row['sseqid']], color='black', lw=5)
				line, = plt.plot([(row['qstart'] + row['qend']) / 2, (row['sstart'] + row['send']) / 2], 
						 [row['qseqid'], row['sseqid']], lw=2, color=my_colors[row['qseqid']])
				legend_handles.append((row['qseqid'], line))

		# Deduplicate legend entries
		seen = {}
		for label, handle in legend_handles:
			if label not in seen:
				seen[label] = handle
		plt.legend(seen.values(), seen.keys(), title="Query sequences", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

		plt.xticks(np.arange(0, max_length+1, 500), rotation=75)
		plt.grid(True, which='both', axis='both', linestyle='--', linewidth=0.5, color='lightgray')
		plt.savefig(output_pdf, bbox_inches='tight')
		plt.close()			
		

	processed_file = check_fasta(fasta_file, outdir)   # Check if the fasta file is valid and free of duplicates
	checked_file = remove_gaps(processed_file)  # Check and remove gaps
	
	os.remove(processed_file)
	
	
	sequences = list(SeqIO.parse(checked_file, 'fasta'))
	seq_names = [record.id for record in sequences]  # Truncate names to 10 characters
	seq_lengths = [len(record.seq) for record in sequences]
	lengths_df = pd.DataFrame({name:length for name, length in zip(seq_names, seq_lengths)}.items(), columns=['qseqid','length'])
	lengths_df['qstart'] = ''
	max_length = max(seq_lengths)
	
	hits = []
	for i in range(len(sequences) - 1):
		for j in range(i + 1, len(sequences)):
			query_file = os.path.join(outdir, f"primo_{i}.fasta")
			subject_file = os.path.join(outdir, f"secondo_{j}.fasta")
			
			SeqIO.write(sequences[i], query_file, 'fasta')
			SeqIO.write(sequences[j], subject_file, 'fasta')
			
			blast_output = run_blast(query_file, subject_file, outdir)
			if blast_output is None:
				continue
			else:
				hits += parse_blast_output(blast_output, length, identity)
			
			os.remove(query_file)
			os.remove(subject_file)
			if os.path.exists(blast_output):
				os.remove(blast_output)
				
	#### Attenzione, non crea un df correttamente. ci sono sseqid che mancano nei qseqid. quindi aggiungili a mano se non li trovi (missing obj). poi c'è da capire come non far cadere il primo qseqid nello 0 dell'y axis 			
				
	hits_df = pd.DataFrame(hits)
	
	hits_df = pd.concat([hits_df, lengths_df], ignore_index=True)
	# Create the plot
	output_pdf = os.path.join(outdir, f"{os.path.basename(fasta_file)}_output.pdf")
	create_plot(seq_names, hits_df, max_length, output_pdf, identity, length)
	print(f"Plot saved to {output_pdf}")
	

def ez_pipe_subcommand(args):

	from Bio.SeqRecord import SeqRecord
	from itaxotools.pygblocks import compute_mask, trim_sequence, Options
	from Bio import  AlignIO
	import re
	
	# Function to write the fixed initial part of the config file
	def write_initial_cfg(path):
		with open(f'{path}/partition_finder.cfg', 'w') as cfg_file:
				cfg_file.write('## ALIGNMENT FILE ##\n')
				cfg_file.write('alignment = infile.phy;\n\n')
				cfg_file.write('## BRANCHLENGTHS: linked | unlinked ##\n')
				cfg_file.write('branchlengths = linked;\n\n')
				cfg_file.write('## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##\n')
				cfg_file.write('models = mrbayes;\n\n')
				cfg_file.write('# MODEL SELECCTION: AIC | AICc | BIC #\n')
				cfg_file.write('model_selection = aicc;\n\n')
				cfg_file.write('## DATA BLOCKS: see manual for how to define ##\n')
				cfg_file.write('[data_blocks]\n')

# Function to write the final part of the config file
	def write_final_cfg(path):
		with open(f'{path}/partition_finder.cfg', 'a') as cfg_file:
				cfg_file.write('\n## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##\n')
				cfg_file.write('[schemes]\n')
				cfg_file.write('search = greedy;\n')
	
	
	def write_charsets(nexus):

		charset_file, filename_extension = os.path.splitext(nexus)
		charset_file = charset_file + '.charsets'
		new_content = ""

		# Open the input file for reading
		with open(nexus, 'r') as f:
				# Read the content and split to find the relevant portion after 'begin sets;'
				charsets = f.read().split('begin sets;')[1]
				
				# Split the content by semicolon to process each part separately
				path = charsets.split(';')
				
				for p in path:
						if 'charpartition combined' not in p:
								# Extract the gene name from the path
								gene_name = p.split('/')[-1].split('.')[0]
								
								gene_name = gene_name.split('_checked')[0]
								
								# Define the regex pattern
								pattern = r"(charset\s+)([\S]+\.nexus)(\s+=\s+\d+-\d+)"
								
								# Replace the matched path while keeping the rest intact and prepend with a tab
								new_text = re.sub(pattern, fr'\tcharset {gene_name}.nexus\3', p)
								
								# Append the modified string to new_content
								new_content += new_text + ';'
						else:
								break

		# Write the new content to the output file
		with open(charset_file, 'w') as f_out:
				f_out.write('BEGIN ASSUMPTIONS;')  # Write the 'begin sets;' line back
				f_out.write(new_content)
				f_out.write('\nEND;')  # Write the 'end sets;' line back 
				
		return charset_file	
				
	# NEXUS to PHYLIP
	def nexus2phylip(nexus):
		output_phylip = nexus.replace('.nexus', '.phy')
		
		# Create a list to store all the records
		records = []
		
		with open(nexus, "r") as nexus_file:
			# Parse the fasta file and annotate each record with molecule type 'DNA'
			for record in SeqIO.parse(nexus_file, "nexus"):
				record.annotations['molecule_type'] = 'DNA'
				records.append(record)
		
		# Write all the records to the Nexus file
		with open(output_phylip, "w") as phylip_file:
			SeqIO.write(records, phylip_file, "phylip-relaxed")
		os.remove(nexus)
		return output_phylip

	def remove_third_codon_position(filename, outdir):
		
		processed_filename, filename_extension = os.path.splitext(filename)
		processed_filename = os.path.join(outdir,os.path.basename(processed_filename)+'_twopos.fasta')
	
		with open(processed_filename, "w") as output_handle:
			for record in SeqIO.parse(filename, "fasta"):
				sequence = str(record.seq)
				trimmed_seq = ''.join([sequence[i:i+2] for i in range(0, len(sequence), 3)])
				trimmed_record = SeqRecord(Seq(trimmed_seq), id=record.id)
				SeqIO.write(trimmed_record, output_handle, "fasta")
				
		return processed_filename
			

	def Gblocks(filename_aa, filename_nt, outdir):
		
		
#		command = f'Gblocks {filename} -t=c -b3=4 -e=gblocks'
#		
#		subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#		
#		output_file = filename+'-gb'
#		htm_file = filename+'-gb.htm'
#		
#		renamed_file, filename_extension = os.path.splitext(filename)
#		renamed_file = os.path.join(outdir,os.path.basename(renamed_file)+'_trimmed.fasta')
#		
#		shutil.copy(output_file, renamed_file)
#		os.remove(htm_file)
#		os.remove(output_file)

		trimmed_file, filename_extension = os.path.splitext(filename_nt)
		trimmed_file = os.path.join(outdir,os.path.basename(trimmed_file)+'_trimmed.fasta')
		
		sequences = AlignIO.read(filename_aa, 'fasta')
		
		options = Options(
				IS = (len(sequences)*0.5)+1,
				FS = len(sequences)*0.85,
				CP = 4,
				BL1 = 10,
				BL2 = 10,
				GT = 0,
				GC = '-'
			)
			
		mask = compute_mask(sequences, options, log=False)
		mask = ''.join([char * 3 for char in mask])
		
		with open(trimmed_file, 'w') as outfile:
			for record in SeqIO.parse(filename_nt, 'fasta'):
				record.seq = Seq(trim_sequence(Seq(record.seq), mask))
				SeqIO.write(record, outfile, 'fasta')
		
		return trimmed_file
		

	def alignment(aa_filename, nt_filename, outdir): 
	
		tmp_filename, filename_extension = os.path.splitext(aa_filename)
		tmp_filename = os.path.join(outdir,os.path.basename(tmp_filename)+'_aligned.fasta')
		
		processed_filename, filename_extension = os.path.splitext(aa_filename)
		processed_filename = os.path.join(outdir,os.path.basename(processed_filename)+'_ntaligned.fasta')
		
		command = f'mafft --quiet {aa_filename} > {tmp_filename}'
		subprocess.run(command, shell=True)
	
		with open(processed_filename,'w') as nt_alignment:
			seq_length=[]
			for al_sequence in SeqIO.parse(tmp_filename,'fasta'):
				seq_length.append(len(al_sequence.seq))
				for nt_sequence in SeqIO.parse(nt_filename, 'fasta'):
					if al_sequence.id == nt_sequence.id:
						seq='' #it is 'refilled' with adjacent codon
						c=0
						for x in al_sequence.seq:
								
							if x != '-':
								x = 3*int(c+1)-3
								codon=nt_sequence.seq[x:x+3]
								seq=seq+str(codon)
								c+=1
							else:
								gap='---'
								seq=seq+gap
									
						nt_alignment.write('>'+nt_sequence.id+'\n'+seq+'\n')
						break

		return processed_filename, tmp_filename

	def translate(filename, genetic_code, outdir):
		processed_filename, filename_extension = os.path.splitext(filename)
		processed_filename = os.path.join(outdir,os.path.basename(processed_filename)+'_translated.fasta')
		
		with open(processed_filename, 'w') as outfile:
			for record in SeqIO.parse(filename, 'fasta'):
				record.seq = record.seq.translate(table=genetic_code, to_stop=False)
				SeqIO.write(record, outfile, 'fasta')
				
		return processed_filename
		
	def is_empty_sequence(fasta_file):
		for record in SeqIO.parse(fasta_file, "fasta"):
			if len(record.seq) == 0 or str(record.seq).strip() == "":
				return True
		return False
	
	def ezpipe_main(path, genetic_code, outdir, positions):
		for filename in os.listdir(path):
			filename = os.path.join(path,filename)
			processed_file = check_fasta(filename, outdir)   # Check if the fasta file is valid and free of duplicates
			check_length(processed_file)
			degapped_file = remove_gaps(processed_file)  # Check and remove gaps
			os.remove(processed_file)
			isIUPAC(degapped_file)
			clean_file = isStopCodon(degapped_file, genetic_code)
			os.remove(degapped_file)
			translated_file = translate(clean_file, genetic_code, outdir)
			alignment_nt_file, alignment_aa_file = alignment(translated_file, clean_file, outdir)
			gblocked = Gblocks(alignment_aa_file, alignment_nt_file, outdir)
			os.remove(alignment_aa_file)
			os.remove(alignment_nt_file)
			os.remove(translated_file)
			os.remove(clean_file)
			if is_empty_sequence(gblocked):
				print(f"WARNING: file '{gblocked}' is empty because GBlocks trimmed more than expected. This gene will be skipped.")
			else:
				nexus_file = fasta2nexus(gblocked, 'DNA')
			os.remove(gblocked)
			
		concatenated = concatenate(outdir, 'infile')
		charsets = write_charsets(concatenated)
		nexus2phylip(concatenated)
				# Write the initial part of the config file
		write_initial_cfg(outdir)
		
		# Open the partitions file for reading and the config file for appending
		with open(charsets, 'r') as pf, open(f'{outdir}/partition_finder.cfg', 'a') as cfg_file:
				for line in pf:
					match = re.search(r"charset\s+(.*)\s+=\s+([0-9]+)-([0-9]+)", line)
					if match:
						part = match.group(1).split('.')[0]  # Partition name without extension
						from_pos = int(match.group(2))  # First position in the partition
						to_pos = int(match.group(3))		# Last position in the partition

						# Option for three codon positions
						if positions == 3:
							for pos in range(3):
								cfg_file.write(f"{part}_p{pos+1} = {from_pos+pos}-{to_pos}\\3;\n")

						# Option for two codon positions
						elif positions == 2:
							for pos in range(2):
								cfg_file.write(f"{part}_p{pos+1} = {from_pos+pos}-{to_pos}\\2;\n")

		# Write the final part of the config file
		write_final_cfg(outdir)
			

			
	genetic_code = args.code
	outdir = args.outdir
	positions = args.positions
	
	replace_dir(outdir)
	
	tmp = os.path.join(outdir, 'tmp')
	replace_dir(tmp)
	
	genes = args.input
	
	print(pyfiglet.figlet_format("EZpipe")) 
	
	print(f"Running EZpipe with the following parameters:\ngenes path {genes}\ngenetic code: {genetic_code}\npositions: {positions}\noutdir: {outdir}\n\n")
	
	ezpipe_main(genes, genetic_code, tmp, positions)
	shutil.copy(os.path.join(tmp,'infile.phy'), os.path.join(outdir,'infile.phy'))
	shutil.copy(os.path.join(tmp,'partition_finder.cfg'), os.path.join(outdir,'partition_finder.cfg'))
	shutil.rmtree(tmp)
	
	
	

def ez_skew_subcommand(args):
	
	from collections import Counter
	from Bio.SeqUtils import GC123

	def ezskew_main(path, genetic_code, strand, outdir):
		for filename in os.listdir(path):
			filename = os.path.join(path,filename)
			processed_file = check_fasta(filename, outdir)   # Check if the fasta file is valid and free of duplicates
			check_length(processed_file)
			degapped_file = remove_gaps(processed_file)  # Check and remove gaps
			os.remove(processed_file)
			isIUPAC(degapped_file)
			clean_file = isStopCodon(degapped_file, genetic_code)
			os.remove(degapped_file)
			padded_file = pad_sequence(clean_file)
			os.remove(clean_file)
			nexus_file = fasta2nexus(padded_file, 'DNA')
		concatenated = concatenate(outdir, strand)
		return concatenated
		
	def first_skew(fasta):
		ATdic = {}
		for record in SeqIO.parse(fasta, 'fasta'):
			gc_tot, gc_1, gc_2, gc_3 = GC123(record.seq)
			at_tot = 100 - gc_tot
			ATdic[record.id] = at_tot
		first_bias =  pd.DataFrame(ATdic.items(), columns=['Species','AT%'])
		return first_bias
	
	def second_skew(fasta, strand):
		ACdic = {}
		GTdic = {}

		for record in SeqIO.parse(fasta, 'fasta'):
			if strand == 'J':
				dic = dict(Counter(record.seq))
				AC = int(dic.get('A')) + int(dic.get('C'))
				ACperc = AC / len(record.seq) * 100
				ACdic[record.id] = ACperc

			
			if strand == 'N':
				dic = dict(Counter(record.seq))
				GT = int(dic.get('G')) + int(dic.get('T'))
				GTperc = GT / len(record.seq) * 100
				GTdic[record.id] = GTperc

			
			if strand == 'JN':
				dic = dict(Counter(record.seq))
				AC = int(dic.get('A')) + int(dic.get('C'))
				ACperc = AC / len(record.seq) * 100
				ACdic[record.id] = ACperc
				dic = dict(Counter(record.seq))
				GT = int(dic.get('G')) + int(dic.get('T'))
				GTperc = GT / len(record.seq) * 100
				GTdic[record.id] = GTperc
				

		Jbias =  pd.DataFrame(ACdic.items(), columns=['Species','AC%'])
		Nbias =  pd.DataFrame(GTdic.items(), columns=['Species','GT%'])  

		return Jbias, Nbias
		
	def skew_main(string, x, y): #calculates at skew
		i=0 # i is A or G
		j=0 # j is T or C
		for nt in string:
			if nt.upper() == x.upper():
				i+=1
			elif nt.upper() == y.upper():
				j+=1
		if i+j == 0:
			tot_skew = 0 #just in case a+t is zero 
		else:
			tot_skew=(i-j)/(i+j) #at skew
		return tot_skew
		
	def third_skew (fasta, strand): #strand = J or N
		diz={} #empty dic
		for element in SeqIO.parse(fasta, 'fasta'): #element in this case correspond to the taxa present in the dataset
			N=[] #list for the first nucleotide position
			NN=[] #list for the second nucleotide position
			NNN=[] #list for the third nucleotide position
			#get the codons
			if len(element.seq) % 3 == 0:
				codons = [element.seq[i:i+3] for i in range(0, len(element.seq), 3)]
				for codon in codons: #loop to parse single codons
					first=N.append(codon[0])
					second=NN.append(codon[1])
					third=NNN.append(codon[2])
					
							
			
			diz[element.id]=[N,NN,NNN] #this dictionary is completed by the taxa id and the corresponig nucleotides 1st, 2nd and 3rd position
		
		skew = pd.DataFrame.from_dict(diz) #dictionary converted to a dataframe
		skew=skew.transpose()
		skew=skew.reset_index()
		skew.columns=['Species','N','NN','NNN'] #new column names

		#new columns (empty for the moment)
		skew['FIRST_AT_SKEW_'+strand]=''
		skew['FIRST_CG_SKEW_'+strand]=''
		skew['SECOND_AT_SKEW_'+strand]=''
		skew['SECOND_CG_SKEW_'+strand]=''
		skew['THIRD_AT_SKEW_'+strand]=''
		skew['THIRD_CG_SKEW_'+strand]=''

	 	#remove extra characters from the df columns
		skew['NNN'] = skew['NNN'].astype(str).str.replace('[','',regex=False).str.replace(']','',regex=False).str.replace("'",'',regex=False).str.replace(', ','',regex=False)
		skew['N'] = skew['N'].astype(str).str.replace('[','',regex=False).str.replace(']','',regex=False).str.replace("'",'',regex=False).str.replace(', ','',regex=False)
		skew['NN'] = skew['NN'].astype(str).str.replace('[','',regex=False).str.replace(']','',regex=False).str.replace("'",'',regex=False).str.replace(', ','',regex=False)

		#AT-CG skew calculations
		length = len(skew)
		for index in range(0,int(length)): #iterate over index
			row=skew.iloc[index]
			first=row[1] #N col
			second=row[2] #NN col
			third=row[3] #NNN col

			#at and cg skew for each N NN or NNN col and associate this value to the corresponding col (eg third_at_skew col)
			at=skew_main(first, 'A', 'T')
			cg=skew_main(first, 'C', 'G')
			skew.at[index,'FIRST_AT_SKEW_'+strand]=at
			skew.at[index,'FIRST_CG_SKEW_'+strand]=cg
			
			at=skew_main(second, 'A', 'T')
			cg=skew_main(second, 'C', 'G')
			skew.at[index,'SECOND_AT_SKEW_'+strand]=at
			skew.at[index,'SECOND_CG_SKEW_'+strand]=cg
			
			at=skew_main(third, 'A', 'T')
			cg=skew_main(third, 'C', 'G')
			skew.at[index,'THIRD_AT_SKEW_'+strand]=at
			skew.at[index,'THIRD_CG_SKEW_'+strand]=cg
			
			
		skew=skew[['Species', 'FIRST_AT_SKEW_'+strand, 'FIRST_CG_SKEW_'+strand, 'SECOND_AT_SKEW_'+strand, 'SECOND_CG_SKEW_'+strand, 'THIRD_AT_SKEW_'+strand, 'THIRD_CG_SKEW_'+strand]]
		return skew


	def plot_codon_skew(df, x_col, y_col, strand_col, title, ax, legend, species_color_map):
		markers = {"J": "o", "N": "^"}  # Markers for different strands
		species_unique = df['Species'].unique()
		
	
		
		# Plotting each species with its respective color and strand with its marker
		for species in species_unique:
			for strand in markers.keys():
				subset = df[(df['Species'] == species) & (df[strand_col] == strand)]
				ax.scatter(subset[x_col], subset[y_col], 
						   label = species if strand == list(markers.keys())[0] else "", 
						   marker=markers.get(strand), 
						   facecolor=species_color_map[species],
						   edgecolor='black',
						   linewidth=0.8)


		# Add horizontal and vertical lines at zero
		ax.axhline(0, color='gray', linestyle='--')
		ax.axvline(0, color='gray', linestyle='--')
		
		# Set axis limits based on the max skew values
		ax.set_xlim(-max_first_AT, max_first_AT)
		ax.set_ylim(-max_first_CG, max_first_CG)
		
		# Axis labels and title
		ax.set_xlabel('AT skew')
		ax.set_ylabel('CG skew')
		ax.set_title(title)
		
		if legend:
			# Species legend
			handles, labels = ax.get_legend_handles_labels()
			by_label = dict(zip(labels, handles))
			species_legend = ax.legend(by_label.values(), by_label.keys(), title="Species", bbox_to_anchor=(1.05, 1), loc='upper left')
			ax.add_artist(species_legend)

			# Strand shape legend
			strand_handles = [
				Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markeredgecolor='black', markersize=8, label='Heavy (J)'),
				Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', markeredgecolor='black', markersize=8, label='Light (N)'),
			]
			ax.legend(handles=strand_handles, title="Strand", bbox_to_anchor=(1.05, 0), loc='lower left')

			
		
	genetic_code = args.code
	outdir = args.outdir
	
	replace_dir(outdir)
	
	tmp = os.path.join(outdir, 'tmp')
	replace_dir(tmp)
	
	plots = os.path.join(outdir, 'plots')
	replace_dir(plots)
	
	tables = os.path.join(outdir, 'tables')
	replace_dir(tables)
	
	print(pyfiglet.figlet_format("EZskew")) 
	
	if args.heavy and not args.light:
		J_path = args.heavy
		J_tmp = os.path.join(tmp, 'J')
		replace_dir(J_tmp)
		print(f"Running EZskew with the following parameters:\nheavy chain fasta files: {J_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	if args.light and not args.heavy:
		N_path = args.light
		N_tmp = os.path.join(tmp, 'N')
		replace_dir(N_tmp)
		print(f"Running EZskew with the following parameters:\nlight chain fasta files: {N_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	if args.light and args.heavy:
		J_path = args.heavy
		N_path = args.light
		J_tmp = os.path.join(tmp, 'J')
		replace_dir(J_tmp)
		N_tmp = os.path.join(tmp, 'N')
		replace_dir(N_tmp)
		JN_tmp = os.path.join(tmp, 'JN')
		replace_dir(JN_tmp)
		print(f"Running EZskew with the following parameters:\nheavy chain fasta files: {J_path}\nlight chain fasta files: {N_path}\ngenetic code: {genetic_code}\noutdir: {outdir}\n\n")
	# Ensure at least one of --heavy or --light is provided
	if not args.heavy and not args.light:
		parser.error("At least one of --heavy (-J) or --light (-N) must be provided.")
		
	
	
	if args.heavy and not args.light:
		concatenated = ezskew_main(J_path, genetic_code, 'J', J_tmp)
		final_file = nexus2fasta(concatenated)
		FSK = first_skew(final_file)
		SSK = second_skew(final_file,'J')[0]
		TSK = third_skew(final_file,'J')


	if args.light and not args.heavy:
		concatenated = ezskew_main(N_path, genetic_code, 'N', N_tmp)
		final_file = nexus2fasta(concatenated)
		FSK = first_skew(final_file)
		SSK = second_skew(final_file,'N')[1]
		TSK = third_skew(final_file,'N')

	if args.light and args.heavy:
		
		concatenatedJ = ezskew_main(J_path, genetic_code, 'J', J_tmp)
		concatenatedN = ezskew_main(N_path, genetic_code, 'N', N_tmp)
		
		# copy the files to the JN folder
		src_files = os.listdir(J_tmp)
		for file_name in src_files:
			full_file_name = os.path.join(J_tmp, file_name)
			if os.path.isfile(full_file_name) and full_file_name.endswith('.nexus') and file_name != 'J.nexus':
				shutil.copy(full_file_name, JN_tmp)
		src_files = os.listdir(N_tmp)
		for file_name in src_files:
			full_file_name = os.path.join(N_tmp, file_name)
			if os.path.isfile(full_file_name) and full_file_name.endswith('.nexus') and file_name != 'N.nexus':
				shutil.copy(full_file_name, JN_tmp)
				
				
		concatenated = concatenate(JN_tmp, 'JN')	
		
		final_file = nexus2fasta(concatenated)		
		final_fileJ = nexus2fasta(concatenatedJ)
		final_fileN = nexus2fasta(concatenatedN)
		
		FSK = first_skew(final_file)
		
		SSKJ = second_skew(final_fileJ,'J')[0]
		SSKN = second_skew(final_fileN,'N')[1]
		
		TSKJ = third_skew(final_fileJ,'J')
		TSKN = third_skew(final_fileN,'N')
		
		
		SSK = SSKJ.merge(SSKN, on = 'Species', how='left')
		TSK = TSKJ.merge(TSKN, on = 'Species', how='left')
		
		
	df_tmp = FSK.merge(SSK, on = 'Species', how='left')
	df_final = df_tmp.merge(TSK, on = 'Species', how='left')
	
	output_table = os.path.join(tables,'Final_table.csv')
	df_final.to_csv(output_table)
	
	# Taking the relevant columns based on gene
	if args.light and args.heavy:
		AT = df_final[['Species', 'AT%']]
		AC = df_final[['Species', 'AC%']]
		GT = df_final[['Species', 'GT%']]
		AT.columns = ['Species', 'value']
		AC.columns = ['Species', 'value']
		GT.columns = ['Species', 'value']
		AT.loc[:, 'strand'] = 'AT%'
		AC.loc[:, 'strand'] = 'AC% (J/heavy strand)'
		GT.loc[:, 'strand'] = 'GT% (N/light strand)'
		df_final2 = pd.concat([AT, AC, GT])
	if args.heavy and not args.light:
		AT = df_final[['Species', 'AT%']]
		AC = df_final[['Species', 'AC%']]
		AT.columns = ['Species', 'value']
		AC.columns = ['Species', 'value']
		AT.loc[:, 'strand'] = 'AT%'
		AC.loc[:, 'strand'] = 'AC% (J/heavy strand)'
		df_final2 = pd.concat([AT, AC])
	if args.light and not args.heavy:
		AT = df_final[['Species', 'AT%']]
		GT = df_final[['Species', 'GT%']]
		AT.columns = ['Species', 'value']
		GT.columns = ['Species', 'value']
		AT.loc[:, 'strand'] = 'AT%'
		GT.loc[:, 'strand'] = 'GT% (N/light strand)'
		df_final2 = pd.concat([AT, GT])
		
	

	# Plotting freqpoly plot
	if len(df_final) > 30:
		fig, ax = plt.subplots(figsize=(10, 6))

		strand_colors = {'AT%': '#e41a1c', 'AC% (J/heavy strand)': '#377eb8', 'GT% (N/light strand)': '#4daf4a'}

		for s in df_final2['strand'].unique():
			subset = df_final2[df_final2['strand'] == s]
			ax.hist(subset['value'], histtype='step', 
					stacked=True, fill=True, bins=10, label=s, color=strand_colors.get(s))

		ax.set_xlabel('%')
		ax.set_ylabel('Frequency')
		ax.legend(title='Bias percentages')

		if args.light and args.heavy == 'JN':
			ax.set_title("First and second bias frequencies: AT%, AC% and GT%")
		elif args.light == 'N':
			ax.set_title("First and second bias frequencies: AT% and GT%")
		elif args.heavy == 'J':
			ax.set_title("First and second bias frequencies: AT% and AC%")

		output_file = os.path.join(plots, 'First_and_second_bias_frequency.pdf')
		plt.savefig(output_file, dpi=500, bbox_inches='tight')
		 
		 
	#plotting barplots	 
	fig, ax = plt.subplots(figsize=(14, 8))

	# Get unique species and strands
	species = df_final2['Species'].unique()
	strands = df_final2['strand'].unique()
	n_strands = len(strands)

	# Set the width of each bar and create offsets
	bar_width = 0.2
	x = np.arange(len(species))

	strand_colors = {'AT%': '#e41a1c', 'AC% (J/heavy strand)': '#377eb8', 'GT% (N/light strand)': '#4daf4a'}


	# Plot bars for each strand with an offset to dodge them
	for i, s in enumerate(strands):
		subset = df_final2[df_final2['strand'] == s]
		ax.bar(x + i * bar_width, subset['value'], width=bar_width, label=s, color=strand_colors.get(s))

	# Adjust the x-ticks and labels
	ax.set_xticks(x + bar_width * (n_strands - 1) / 2)
	ax.set_xticklabels(species, rotation=90)

	# Set labels and title based on the strand
	if 'JN' in strands:
		ax.set_title("First and second bias: AT%, AC% and GT%")
	elif 'N' in strands:
		ax.set_title("First and second bias: AT% and GT%")
	elif 'J' in strands:
		ax.set_title("First and second bias: AT% and AC%")

	# Set axis labels
	ax.set_ylabel('%')
	ax.set_xlabel('Species')

	# Add legend
	ax.legend(labels=strands)
	
	output_file = os.path.join(plots, 'First_and_second_bias.pdf')
	plt.savefig(output_file, dpi=500, bbox_inches='tight')
	
	# Plotting the skews
	if args.light and args.heavy:
		J = df_final[['Species','FIRST_AT_SKEW_J', 'FIRST_CG_SKEW_J', 'SECOND_AT_SKEW_J', 'SECOND_CG_SKEW_J', 'THIRD_AT_SKEW_J', 'THIRD_CG_SKEW_J']]
		N = df_final[['Species','FIRST_AT_SKEW_N', 'FIRST_CG_SKEW_N', 'SECOND_AT_SKEW_N', 'SECOND_CG_SKEW_N', 'THIRD_AT_SKEW_N', 'THIRD_CG_SKEW_N']]
		J.columns = ['Species', 'FAT', 'FCG', 'SAT', 'SCG', 'TAT', 'TCG']
		N.columns = ['Species', 'FAT', 'FCG', 'SAT', 'SCG', 'TAT', 'TCG']
		J['strand'] = 'J'
		N['strand'] = 'N'
		JN = pd.concat([J, N])
	if args.heavy and not args.light:
		JN = df_final[['Species','FIRST_AT_SKEW_J', 'FIRST_CG_SKEW_J', 'SECOND_AT_SKEW_J', 'SECOND_CG_SKEW_J', 'THIRD_AT_SKEW_J', 'THIRD_CG_SKEW_J']]
		JN.columns = ['Species', 'FAT', 'FCG', 'SAT', 'SCG', 'TAT', 'TCG']
		JN['strand'] = 'J'
	if args.light and not args.heavy:
		JN = df_final[['Species','FIRST_AT_SKEW_N', 'FIRST_CG_SKEW_N', 'SECOND_AT_SKEW_N', 'SECOND_CG_SKEW_N', 'THIRD_AT_SKEW_N', 'THIRD_CG_SKEW_N']]
		JN.columns = ['Species', 'FAT', 'FCG', 'SAT', 'SCG', 'TAT', 'TCG']
		JN['strand'] = 'N'

	# Dividing data for codon positions
	first = JN[['Species', 'FAT', 'FCG', 'strand']]
	second = JN[['Species', 'SAT', 'SCG', 'strand']]
	third = JN[['Species', 'TAT', 'TCG', 'strand']]

	# Calculate max values for axis limits
	max_first_AT = max(abs(first['FAT']).max(), abs(second['SAT']).max(), abs(third['TAT']).max())
	max_first_CG = max(abs(first['FCG']).max(), abs(second['SCG']).max(), abs(third['TCG']).max())

	# Unique species list for color assignment
	n_colors = len(first['Species'].unique())
	species_unique = first['Species'].unique()
	if n_colors <= 10:
		cmap = plt.get_cmap('tab10')
		species_color_map = {species: cmap(i) for i, species in enumerate(species_unique)}
	elif n_colors <= 20:
		cmap = plt.get_cmap('tab20')
		species_color_map = {species: cmap(i) for i, species in enumerate(species_unique)}
	else:
		species_color_map = {
			species: (1.0, 1.0, 1.0, 1.0)
			for species in species_unique
		}

	fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))
		
	plot_codon_skew(first, 'FAT', 'FCG', 'strand', 'First codon position', axes[0], legend=False, species_color_map=species_color_map)
	plot_codon_skew(second, 'SAT', 'SCG', 'strand', 'Second codon position', axes[1], legend=False, species_color_map=species_color_map)
	if n_colors <= 20:
		plot_codon_skew(third, 'TAT', 'TCG', 'strand', 'Third codon position', axes[2], True, species_color_map)
	else:
		plot_codon_skew(third, 'TAT', 'TCG', 'strand', 'Third codon position', axes[2], False, species_color_map)

	# Show plots
	plt.tight_layout()

	output_file = os.path.join(plots, 'Third_bias.pdf')
	plt.savefig(output_file, dpi=500)
			
	shutil.rmtree(tmp)

def ez_split_subcommand(args):
	
	from BCBio.GFF import GFFExaminer
	from BCBio import GFF

	fasta_file = args.input
	gff_file = args.gff
	outdir = args.outdir

	replace_dir(outdir)
	
	print(pyfiglet.figlet_format("EZsplit")) 
	
	print(f"Running EZsplit with the following parameters:\nfasta file: {fasta_file}\nGFF file: {gff_file}\noutdir: {outdir}\n\n")
	
	is_gff3(gff_file)
	is_fasta(fasta_file)
	
	examiner = GFFExaminer()
	in_handle = open(gff_file)

	ids = [list(i)[0] for i in examiner.available_limits(in_handle).get('gff_id').keys()]
		
	list_of_genes = list(names.keys())
	d_found_genes = {}
	for rec in GFF.parse(gff_file):
		rec_id = rec.id
		genes_found = []
		for feature in rec.features:
			try:
				if feature.qualifiers.get('gene_biotype')[0] == 'protein_coding':
					start, end, strand = int(str(feature.location.start)), int(str(feature.location.end)), str(feature.location.strand)
					gene = feature.qualifiers.get('Name')[0].upper()
					for k, v in names.items():
						if gene in v:
							gene_syn = k
						else:
							pass
					for record in SeqIO.parse(fasta_file, 'fasta'):
						if rec_id == record.id:
							record.description = record.description
							record.seq = record.seq[start:end]
							if strand == '-1':
								record.seq = record.seq.reverse_complement()

							genes_found.append(gene_syn)
							with open(f'{outdir}/{gene_syn.lower()}.fasta', 'a') as outfile:
								outfile.write('>'+str(record.description)+'\n'+str(record.seq)+'\n')

			except:
				pass
			d_found_genes[rec_id] = genes_found

	with open(f'{outdir}/missing_genes.txt', 'w') as outfile:
		for k,v in d_found_genes.items():
			missing_genes = list(set(list_of_genes) - set(v))
			if len(missing_genes) > 0:
				outfile.write(f'{k}\t{missing_genes}\n')
				
				
				
				

def ez_trampo_subcommand(args):	
	import os
	import subprocess
	import shutil
	import re
	from collections import defaultdict, Counter
	from operator import itemgetter
	from Bio.Nexus import Nexus
	from Bio.SeqUtils import GC123
	from Bio import SeqIO
	import numpy as np
	import pandas as pd
	import plotly.express as px
	from cai2.CAI import RSCU
	import scipy.stats as spstats
	from sklearn.decomposition import PCA
	from sklearn.preprocessing import StandardScaler
	from statsmodels.stats.multitest import multipletests
	
	fasta_files  = args.path
	gene_order   = args.gene_order
	genetic_code = args.code
	model		= args.model
	sequence	 = args.sequence
	thmm_tables  = args.tables
	outdir	   = args.outdir
	threads	  = args.threads
	
	replace_dir(outdir)
	
	print(pyfiglet.figlet_format("EZtrampo"))
	
	print(
		f"Running EZtrampo with the following parameters:\n"
		f"FASTA directory:	{fasta_files}\n"
		f"Genetic code:	   {genetic_code}\n"
		f"Model organism:	 {model}\n"
		f"Gene order:		 {gene_order}\n"
		f"Custom sequence:	{sequence}\n"
		f"Custom TMHMM tables:{thmm_tables}\n"
		f"Threads:			{threads}\n"
		f"Output directory:   {outdir}\n\n"
	)
	
	
	def revise_user_sequences(file, outdir):

		revised_records = []

		for record in SeqIO.parse(file, "fasta"):

			stem = record.id.split("_")[0].upper()

			standard_name = ALIAS_TO_STANDARD.get(stem).lower()

			if standard_name is None:
				print(
					f"Invalid gene name '{stem}' in user reference file."
				)
				raise Exception

			# rename sequence id
			record.id = standard_name + '_user'
			record.description = ""

			revised_records.append(record)

		output_file = os.path.join(outdir, "user.fas")

		SeqIO.write(revised_records, output_file, "fasta")

		return output_file

	def create_partitions(table,domain,molecule):
		# subset the table depending on domain name (the table passed must be a subsetted one from the model organism table. drop the index!
		dataframe=table[table['Domain'] == domain].reset_index(drop=True)
		# make a dictionary using starts and ends
		start_end_dict=dataframe.set_index('Start').to_dict()['End']
		string=domain+' ='
		# what to trim will be filled with minimum and maximum starts and ends
		what_to_trim=[]
		for start,end in start_end_dict.items():
			if molecule == 'protein':
				string=string+' '+str(start)+'-'+str(end)
				what_to_trim.extend([start,end])
			else:
				string=string+' '+str((3*start)-2)+'-'+str(3*end)
				what_to_trim.extend([(3*start)-2,3*end])
		trimstart=min(what_to_trim)
		trimend=max(what_to_trim)
		string=string+';'
		# returns the string that must be written in nexus file and the start and end to be trimmed
		return string, trimstart, trimend

	def move_domains_and_save_nexus(alignment, table, out_nexus_aa, out_nexus_nt, organism):
	# pick up the gene name, it will be necessary to subset the table
		gene_name = os.path.basename(alignment).split('_')[0].upper()
		# read the model organism table as a dataframe
		df=pd.read_csv(table, sep='\t')
		# subset the table depending on gene name
		subsetted_table = df[df['Gene'] == gene_name]
		# create output files : nexus nt and aa
		with open(out_nexus_aa, 'a') as nexus_file_aa:
			with open(out_nexus_nt,'a') as nexus_file_nt:
				# parse the fasta sequences until arrive to the gene of interest (it's the model organism one - eg atp8_hsa)
				for sequence in SeqIO.parse(alignment, 'fasta'):
					if organism in sequence.id:
						for index,row in subsetted_table.iterrows():
							# pick up the domain name, start and end of each domain from the model organism table (remeber, it's reduced to the gene of interest)
							start=int(row['Start'])
							end=int(row['End'])
							domain=row['Domain']
							# slice the sequence depending on the start and end
							seq_slide=sequence.seq[start-1:end]
							# count the gaps
							gaps=seq_slide.count('-')
							# if gaps are present then shift the whole portions and creates the final table
							if gaps > 0:
								df.at[index,'Start']+=gaps
								df.at[index,'End']+=gaps

						#each obj created here will contain three variables: the string to be written in the nexus file, the start and end to be trimmed.
						MA_AA=create_partitions(subsetted_table,'MA','protein')
						TM_AA=create_partitions(subsetted_table,'TM','protein')
						IM_AA=create_partitions(subsetted_table,'IM','protein')

						MA_NT=create_partitions(subsetted_table,'MA','dna')
						TM_NT=create_partitions(subsetted_table,'TM','dna')
						IM_NT=create_partitions(subsetted_table,'IM','dna')
						
						# group trim starts and ends
						trimstart_aa=[MA_AA[1],TM_AA[1],IM_AA[1]]
						trimstart_aa=min(trimstart_aa)
						trimend_aa=[MA_AA[2],TM_AA[2],IM_AA[2]]
						trimend_aa=max(trimend_aa)

						trimstart_nt=[MA_NT[1],TM_NT[1],IM_NT[1]]
						trimstart_nt=min(trimstart_nt)
						trimend_nt=[MA_NT[2],TM_NT[2],IM_NT[2]]
						trimend_nt=max(trimend_nt)
						

						#write the partitions files
						nexus_file_aa.write('begin sets;\n'+
											'\tcharset '+ MA_AA[0]+'\n'+
											'\tcharset '+ TM_AA[0]+'\n'+
											'\tcharset '+ IM_AA[0]+'\n'+
											'end;')

						nexus_file_nt.write('begin sets;\n'+
											'\tcharset '+ MA_NT[0]+'\n'+
											'\tcharset '+ TM_NT[0]+'\n'+
											'\tcharset '+ IM_NT[0]+'\n'+
											'end;')
						break

				else:
					if organism != 'cel':
						raise Exception('No association with gene name found. Please consider to rename your file as described in the manual')
		return trimstart_aa, trimend_aa, trimstart_nt, trimend_nt


	def remove_gaps_and_modelorganism(nexfile, trimstart, trimend, taxon_to_exclude):
		# Read the alignment
		align = AlignIO.read(nexfile, 'nexus')
		length = align.get_alignment_length()

		# Get gap positions
		nexus = Nexus.Nexus(nexfile)
		gaps = nexus.gaponly()
		gaps.extend(range(1, int(trimstart)))
		gaps.extend(range(trimend, length + 1))

		# Print information about removed characters
		if '_aa_aligned' in nexfile:
			print('The original alignment of ' + os.path.basename(nexfile) + ' lost ' + str(len(gaps)) + '/' + str(length) + ' - ' +
				  str(round(100 * len(gaps) / length, 2)) + ('% of amino acid characters'))
		else:
			print('The original alignment of ' + os.path.basename(nexfile) + ' lost ' + str(len(gaps)) + '/' + str(length) + ' - ' +
				  str(round(100 * len(gaps) / length, 2)) + ('% of nucleotide characters'))

		try:
			# Write Nexus data excluding gaps and specified taxon
			nexus.write_nexus_data(nexfile, exclude=gaps, delete=[taxon_to_exclude])
		except Exception:
			# Write Nexus data excluding gaps only if taxon_to_exclude is not present
			nexus.write_nexus_data(nexfile, exclude=gaps)

		return nexfile


	def partitions(alignment, table, model, outdir):
		out_nexus_aa = alignment.replace('.fasta','.nexus')
		out_nexus_nt = out_nexus_aa.replace('aa', 'nt')
		trimstart_aa, trimend_aa, trimstart_nt, trimend_nt = move_domains_and_save_nexus(alignment, table, out_nexus_aa, out_nexus_nt, model)
		remove_gaps_and_modelorganism(out_nexus_aa, trimstart_aa, trimend_aa, model)
		remove_gaps_and_modelorganism(out_nexus_nt, trimstart_nt, trimend_nt, model)



	def transform_dataframe(df, gene):
		# put some fixed variables to parse the tmhmm tables
		GENE_COLUMN_NAME = 'Gene'
		START_COLUMN_NAME = 'Start'
		END_COLUMN_NAME = 'End'
		DOMAIN_COLUMN_NAME = 'Domain'
		COLUMNS_TO_DROP = ['Species', 'TMHMM', 'Range']
		# renames column names
		df[['Start', 'End']] = df['Range'].str.split(n=1, expand=True)
		df[DOMAIN_COLUMN_NAME] = df[DOMAIN_COLUMN_NAME].str.upper()
		df[DOMAIN_COLUMN_NAME] = df[DOMAIN_COLUMN_NAME].str.replace('TMHELIX', 'TM')
		df[DOMAIN_COLUMN_NAME] = df[DOMAIN_COLUMN_NAME].str.replace('INSIDE', 'MA')
		df[DOMAIN_COLUMN_NAME] = df[DOMAIN_COLUMN_NAME].str.replace('OUTSIDE', 'IM')
		df = df.drop(COLUMNS_TO_DROP, axis=1)
		df.insert(loc=0, column=GENE_COLUMN_NAME, value=gene)
		df = df[[GENE_COLUMN_NAME, START_COLUMN_NAME, END_COLUMN_NAME, DOMAIN_COLUMN_NAME]]
		return df

	# check THMM tables
	def check_tables(tables, outdir):
		all_tables = pd.DataFrame()

		for file in tables:
			try:
				#check tables names
				gene = [
					key
					for key, value in names.items()
					for v in value
					if os.path.basename(file).split('.')[0].upper() == v.upper()
				][0]

			except:
				print('***CODE ERROR: Name of TMHMM table is incorrect. Check the manual to know how to rename them')
				raise ValueError
			
			df = pd.read_csv(file, sep='\t', comment='#', names=['Species', 'TMHMM', 'Domain', 'Range', 'Start', 'End'])
			df = transform_dataframe(df, gene)
			all_tables = pd.concat([all_tables, df], ignore_index=True)
			
		out_table = os.path.join(outdir, 'user.tsv')
		all_tables.to_csv(out_table, sep='\t',  index=False, header=True)
		return out_table

	# change TRAMPO file names
	def revise_names(files, tmp):
	# revise file names
		
		new_names = []

		for file in files:
			if not os.path.isfile(file):
				continue

			stem = os.path.splitext(os.path.basename(file))[0].upper()
			standard_name = ALIAS_TO_STANDARD.get(stem)
			
			if standard_name is None:
				continue

			new_file_path = os.path.join(tmp, standard_name.lower() + '.fa')

			shutil.copy(file, new_file_path)
			new_names.append(new_file_path)
		
		return new_names


	# TRAMPO n° of files
	def num_of_files(fasta_list):
		if len(fasta_list) < 13:
			print('Warning:', len(fasta_list), 'files detected in your input FASTA')

		
	def alignment(aa_filename, outdir):
		aligned_file, _ = os.path.splitext(aa_filename)
		aligned_file = os.path.join(outdir, os.path.basename(aligned_file) + '_ali_aa.fasta')
		command = f'mafft --quiet --thread {threads} {aa_filename} > {aligned_file}'
		subprocess.run(command, shell=True)
		return aligned_file

	def alignment_ref(aa_filename, nt_filename, ref, outdir):
		aligned_file	   = aa_filename.replace('_ali_aa.fasta', '_ali_aa_ref.fasta')
		processed_filename = aa_filename.replace('_ali_aa.fasta', '_ali_nt_ref.fasta')
		command = f'mafft --quiet --add {ref} {aa_filename} > {aligned_file}'
		subprocess.run(command, shell=True)
		with open(processed_filename, 'w') as nt_alignment:
			for al_sequence in SeqIO.parse(aligned_file, 'fasta'):
				for nt_sequence in SeqIO.parse(nt_filename, 'fasta'):
					if al_sequence.id == nt_sequence.id:
						seq = ''
						c = 0
						for x in al_sequence.seq:
							if x != '-':
								x = 3 * int(c + 1) - 3
								codon = nt_sequence.seq[x:x + 3]
								seq += str(codon)
								c += 1
							else:
								seq += '---'
						nt_alignment.write('>' + nt_sequence.id + '\n' + seq + '\n')
						break
		return processed_filename, aligned_file
		
	
	# Check IUPAC protein ambiguities
	def isIUPAC_prot(filename):
		allowed_alphabet = 'ACDEFGHIKLMNPQRSTVWYBXZ'
		modified = False
		records = []

		for sequence in SeqIO.parse(filename, 'fasta'):
			seq = Seq(str(sequence.seq).upper())
			test = _verify_alphabet(seq, allowed_alphabet)

			if not test:
				print(f'Warning: Found a not recognized amino acid in file {filename}, sequence ID {sequence.id}. It will be switched to X to continue the analysis. Otherwise, check it and re-submit a new job.\n')
				# sanitize sequence
				cleaned_seq = "".join([c if _verify_alphabet(c, allowed_alphabet) else "X" for c in seq])
				sequence.seq = Seq(cleaned_seq)
				modified = True

			records.append(sequence)

		# If any sequence was modified, overwrite the file
		if modified:
			with open(filename, "w") as out:
				SeqIO.write(records, out, "fasta")
			return False

		return True

	# -----------------------------
	# 1. Extract old charsets
	# -----------------------------
	def extract_old_charsets(path):
		charsets_dict = {}
		for file in os.listdir(path):
			if file.endswith('_checked_translated_ali_nt_ref.nexus'):
				with open(os.path.join(path, file)) as f:
					for line in f:
						if "charset" in line:
							line = str(line.strip("charset"))
							line = str(line.replace("MA =","MA = [").replace("TM =","TM = [").replace("IM =","IM = [").replace(";","]").replace(" ",",").replace(",=,","=").replace(",MA=[,","_MA=[").replace(",TM=[,","_TM=[").replace(",IM=[,","_IM=["))
							clean_file = str(file.replace("_checked_translated_ali_nt_ref.nexus","").replace(" ",""))
							key = clean_file + line.split("=")[0].strip(",")
							value = line.split("=")[1].strip(",\n") if "=" in line else line.strip(",\n")
							charsets_dict[key] = value
		return charsets_dict

	# -----------------------------
	# 2. Extract sequence lengths
	# -----------------------------
	def extract_lengths(path):
		lengths_dict = {}
		for file in os.listdir(path):
			if file.endswith('_checked_translated_ali_nt_ref.nexus'):
				with open(os.path.join(path, file)) as f:
					for line in f:
						if "nchar=" in line:
							start = "nchar"
							end = ";"
							length = str(line.split(start)[1].split(end)[0]).strip().strip("=")
							clean_file = str(file.replace("_checked_translated_ali_nt_ref.nexus", "").replace(" ", ""))
							lengths_dict[clean_file] = int(length)
		return lengths_dict

	# -----------------------------
	# 3. Build new intervals
	# -----------------------------
	def build_new_intervals(lengths_dict, old_charsets_dict, outdir):
		genes = ("atp6","atp8","cob","cox1","cox2","cox3","nd1","nd2","nd3","nd4","nd5","nd6","ndl")
		regions = ("MA", "TM", "IM")
		new_dict = {}
		current_sum = 0
		for key in sorted(lengths_dict.keys()):
			new_dict[key] = current_sum
			current_sum += int(lengths_dict[key])
		outfile = os.path.join(outdir, '39p_gen_dom.nex')
		with open(outfile, 'w') as p39:
			p39.write('#nexus\n')
			p39.write('begin sets;\n')
			for gene in genes:
				for region in regions:
					key = f"{gene}_{region}"
					if key in old_charsets_dict:
						intervals_str = old_charsets_dict[key].replace('[','').replace(']','')
						intervals = intervals_str.split(',')
						offset = new_dict.get(gene, 0)
						p39.write(f"\tcharset {key} = ")
						for interval in intervals:
							left, right = interval.split('-')
							new_left = int(left) + offset
							new_right = int(right) + offset
							p39.write(f"{new_left}-{new_right} ")
						p39.write(";\n")
			p39.write('end;\n')
		return outfile

	# -----------------------------
	# 4. Add codon sets
	# -----------------------------
	def add_codonsets2charsets(p39):
		outfile = p39.replace('39p_gen_dom', '117p_gen_dom_cod')
		with open(outfile, 'w') as p117_file:
			p117_file.write('#nexus\n')
			p117_file.write('begin sets;\n')
			with open(p39) as infile:
				for line in infile:
					if '#nexus' not in line and 'begin sets;' not in line and 'end;' not in line:
						input_string = line.replace('\t','')
						nom, _ = input_string.strip().split("=")
						nom = re.sub(r' +', ' ', nom).replace("charset ", "").replace(" ", "")
						pattern = r'(\w+)-(\w+)'
						key_value_pairs = re.findall(pattern, input_string)
						my_dict = dict(key_value_pairs)
						p117_file.write(f"\tcharset {nom}_cod1 = ")
						for i in range(len(my_dict)):
							key = list(my_dict.keys())[i]
							value = list(my_dict.values())[i]
							p117_file.write(f"{key}-{value}\\3 ")
						p117_file.write(";\n")
						p117_file.write(f"\tcharset {nom}_cod2 = ")
						for i in range(len(my_dict)):
							key2 = int(list(my_dict.keys())[i]) + 1
							value = list(my_dict.values())[i]
							p117_file.write(f"{key2}-{value}\\3 ")
						p117_file.write(";\n")
						p117_file.write(f"\tcharset {nom}_cod3 = ")
						for i in range(len(my_dict)):
							key3 = int(list(my_dict.keys())[i]) + 2
							value = list(my_dict.values())[i]
							p117_file.write(f"{key3}-{value}\\3 ")
						p117_file.write(";\n")
			p117_file.write('end;\n')
		return outfile

	# -----------------------------
	# 5. Process domain partitions
	# -----------------------------
	def process_domain_rc1(p117):
		MA = []
		TM = []
		IM = []
		outfile = p117.replace('117p_gen_dom_cod', '3p_dom')
		with open(outfile, 'w') as p3:
			with open(p117, 'r') as infile:
				for line in infile:
					line = line.replace("\3", "")
					if 'charset' in line:
						if 'MA_cod1' in line:
							processed_line = line.split('=')[1].replace(';', '').replace('\\3', '').strip()
							MA.extend(processed_line.split())
						elif 'TM_cod1' in line:
							processed_line = line.split('=')[1].replace(';', '').replace('\\3', '').strip()
							TM.extend(processed_line.split())
						elif 'IM_cod1' in line:
							processed_line = line.split('=')[1].replace(';', '').replace('\\3', '').strip()
							IM.extend(processed_line.split())
			p3.write("#nexus\n")
			p3.write("begin sets;\n")
			p3.write(f"\tcharset MA = {' '.join(MA)};\n")
			p3.write(f"\tcharset TM = {' '.join(TM)};\n")
			p3.write(f"\tcharset IM = {' '.join(IM)};\n")
			p3.write("end;\n")
		return outfile

	# -----------------------------
	# 6. Merge by region and codon
	# -----------------------------
	def merge_new_charsets_regions_codons_rc1(p117):
		genes = ("atp6","atp8","cob","cox1","cox2","cox3","nd1","nd2","nd3","nd4","nd5","nd6","ndl")
		domains = ("MA", "TM", "IM")
		codons = (1, 2, 3)
		outfile = p117.replace('117p_gen_dom_cod', '9p_dom_cod')
		with open(outfile, 'w') as p9:
			p9.write('#nexus\n')
			p9.write('begin sets;\n')
			for domain in domains:
				for codon in codons:
					tag = f"_{domain}_cod{codon}"
					p9.write(f"\tcharset {domain}_cod{codon} =")
					with open(p117, 'r') as f:
						for line in f:
							line = line.replace('\t', '')
							if '#nexus' not in line and 'begin sets;' not in line and 'end;' not in line:
								if tag in line:
									for word in genes:
										if word in line:
											line = line.replace("\n","").replace(word,"").replace("=","").replace(";","").replace("charset","").replace("updated_","").replace("_cod1","").replace("_cod2","").replace("_cod3","").replace(f"_{domain}","")
											line = re.sub(r'\s+', ' ', line).strip()
											p9.write(f" {line}")
					p9.write(";\n")
			p9.write('end;\n')

	# -----------------------------
	# 7. Merge by gene
	# -----------------------------
	def merge_new_charsets_genes_rc1(p39):
		genes = ["atp6","atp8","cob","cox1","cox2","cox3","nd1","nd2","nd3","nd4","nd5","nd6","ndl"]
		name_lines = {name: [] for name in genes}
		with open(p39, 'r') as f:
			for line in f:
				if '#nexus' not in line and 'begin sets;' not in line and 'end;' not in line:
					line = line.replace('\t', '').strip()
					for name in genes:
						if name in line:
							line = line.replace(name,"").replace("\n","").replace("=","").replace(";","").replace("charset","").replace("updated_","").replace("_MA","").replace("_TM","").replace("_IM","")
							name_lines[name].append(line)
		outfile = p39.replace('39p_gen_dom', '13p_gen')
		with open(outfile, 'w') as p13:
			p13.write('#nexus\n')
			p13.write('begin sets;\n')
			for name, lines in name_lines.items():
				if lines:
					p13.write(f"\tcharset {name} =")
					for line in lines:
						line = re.sub(r'\s+', ' ', line)
						p13.write(line)
					p13.write(";\n")
			p13.write('end;\n')

	# -----------------------------
	# 8. Merge by strand + codon + gene order
	# -----------------------------
	def merge_new_charsets_codon_strand(gene_order, p, num):
		if num == 117:
			outfile = p.replace('117p_gen_dom_cod', '18p_str_dom_cod')
		if num == 39:
			outfile = p.replace('39p_gen_dom', '6p_str_dom')
		dict_go = {}
		with open(gene_order, 'r') as gene_order_file:
			for line in gene_order_file:
				if line.startswith('positive:'):
					line = line.replace('positive:', '')
					for gene in line.split(','):
						dict_go[gene.replace('\n', '').strip()] = 'pos'
				if line.startswith('negative:'):
					line = line.replace('negative:', '')
					for gene in line.split(','):
						dict_go[gene.replace('\n', '').strip()] = 'neg'
		dict_go = dict(sorted(dict_go.items(), key=lambda x: (x[1], x[0])))
		new_partition_dict = {}
		for gene, chain in dict_go.items():
			with open(p, 'r') as partition_file:
				for line in partition_file:
					if line.startswith('\t') and gene in line:
						line = line.replace('\n', '').replace('\t', '')
						charset_name, position = line.split('=')
						charset_name = charset_name.replace('charset ', '').strip()
						new_charset_name = chain + charset_name.replace(gene, '')
						position = position.replace(';', '').strip()
						if new_charset_name not in new_partition_dict:
							new_partition_dict[new_charset_name] = [position]
						else:
							new_partition_dict[new_charset_name].append(position)
		with open(outfile, 'w') as out:
			out.write('#nexus\n')
			out.write('begin sets;\n')
			for k, v in new_partition_dict.items():
				out.write(f'\tcharset {k} = {" ".join(v)};\n')
			out.write('end;\n')
		return outfile

	# -----------------------------
	# 9. Merge codon positions
	# -----------------------------
	def process_codon_rc1(p18):
		cod1 = []
		cod2 = []
		cod3 = []
		with open(p18, 'r') as infile:
			for line in infile:
				if 'charset' in line:
					if 'cod1' in line:
						processed_line = line.split('=')[1].replace(';', '').strip()
						cod1.extend(processed_line.split())
					elif 'cod2' in line:
						processed_line = line.split('=')[1].replace(';', '').strip()
						cod2.extend(processed_line.split())
					elif 'cod3' in line:
						processed_line = line.split('=')[1].replace(';', '').strip()
						cod3.extend(processed_line.split())
		outfile = p18.replace('18p_str_dom_cod', '3p_cod')
		with open(outfile, 'w') as p3:
			p3.write("#nexus\n")
			p3.write("begin sets;\n")
			p3.write(f"\tcharset cod1 = {' '.join(cod1)};\n")
			p3.write(f"\tcharset cod2 = {' '.join(cod2)};\n")
			p3.write(f"\tcharset cod3 = {' '.join(cod3)};\n")
			p3.write("end;\n")
		return outfile

	# -----------------------------
	# 10. Merge by strand and codon
	# -----------------------------
	def process_codon_strand_rc1(p18):
		pos_cod1 = []
		pos_cod2 = []
		pos_cod3 = []
		neg_cod1 = []
		neg_cod2 = []
		neg_cod3 = []
		with open(p18, 'r') as infile:
			for line in infile:
				if 'charset' in line:
					if 'pos_' in line and '_cod1' in line:
						processed_line = line.split('=')[1].replace(';', '').strip()
						pos_cod1.extend(processed_line.split())
					elif 'pos_' in line and '_cod2' in line:
						processed_line = line.split('=')[1].replace(';', '').strip()
						pos_cod2.extend(processed_line.split())
					elif 'pos_' in line and '_cod3' in line:
						processed_line = line.split('=')[1].replace(';', '').strip()
						pos_cod3.extend(processed_line.split())
					elif 'neg_' in line and '_cod1' in line:
						processed_line = line.split('=')[1].replace(';', '').strip()
						neg_cod1.extend(processed_line.split())
					elif 'neg_' in line and '_cod2' in line:
						processed_line = line.split('=')[1].replace(';', '').strip()
						neg_cod2.extend(processed_line.split())
					elif 'neg_' in line and '_cod3' in line:
						processed_line = line.split('=')[1].replace(';', '').strip()
						neg_cod3.extend(processed_line.split())
		outfile = p18.replace('18p_str_dom_cod', '6p_cod_str')
		with open(outfile, 'w') as p6:
			p6.write("#nexus\n")
			p6.write("begin sets;\n")
			p6.write(f"\tcharset pos_cod1 = {' '.join(pos_cod1)};\n")
			p6.write(f"\tcharset pos_cod2 = {' '.join(pos_cod2)};\n")
			p6.write(f"\tcharset pos_cod3 = {' '.join(pos_cod3)};\n")
			p6.write(f"\tcharset neg_cod1 = {' '.join(neg_cod1)};\n")
			p6.write(f"\tcharset neg_cod2 = {' '.join(neg_cod2)};\n")
			p6.write(f"\tcharset neg_cod3 = {' '.join(neg_cod3)};\n")
			p6.write("end;\n")

	# -----------------------------
	# 11. Merge by strand, TM vs MA+IM, and codon
	# -----------------------------
	def process_12part_MAIM_merged_rc1(p18):
		partitions = {
			'neg_TM_cod1': [], 'neg_TM_cod2': [], 'neg_TM_cod3': [],
			'pos_TM_cod1': [], 'pos_TM_cod2': [], 'pos_TM_cod3': [],
			'neg_NO_cod1': [], 'neg_NO_cod2': [], 'neg_NO_cod3': [],
			'pos_NO_cod1': [], 'pos_NO_cod2': [], 'pos_NO_cod3': []
		}
		with open(p18, 'r') as infile:
			for line in infile:
				if 'charset' in line:
					processed_line = line.split('=')[1].replace(';', '').strip()
					for strand in ('neg', 'pos'):
						if f'{strand}_' in line:
							for codon in ('cod1', 'cod2', 'cod3'):
								if f'_{codon}' in line:
									if '_TM_' in line:
										partitions[f'{strand}_TM_{codon}'].extend(processed_line.split())
									elif '_MA_' in line or '_IM' in line:
										partitions[f'{strand}_NO_{codon}'].extend(processed_line.split())
		outfile = p18.replace('18p_str_dom_cod', '12p_MAIM_merged')
		with open(outfile, 'w') as p12:
			p12.write("#nexus\n")
			p12.write("begin sets;\n")
			for name, positions in partitions.items():
				p12.write(f"\tcharset {name} = {' '.join(positions)};\n")
			p12.write("end;\n")

	def split_partitions(infile, outfile):
		aln = Nexus.Nexus()
		aln.read(infile)
		aln.write_nexus_data_partitions(filename=outfile, charpartition=aln.charsets)
		return outfile

	# ── Module-level constants ─────────────────────────────────────────────────
	ACCEPTED_FILES = [
		'nt_combined_regions_IM',
		'nt_combined_regions_MA',
		'nt_combined_regions_TM',
		'nt_combined_regions.nex',
		'nt_combined_strand_regions_neg_IM',
		'nt_combined_strand_regions_neg_MA',
		'nt_combined_strand_regions_neg_TM',
		'nt_combined_strand_regions_pos_IM',
		'nt_combined_strand_regions_pos_MA',
		'nt_combined_strand_regions_pos_TM',
	]
	AA_GROUPS = {
		'A': 'G1', 'G': 'G1', 'P': 'G1', 'S': 'G1', 'T': 'G1',
		'D': 'G2', 'E': 'G2', 'N': 'G2', 'Q': 'G2',
		'H': 'G3', 'K': 'G3', 'R': 'G3',
		'I': 'G4', 'L': 'G4', 'M': 'G4', 'V': 'G4',
		'F': 'G5', 'W': 'G5', 'Y': 'G5',
		'C': 'G6',
	}
	DOMAIN_CHARSETS = ['IM', 'MA', 'TM']
	CHAIN_CHARSETS  = ['pos_TM', 'pos_MA', 'pos_IM', 'neg_TM', 'neg_MA', 'neg_IM']
	DOMAIN_COLORS = {'IM': '#E69F00', 'MA': '#56B4E9', 'TM': '#009E73'}
	CHAIN_COLORS  = {
		'pos_TM': '#E8F086', 'pos_MA': '#6FDE6E', 'pos_IM': '#235FA4',
		'neg_TM': '#FF4242', 'neg_MA': '#A691AE', 'neg_IM': '#0A284B',
	}
	CODON_POSITIONS = ['First', 'Second', 'Third']

	# ── Helper ─────────────────────────────────────────────────────────────────
	def _safe_write_image(fig, path: str) -> None:
		try:
			fig.write_image(path)
		except Exception as e:
			print(f'Warning: could not write {path} ({type(e).__name__}: {e}). '
				  f'HTML version is still available.')

	def ATCG_skew(seq, n1: str, n2: str) -> float:
		a = seq.upper().count(n1)
		b = seq.upper().count(n2)
		try:
			return round(a / float(a + b), 4)
		except ZeroDivisionError:
			return 0.0

	# ── 1. Codon-usage analysis ────────────────────────────────────────────────
	def analyze_codon_usage(folder, output_RSCU, output_AA, output_1_SKEW, output_23_SKEW, code):
		df_AA	  = pd.DataFrame()
		df_RSCU	= pd.DataFrame()
		df_1_SKEW  = pd.DataFrame()
		df_23_SKEW = pd.DataFrame()
		for partition in os.listdir(folder):
			if partition not in ACCEPTED_FILES:
				continue
			charset = (
				'Combined' if 'nex' in partition
				else partition
					.replace('nt_combined_genes_regions_', '')
					.replace('nt_combined_strand_regions_', '')
					.replace('nt_combined_regions_', '')
			)
			for record in SeqIO.parse(os.path.join(folder, partition), 'nexus'):
				codons = [record.seq[i:i + 3] for i in range(0, len(record.seq), 3)]
				translated	  = record.seq.replace('?', '-').translate(code)
				aminoacid_count = Counter(translated)
				aminoacid_count.pop('-', None)
				tmp_df_aa			= pd.DataFrame(aminoacid_count.items(), columns=['AA', 'Count'])
				tmp_df_aa['Species'] = record.id
				tmp_df_aa['Charset'] = charset
				df_AA = pd.concat([df_AA, tmp_df_aa], ignore_index=True)
				whole	 = ''.join(str(c) for c in codons)
				rscu_dict = RSCU([whole], code)
				aa_rscu   = {
					codon: [str(Seq(codon).translate(table=code)), rscu_val]
					for codon, rscu_val in rscu_dict.items()
				}
				temp_df_RSCU			= pd.DataFrame(
					([k] + v for k, v in aa_rscu.items()),
					columns=['Codon', 'AA', 'RSCU'],
				)
				temp_df_RSCU['Species'] = record.id
				temp_df_RSCU['Charset'] = charset
				df_RSCU = pd.concat([df_RSCU, temp_df_RSCU], ignore_index=True)
				try:
					GC = GC123(record.seq)
				except ZeroDivisionError:
					print(
						f'Warning: unable to calculate the first skew for '
						f'{record.id} {charset}. Check {partition} for details.'
					)
					GC = (0, 0, 0, 0)
				gc_positions = {
					'All': GC[0], 'First': GC[1], 'Second': GC[2], 'Third': GC[3],
				}
				gc_rows	   = {pos: [round(v, 4), round(100 - v, 4)] for pos, v in gc_positions.items()}
				temp_df_first = pd.DataFrame(
					([k] + v for k, v in gc_rows.items()),
					columns=['CodonPosition', 'GC_frequency', 'AT_frequency'],
				)
				temp_df_first['Species'] = record.id
				temp_df_first['Charset'] = charset
				df_1_SKEW = pd.concat([df_1_SKEW, temp_df_first], ignore_index=True)
				seq_str	= str(record.seq)
				first_pos  = Seq(seq_str[0::3])
				second_pos = Seq(seq_str[1::3])
				third_pos  = Seq(seq_str[2::3])
				skew_rows = [
					[record.id, charset, 'All',	ATCG_skew(record.seq, 'G', 'C'), ATCG_skew(record.seq, 'A', 'T')],
					[record.id, charset, 'First',  ATCG_skew(first_pos,  'G', 'C'), ATCG_skew(first_pos,  'A', 'T')],
					[record.id, charset, 'Second', ATCG_skew(second_pos, 'G', 'C'), ATCG_skew(second_pos, 'A', 'T')],
					[record.id, charset, 'Third',  ATCG_skew(third_pos,  'G', 'C'), ATCG_skew(third_pos,  'A', 'T')],
				]
				temp_df_23 = pd.DataFrame(skew_rows, columns=['Species', 'Charset', 'Position', 'GCskew', 'ATskew'])
				df_23_SKEW = pd.concat([df_23_SKEW, temp_df_23], ignore_index=True)
		df_AA	  = df_AA[['Species', 'Charset', 'AA', 'Count']]
		df_RSCU	= df_RSCU[['Species', 'Charset', 'Codon', 'AA', 'RSCU']]
		df_1_SKEW  = df_1_SKEW[['Species', 'Charset', 'CodonPosition', 'GC_frequency', 'AT_frequency']]
		df_AA.to_csv(output_AA, index=False, sep='\t')
		df_RSCU.to_csv(output_RSCU, index=False, sep='\t')
		df_1_SKEW.to_csv(output_1_SKEW, index=False, sep='\t')
		df_23_SKEW.to_csv(output_23_SKEW, index=False, sep='\t')
		return df_AA, df_RSCU, df_1_SKEW, df_23_SKEW

	# ── 2. Amino acid frequency ────────────────────────────────────────────────
	def perform_hist(data, title, output_stats_folder):
		aa_groups_df			 = pd.DataFrame(AA_GROUPS.items(), columns=['AA', 'Groups'])
		new_data				 = data.merge(aa_groups_df, on='AA', how='left')
		new_data['Groups_by_Domain'] = new_data['Charset'] + '_' + new_data['Groups']
		new_data.dropna(subset=['Groups'], inplace=True)
		group_domains			= new_data['Groups_by_Domain'].unique()
		group_counts	  = new_data.groupby(['Charset', 'Groups'])['Count'].sum().unstack(fill_value=0)
		group_percentages = group_counts.div(group_counts.sum(axis=1), axis=0) * 100
		final_result	  = group_percentages.stack().reset_index()
		final_result.columns = ['Charset', 'Group', 'Percentage']
		split_idx = 1 if title == 'by domain' else 2
		results = []
		for i in range(len(group_domains)):
			for j in range(i + 1, len(group_domains)):
				parts_i = group_domains[i].split('_')
				parts_j = group_domains[j].split('_')
				if (
					len(parts_i) > split_idx
					and len(parts_j) > split_idx
					and parts_i[split_idx] == parts_j[split_idx]
					and parts_i[split_idx] in AA_GROUPS.values()
				):
					data_1 = new_data[new_data['Groups_by_Domain'] == group_domains[i]]['Count']
					data_2 = new_data[new_data['Groups_by_Domain'] == group_domains[j]]['Count']
					stat, p_val = spstats.mannwhitneyu(data_1, data_2)
					results.append({
						'Group_Domain_1':		 group_domains[i],
						'Group_Domain_2':		 group_domains[j],
						'MannWhitneyU_statistic': stat,
						'p_value':				p_val,
					})
		results_df					= pd.DataFrame(results)
		results_df['corrected_p_value'] = multipletests(results_df['p_value'], method='bonferroni')[1]
		mean_stdev		 = new_data.groupby('Groups_by_Domain').agg({'Count': ['mean', 'std']})
		mean_stdev.columns = ['mean_Freq', 'stdev_Freq']
		mean_stdev		 = mean_stdev.reset_index()
		modified_title = title.replace(' ', '')
		results_df.to_csv(
			os.path.join(output_stats_folder, f'AminoAcid_frequency_{modified_title}_wilcoxon-test.tsv'),
			sep='\t', index=False,
		)
		mean_stdev.to_csv(
			os.path.join(output_stats_folder, f'AminoAcid_frequency_{modified_title}_mean-stdev-test.tsv'),
			sep='\t', index=False,
		)
		fig = px.histogram(final_result, x='Charset', y='Percentage', color='Group')
		fig.update_xaxes(tickangle=0, title=None, ticks='outside', categoryorder='array')
		fig.update_yaxes(title='Relative Frequencies (%)')
		fig.update_traces(textposition='outside')
		fig.update_layout(
			{'plot_bgcolor': '#ffffff', 'paper_bgcolor': '#ffffff'},
			uniformtext_minsize=8, uniformtext_mode='show',
			autosize=False, width=2000, height=1000,
			title={'text': 'Amino acid frequency ' + title},
		)
		return fig

	# ── 3. GC frequency ────────────────────────────────────────────────────────
	def perform_boxplot(data, title, colors, output_stats_folder):
		fig = px.box(
			data, x='Charset', y='GC_frequency', color='Charset',
			facet_col='CodonPosition', facet_col_wrap=1,
			color_discrete_map=colors, template='simple_white',
		)
		fig.update_traces(width=0.5)
		fig.update_layout(height=1000, width=1000)
		fig.for_each_annotation(lambda a: a.update(text=a.text.split('=')[1] + ' codon positions'))
		if title == 'by chain':
			fig.update_xaxes(
				tickangle=0, title=None, ticks='outside', categoryorder='array',
				categoryarray=sorted(
					data['Charset'].unique(),
					key=lambda x: (x.split('_')[1], x.split('_')[0]),
				),
			)
		data = data.copy()
		data['TEST'] = data['Charset'] + '_' + data['CodonPosition']
		group_codons = data['TEST'].unique()
		results = []
		for i in range(len(group_codons)):
			for j in range(i + 1, len(group_codons)):
				parts_i = group_codons[i].split('_')
				parts_j = group_codons[j].split('_')
				if parts_i[-1] == parts_j[-1]:
					data_1 = data[data['TEST'] == group_codons[i]]['GC_frequency']
					data_2 = data[data['TEST'] == group_codons[j]]['GC_frequency']
					stat, p_val = spstats.mannwhitneyu(data_1, data_2)
					results.append({
						'Group_Domain_1':		 group_codons[i],
						'Group_Domain_2':		 group_codons[j],
						'MannWhitneyU_statistic': stat,
						'p_value':				p_val,
					})
		results_df					= pd.DataFrame(results)
		results_df['corrected_p_value'] = multipletests(results_df['p_value'], method='bonferroni')[1]
		mean_stdev		 = data.groupby('TEST').agg({'GC_frequency': ['mean', 'std']})
		mean_stdev.columns = ['mean_Freq', 'stdev_Freq']
		mean_stdev		 = mean_stdev.reset_index()
		modified_title = title.replace(' ', '')
		results_df.to_csv(
			os.path.join(output_stats_folder, f'GC_frequency_{modified_title}_wilcoxon-test.tsv'),
			sep='\t', index=False,
		)
		mean_stdev.to_csv(
			os.path.join(output_stats_folder, f'GC_frequency_{modified_title}_mean-stdev-test.tsv'),
			sep='\t', index=False,
		)
		return fig

	# ── 4. PCA and scatter plots ───────────────────────────────────────────────
	def scatter_plot(data, title, colors):
		data			 = data.copy()
		data['dummy_size'] = 1
		fig = px.scatter(
			data, x='ATskew', y='GCskew', color='Charset', symbol='Position',
			opacity=0.7, symbol_sequence=['circle', 'triangle-up', 'x'],
			color_discrete_map=colors, size='dummy_size', size_max=15,
		)
		fig.update_layout(
			legend=dict(x=1.02, y=1, bordercolor='white', borderwidth=1, font=dict(size=20)),
			title='AT-GC skews ' + title, template='plotly_white',
			autosize=False, width=2000, height=1000,
			xaxis=dict(tickfont=dict(size=30)),
			yaxis=dict(tickfont=dict(size=30)),
		)
		return fig

	def pca_plot(data, title, colors):
		features  = ['ATskew', 'GCskew']
		x_scaled  = StandardScaler().fit_transform(data[features].values)
		pca_model = PCA(n_components=2).fit(x_scaled)
		per_var   = np.round(pca_model.explained_variance_ratio_ * 100, decimals=1)
		coords	= pca_model.transform(x_scaled)
		principal_df = pd.DataFrame(coords, columns=['PCA1', 'PCA2'])
		final_df	 = pd.concat(
			[principal_df, data[['Charset', 'Position']].reset_index(drop=True)], axis=1
		)
		final_df['dummy_size'] = 1
		fig = px.scatter(
			final_df, x='PCA1', y='PCA2', color='Charset', symbol='Position',
			opacity=0.7, symbol_sequence=['circle', 'triangle-up', 'x'],
			color_discrete_map=colors, size='dummy_size', size_max=15,
		)
		fig.update_layout(
			legend=dict(x=1.02, y=1, bordercolor='white', borderwidth=1, font=dict(size=20)),
			title='AT-GC skews ' + title, template='plotly_white',
			autosize=False, width=2000, height=1000,
		)
		fig.update_xaxes(title=f'PC1 - {per_var[0]}%', title_font=dict(size=15))
		fig.update_yaxes(title=f'PC2 - {per_var[1]}%', title_font=dict(size=15))
		return fig

	# ── 5. RSCU statistics ─────────────────────────────────────────────────────
	def perform_RSCU_stats(data, title, output_stats_folder):
		data		  = data.copy()
		data['TEST']  = data['Charset'] + '_' + data['Codon']
		group_domains = data['TEST'].unique()
		results = []
		for i in range(len(group_domains)):
			for j in range(i + 1, len(group_domains)):
				if group_domains[i].split('_')[-1] == group_domains[j].split('_')[-1]:
					data_1 = data[data['TEST'] == group_domains[i]]['RSCU']
					data_2 = data[data['TEST'] == group_domains[j]]['RSCU']
					stat, p_val = spstats.mannwhitneyu(data_1, data_2)
					results.append({
						'Group_Domain_1':		 group_domains[i],
						'Group_Domain_2':		 group_domains[j],
						'MannWhitneyU_statistic': stat,
						'p_value':				p_val,
					})
		results_df					= pd.DataFrame(results)
		results_df['corrected_p_value'] = multipletests(results_df['p_value'], method='bonferroni')[1]
		mean_stdev		 = data.groupby('TEST').agg({'RSCU': ['mean', 'std']})
		mean_stdev.columns = ['mean_Freq', 'stdev_Freq']
		mean_stdev		 = mean_stdev.reset_index()
		modified_title = title.replace(' ', '')
		results_df.to_csv(
			os.path.join(output_stats_folder, f'RSCU_{modified_title}_wilcoxon-test.tsv'),
			sep='\t', index=False,
		)
		mean_stdev.to_csv(
			os.path.join(output_stats_folder, f'RSCU_{modified_title}_mean-stdev-test.tsv'),
			sep='\t', index=False,
		)

	# ── 6. Full pipeline ───────────────────────────────────────────────────────
	def run_full_analysis(folder, output_tables_folder, output_graphics_folder, output_stats_folder, code, gene_order):
		df_AA, df_RSCU, df_1_SKEW, df_23_SKEW = analyze_codon_usage(
			folder=folder,
			output_RSCU=os.path.join(output_tables_folder, 'RSCU.tsv'),
			output_AA=os.path.join(output_tables_folder, 'AA_freq.tsv'),
			output_1_SKEW=os.path.join(output_tables_folder, 'First_skew.tsv'),
			output_23_SKEW=os.path.join(output_tables_folder, 'Second_Third_skew.tsv'),
			code=code,
		)
		fig_domains = perform_hist(df_AA[df_AA['Charset'].isin(DOMAIN_CHARSETS)], 'by domain', output_stats_folder)
		fig_domains.write_html(os.path.join(output_graphics_folder, 'AminoAcid_frequency_by_domain.html'))
		_safe_write_image(fig_domains, os.path.join(output_graphics_folder, 'AminoAcid_frequency_by_domain.pdf'))
		if gene_order:
			fig_chains = perform_hist(df_AA[df_AA['Charset'].isin(CHAIN_CHARSETS)], 'by chain', output_stats_folder)
			fig_chains.write_html(os.path.join(output_graphics_folder, 'AminoAcid_frequency_by_chain.html'))
			_safe_write_image(fig_chains, os.path.join(output_graphics_folder, 'AminoAcid_frequency_by_chain.pdf'))
		df_1_domains = (
			df_1_SKEW[df_1_SKEW['Charset'].isin(DOMAIN_CHARSETS)]
			.pipe(lambda d: d[d['CodonPosition'].isin(CODON_POSITIONS)])
			.reset_index(drop=True)
		)
		fig_gc_domains = perform_boxplot(df_1_domains, 'by domain', DOMAIN_COLORS, output_stats_folder)
		fig_gc_domains.write_html(os.path.join(output_graphics_folder, 'GC_frequency_by_domain.html'))
		_safe_write_image(fig_gc_domains, os.path.join(output_graphics_folder, 'GC_frequency_by_domain.pdf'))
		if gene_order:
			df_1_chain = (
				df_1_SKEW[df_1_SKEW['Charset'].isin(CHAIN_CHARSETS)]
				.pipe(lambda d: d[d['CodonPosition'].isin(CODON_POSITIONS)])
				.reset_index(drop=True)
			)
			fig_gc_chains = perform_boxplot(df_1_chain, 'by chain', CHAIN_COLORS, output_stats_folder)
			fig_gc_chains.write_html(os.path.join(output_graphics_folder, 'GC_frequency_by_chain.html'))
			_safe_write_image(fig_gc_chains, os.path.join(output_graphics_folder, 'GC_frequency_by_chain.pdf'))
		df_skews   = df_23_SKEW[df_23_SKEW['Position'] != 'All']
		df_domains = df_skews[df_skews['Charset'].isin(DOMAIN_CHARSETS)].reset_index(drop=True)
		fig_pca_dom = pca_plot(df_domains, 'by Domains', DOMAIN_COLORS)
		fig_pca_dom.write_html(os.path.join(output_graphics_folder, 'Skew-PCA_by_domain.html'))
		_safe_write_image(fig_pca_dom, os.path.join(output_graphics_folder, 'Skew-PCA_by_domain.pdf'))
		fig_sct_dom = scatter_plot(df_domains, 'by Domains', DOMAIN_COLORS)
		fig_sct_dom.write_html(os.path.join(output_graphics_folder, 'Skew-scatter-plot_by_domain.html'))
		_safe_write_image(fig_sct_dom, os.path.join(output_graphics_folder, 'Skew-scatter-plot_domain.pdf'))
		if gene_order:
			df_chain = df_skews[df_skews['Charset'].isin(CHAIN_CHARSETS)].reset_index(drop=True)
			fig_pca_chn = pca_plot(df_chain, 'by Chains', CHAIN_COLORS)
			fig_pca_chn.write_html(os.path.join(output_graphics_folder, 'Skew-PCA_by_chain.html'))
			_safe_write_image(fig_pca_chn, os.path.join(output_graphics_folder, 'Skew-PCA_by_chain.pdf'))
			fig_sct_chn = scatter_plot(df_chain, 'by Chains', CHAIN_COLORS)
			fig_sct_chn.write_html(os.path.join(output_graphics_folder, 'Skew-scatter-plot_by_chain.html'))
			_safe_write_image(fig_sct_chn, os.path.join(output_graphics_folder, 'Skew-scatter-plot_chain.pdf'))
		perform_RSCU_stats(df_RSCU[df_RSCU['Charset'].isin(DOMAIN_CHARSETS)], 'by domain', output_stats_folder)
		if gene_order:
			perform_RSCU_stats(
				df_RSCU[df_RSCU['Charset'].isin(CHAIN_CHARSETS + ['Combined'])],
				'by chain', output_stats_folder,
			)

	# ── eztrampo_main ──────────────────────────────────────────────────────────
	def eztrampo_main(files, gene_order, genetic_code, model, thmm_tables, outdir, tmp, plots, tables, stats, partitions_dir):
		files = [os.path.join(files, f) for f in os.listdir(files)]
		num_of_files(files)
		new_names = revise_names(files, tmp)
		for fasta in new_names:
			processed_file = check_fasta(fasta, tmp)
			os.remove(fasta)
			check_length(processed_file)
			degapped_file = remove_gaps(processed_file)
			os.remove(processed_file)
			isIUPAC(degapped_file)
			clean_file = isStopCodon(degapped_file, genetic_code)
			os.remove(degapped_file)
			translated_file = translate(clean_file, genetic_code, tmp)
			isIUPAC_prot(translated_file)
			gene = os.path.splitext(os.path.basename(fasta))[0]
			alignment_aa_file = alignment(translated_file, tmp)
			with open(os.path.join(tmp, 'tmp.fas'), 'w') as outfile:
				for record in SeqIO.parse(model, 'fasta'):
					if gene in record.id:
						SeqIO.write(record, outfile, 'fasta')
			alignment_nt_ref_file, alignment_aa_ref_file = alignment_ref(alignment_aa_file, clean_file, os.path.join(tmp, 'tmp.fas'), tmp)
			fasta2nexus(alignment_nt_ref_file, 'DNA')
			fasta2nexus(alignment_aa_ref_file, 'protein')
			partitions(alignment_aa_ref_file, thmm_tables, os.path.splitext(os.path.basename(model))[0], tmp)
			for f in [clean_file, translated_file, os.path.join(tmp, 'tmp.fas'),
					  alignment_aa_file, alignment_nt_ref_file, alignment_aa_ref_file]:
				os.remove(f)
		old_chars = extract_old_charsets(tmp)
		lengths   = extract_lengths(tmp)
		p39  = build_new_intervals(lengths, old_chars, partitions_dir)
		p117 = add_codonsets2charsets(p39)
		p3   = process_domain_rc1(p117)
		merge_new_charsets_regions_codons_rc1(p117)
		merge_new_charsets_genes_rc1(p39)
		checkfiles = ['117p_gen_dom_cod.nex', '13p_gen.nex', '39p_gen_dom.nex', '3p_dom.nex', '9p_dom_cod.nex']
		if gene_order:
			p18 = merge_new_charsets_codon_strand(gene_order, p117, 117)
			process_codon_rc1(p18)
			process_codon_strand_rc1(p18)
			process_12part_MAIM_merged_rc1(p18)
			p6  = merge_new_charsets_codon_strand(gene_order, p39, 39)
			checkfiles = ['117p_gen_dom_cod.nex', '12p_MAIM_merged.nex', '13p_gen.nex',
					  '18p_str_dom_cod.nex', '39p_gen_dom.nex', '3p_cod.nex',
					  '3p_dom.nex', '6p_cod_str.nex', '6p_str_dom.nex', '9p_dom_cod.nex']
			
		missing = [os.path.join(partitions_dir, f) for f in checkfiles
				   if not os.path.exists(os.path.join(partitions_dir, f))]
		if missing:
			raise ValueError(f"Missing partition files: {missing}")
		NT = os.path.join(tmp, "NT")
		os.makedirs(NT, exist_ok=True)
		for fname in os.listdir(tmp):
			if fname.endswith("ali_nt_ref.nexus"):
				src	 = os.path.join(tmp, fname)
				tmpfile = src + ".tmp"
				with open(src) as fin, open(tmpfile, "w") as fout:
					for line in fin:
						if line.strip().lower().startswith("begin sets;"):
							break
						fout.write(line)
				shutil.move(tmpfile, os.path.join(NT, fname))
		concat_file  = concatenate(NT, 'infile')
		concat_file2 = concat_file.replace('infile', 'infile2')
		shutil.copyfile(concat_file, concat_file2)
		main_nexus = os.path.join(outdir, 'nt_combined.nex')
		with open(concat_file) as src, open(main_nexus, 'w') as dst:
			for line in src:
				if line.strip().lower().startswith("begin sets;"):
					break
				dst.write(line)
		with open(p39) as src, open(concat_file, "a") as dst:
			dst.write(src.read())
		with open(p3) as src, open(concat_file2, "a") as dst:
			dst.write(src.read())
		split_partitions(concat_file,  os.path.join(tmp, 'nt_combined_genes_regions'))
		split_partitions(concat_file2, os.path.join(tmp, 'nt_combined_regions'))
		if gene_order :
			concat_file3 = concat_file.replace('infile', 'infile3')
			shutil.copyfile(concat_file, concat_file3)
			with open(p6) as src, open(concat_file3, "a") as dst:
				dst.write(src.read())
			split_partitions(concat_file3, os.path.join(tmp, 'nt_combined_strand_regions'))
		run_full_analysis(
			folder				 = tmp,
			output_tables_folder   = tables,
			output_graphics_folder = plots,
			output_stats_folder	= stats,
			code				   = genetic_code,
			gene_order			 = gene_order,
		)

	# ── output folders ─────────────────────────────────────────────────────────
	tmp			= os.path.join(outdir, 'tmp')
	plots		  = os.path.join(outdir, 'plots')
	tables		 = os.path.join(outdir, 'tables')
	stats		  = os.path.join(outdir, 'stats')
	partitions_dir = os.path.join(outdir, 'partitions')

	for folder in [tmp, plots, tables, stats, partitions_dir]:
		os.makedirs(folder, exist_ok=True)

	# ── validate optional user files ───────────────────────────────────────────
	if sequence and not thmm_tables:
		raise ValueError('If you pass a custom sequence file you must also provide the TMHMM tables (-t)')
	if thmm_tables and not sequence:
		raise ValueError('If you pass custom TMHMM tables you must also provide the sequence file (-s)')

	if thmm_tables and sequence:
		if model.lower() != 'user':
			raise ValueError('If you pass custom table and sequence files, select "user" as the model organism (-m user)')
		user = os.path.join(outdir, 'user')
		os.makedirs(user, exist_ok=True)
		processed_file = check_fasta(sequence, user)
		degapped_file  = remove_gaps(processed_file)
		os.remove(processed_file)
		isIUPAC_prot(degapped_file)
		clean_file = isTerminalStop(degapped_file)
		os.remove(degapped_file)
		new_name = revise_user_sequences(clean_file, user)
		os.remove(clean_file)
		model	   = new_name
		thmm_tables = check_tables(thmm_tables, user)
	else:
		model	   = os.path.join('templates/sequences/',			 args.model + '.fas')
		thmm_tables = os.path.join('templates/model_organism_tables/', args.model + '.tsv')

	if model == 'cel':
		print("Warning: the selected model organism (C. elegans) lacks the ATP8 gene — it will be skipped")
	if gene_order:
		gene_order = os.path.join('templates/go', gene_order)

	eztrampo_main(fasta_files, gene_order, genetic_code, model, thmm_tables,
				  outdir, tmp, plots, tables, stats, partitions_dir)
				  
	if os.path.exists(tmp):
		shutil.rmtree(tmp)


def main():
	parser = argparse.ArgumentParser(
		description="Choose the EZmito program you want to use.",
		formatter_class=argparse.RawTextHelpFormatter
	)
	subparsers = parser.add_subparsers(
		title="Available commands",
		help="Use one of these subcommands for the corresponding analysis.",
		dest="command"
	)

	# EZcircular subcommand
	parser_ez_circular = subparsers.add_parser(
		'ezcircular',
		help='Circularize a sequence from another starting gene',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_circular.add_argument("-f", "--feature", help="Genome feature. [circular or linear]", default='circular', type=str)
	parser_ez_circular.add_argument("-s", "--start", help="Starting gene", default='cox1', type=str)
	parser_ez_circular.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	required_ez_circular = parser_ez_circular.add_argument_group('required named arguments')
	required_ez_circular.add_argument("-i", "--input", required=True, help="FASTA input file")
	required_ez_circular.add_argument("-b", "--bed", required=True, help="BED input file")
	parser_ez_circular.set_defaults(func=ez_circular_subcommand)

	# EZcodon subcommand
	parser_ez_codon = subparsers.add_parser(
		'ezcodon',
		help='Analyze codon usage of J (heavy) and N (light) chains',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_codon.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	required_ez_codon = parser_ez_codon.add_argument_group('required named arguments')
	required_ez_codon.add_argument("-J", "--heavy", help="Path to J chain FASTA files")
	required_ez_codon.add_argument("-N", "--light", help="Path to N chain FASTA files")
	required_ez_codon.add_argument("-c", "--code", required=True, help=code_help, type=int)
	parser_ez_codon.set_defaults(func=ez_codon_subcommand)

	# EZmap subcommand
	parser_ez_map = subparsers.add_parser(
		'ezmap',
		help='Create a custom map plot of your mitogenome',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_map.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	parser_ez_map.add_argument("-f", "--feature", help="Genome feature. [circular or linear]", default='circular', type=str)
	parser_ez_map.add_argument("-colJ", "--colorJ", help="J (heavy) strand color", default='#add8e6', type=str)
	parser_ez_map.add_argument("-colN", "--colorN", help="N (light) strand color", default='#B22222', type=str)
	required_ez_map = parser_ez_map.add_argument_group('required named arguments')
	required_ez_map.add_argument("-g", "--gff", required=True, help="GFF3 input file")
	parser_ez_map.set_defaults(func=ez_map_subcommand)

	# EZmix subcommand
	parser_ez_mix = subparsers.add_parser(
		'ezmix',
		help='Check for possible chimeras occurred during your mitogenome assemblies',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_mix.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	parser_ez_mix.add_argument("-bn", "--blastn", help="BLASTn path", default='', type=str)
	parser_ez_mix.add_argument("-id", "--identity", help="Identity threshold [0.5-1]", default=0.95, type=float)
	parser_ez_mix.add_argument("-len", "--length", help="Length threshold (bp)", default=200, type=int)
	required_ez_mix = parser_ez_mix.add_argument_group('required named arguments')
	required_ez_mix.add_argument("-i", "--input", required=True, help="MultiFASTA input file")
	parser_ez_mix.set_defaults(func=ez_mix_subcommand)

	# EZpipe subcommand
	parser_ez_pipe = subparsers.add_parser(
		'ezpipe',
		help='Prepare J and N (heavy and light) FASTA files for phylogenetic analysis',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_pipe.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	parser_ez_pipe.add_argument("-p", "--positions", help="Number of codon positions to analyze (should the program remove 3rd codon positions?) [2,3]", default=3, type=int)
	required_ez_pipe = parser_ez_pipe.add_argument_group('required named arguments')
	required_ez_pipe.add_argument("-i", "--input", required=True, help="Path to FASTA files")
	required_ez_pipe.add_argument("-c", "--code", required=True, help=code_help, type=int)
	parser_ez_pipe.set_defaults(func=ez_pipe_subcommand)

	# EZskew subcommand
	parser_ez_skew = subparsers.add_parser(
		'ezskew',
		help='Analyze codon skew biases of J (heavy) and N (light) chains',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_skew.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	required_ez_skew = parser_ez_skew.add_argument_group('required named arguments')
	required_ez_skew.add_argument("-J", "--heavy", help="Path to J chain FASTA files")
	required_ez_skew.add_argument("-N", "--light", help="Path to N chain FASTA files")
	required_ez_skew.add_argument("-c", "--code", required=True, help=code_help, type=int)
	parser_ez_skew.set_defaults(func=ez_skew_subcommand)

	# EZsplit subcommand
	parser_ez_split = subparsers.add_parser(
		'ezsplit',
		help='Split multiple Genbank downloaded multifasta to a set of PCGs',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_split.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	required_ez_split = parser_ez_split.add_argument_group('required named arguments')
	required_ez_split.add_argument("-g", "--gff", required=True, help="GFF3 input file")
	required_ez_split.add_argument("-i", "--input", required=True, help="MultiFASTA input file of complete mitogenomes")
	parser_ez_split.set_defaults(func=ez_split_subcommand)

	# EZtrampo subcommand
	parser_ez_trampo = subparsers.add_parser(
		'eztrampo',
		help='Partition mitochondrial genomes depending on membrane domains',
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser_ez_trampo.add_argument("-o", "--outdir", help="Output directory", default='outdir', type=str)
	parser_ez_trampo.add_argument("-s", "--sequence", help="User's model organism genes file in FASTA (amino acid) format", type=str)
	parser_ez_trampo.add_argument("-t", "--tables", help="Path to user's model organism table files in TMHMM format", type=str)
	parser_ez_trampo.add_argument("-g", "--gene_order", help="Gene order model nickname to be employed during the analysis\n"
								  "\tAll vertebrates: vert\n"
								  "\tArthropods: panc\n"
								  "\tLumbricus terrestris, Caenorhabditis elegans, Metridium senile: ances\n"
								  "\tAlbinaria caerulea: albin\n"
								  "\tMetacangronyx: meta", type=str)
	parser_ez_trampo.add_argument("-n", "--threads", help="Number of threads to employ in MAFFT", default=1, type=int)
	required_ez_trampo = parser_ez_trampo.add_argument_group('required named arguments')
	required_ez_trampo.add_argument("-p", "--path", required=True, help="Path to FASTA files directory")
	required_ez_trampo.add_argument("-c", "--code", required=True, help=code_help, type=int)
	required_ez_trampo.add_argument("-m", "--model", required=True, help="Nickname of model organism\n"
								  "\tUser table (if you pass -s and -t): user\n"
								  "\tHomo sapiens (Chordata): hsa\n"
								  "\tPatiria pectinifera (Echinodermata): ppe\n"
								  "\tDrosophila melanogaster (Pancrustacea+Chelicerata): dme\n"
								  "\tAlbinaria caerulea (Mollusca): aca\n"
								  "\tLumbricus terrestris (Annelida): lte\n"
								  "\tCaenorhabditis elegans (Nematoda): cel\n"
								  "\tMetridium senile (Cnidaria): mse", type=str)
	parser_ez_trampo.set_defaults(func=ez_trampo_subcommand)

	# ── Parse and dispatch ─────────────────────────────────────────────────────
	args = parser.parse_args()
	if args.command is None:
		parser.print_help()
	else:
		args.func(args)


if __name__ == "__main__":
	main()
	end_time = time.time()
	runtime = round(end_time - start_time, 2)
	print(f"\n\n-------------------The process was correctly completed in {runtime} seconds-------------------")
	print("Thank you for using this code. If it helped you, please cite:")
	print("Cucini C., Leo C., Iannotti N., Boschi S., Brunetti C., Pons J., Fanciulli P. P., Frati F., Carapelli A., & Nardi F. (2021)")
	print("EZmito: a simple and fast tool for multiple mitogenome analyses, Mitochondrial DNA Part B, 6(3), 1101-1109.")
	print("Doi: 10.1080/23802359.2021.1899865")

