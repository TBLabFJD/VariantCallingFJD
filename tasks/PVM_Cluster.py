import sys
import os
from glob import glob
import subprocess
import argparse
from collections import Counter
from collections import defaultdict
import csv
import datetime
import time
import re


start_time = time.time()


parser = argparse.ArgumentParser(description="Post VEP Annotation FJD")
parser.add_argument('input_file', help='\t\t Text file after vep2tsv')
parser.add_argument('-k', '--plinkfile', help='\t\t LOH file', required=False)
parser.add_argument('-o', '--output', dest="output",help='\t\tOutput file name.', required=True)
parser.add_argument('-P', '--pathology', help="Disease group name: resp, digest, skin, musc, gu, pregnancy, perinatal, cong_mal, clinic, infectious, other, neoplasm, blood, endoc, mental, CNS, eye, ear, circ,healthy.", default='healthy')
parser.add_argument('-l', '--local', help="Run in local.", required=True)


args = parser.parse_args()

if args.local == "True":
	Spanish_freq_disease="/mnt/genetica2/NGS_data/CSVS_frequencies/dis_group.csv"
	Spanish_freq_folder="/mnt/genetica2/NGS_data/CSVS_frequencies"
	Region_dict="/mnt/genetica3/raquel/VEP/dict_region.csv"
	Gene_dict="/mnt/genetica3/marius/pipeline_practicas_marius/software/variant_effect_predictor/.vep/dbs/dbNSFP3.5_gene"

else:	
	Spanish_freq_disease="/home/proyectos/bioinfo/references/CSVS_frequencies/dis_group.csv"
	Spanish_freq_folder="/home/proyectos/bioinfo/references/CSVS_frequencies"
	Region_dict="/home/proyectos/bioinfo/references/other/dict_region.csv"
	Gene_dict="/home/proyectos/bioinfo/references/VEPdbs/dbNSFP3.5_gene"
# RAQUEL

# Upload eq table of each freq table with pathology group:

print("Loading Spanish Frequency Dictionary")

with open(Spanish_freq_disease) as csvfile:
    Pathology=[]
    data=[]
    reader = csv.reader(csvfile, delimiter='\t')
    for x in reader:
        data.append( x[0])
        Pathology.append(x[1])

    path_dict = dict(zip(Pathology, data))

if args.pathology in Pathology:
	 file = path_dict[args.pathology]
	 dir=Spanish_freq_folder
	 full_path = os.path.join(dir, file)
else:
	 sys.stderr.write("ERROR IN PATHOLOGY OPTION")


with open( full_path, "r") as f:
# with open( '/mnt/genetica2/NGS_data/CSVS_frequencies/all_less_group1.csv', "r") as f:

	Keylist=[]
	Values=[]
	f.readline
	for x in f:
		myxList=x.split('\t')
		Keylist.append('chr' + myxList[0] + '|' + str(myxList[1]) + '|' + myxList[2] + '|' + myxList[3])
		Values.append( myxList[10])

	freq_dict = dict(zip(Keylist, Values))

print("Done with Spanish Frequency Dictionary For Pathology: " + args.pathology)

# Dictionary of genomic regions (RAQUEL)


with open(Region_dict) as dict_file:
	key=[]
	region=[]
	reader = csv.reader(dict_file, delimiter=',')
	for line in reader:
		key.append(line[0])
		region.append(line[1])

	region_dict=dict(zip(key, region))

print("Done with Genomic Regions Dictionary")

# Dictionary of Gene info

with open(Gene_dict) as gene_file:
    Key=[]
    Values=[]
    gene_file.readline()
    for line in gene_file:
        line = line.split('\t')
        Genekey=line[0]
        Key.append(Genekey)
        fullname=line[12]
        nim_dis=line[22]
        Value = (fullname,  nim_dis)
        Values.append(Value)
    
    gene_dict= dict(zip(Key, Values))
print("Done with Genomic Info Dictionary")

# RAQUEL

# Getting Spanish frequencies funtion:

def add_spanish_freq( chr, coord, ref, alt, freq_dict):
	
	Variations = alt.split(',')
	spanish_freqs=[]
	for i in range(len(Variations)) :

		# Create two keys for each INDEL

		Multi_key = None # Default value
		Multi_key_indel = None # Default value

		# Insertions:
		if len(Variations[i]) > 1 and len(str(ref)) == 1 :
			Multi_key_indel = str(chr) + '|' + str(int(coord) + 1) + '|.|' + Variations[i][1:]
			Multi_key = str(chr) + '|' + str(coord) + '|' + str(ref) + '|' + Variations[i]

		# Delections:
		elif len(Variations[i]) == 1 and len(str(ref)) > 1 :
			Multi_key_indel = str(chr) + '|' + str(int(coord) + 1) + '|' + str(ref)[1:] + '|.'
			Multi_key = str(chr) + '|' + str(coord) + '|' + str(ref) + '|' + Variations[i]

		# SNP:
		else:
			Multi_key = str(chr) + '|' + str(coord) + '|' + str(ref) + '|' + Variations[i]


		# Check if the keys are in the dictionary:
		
		# spanish_freq = None
		if Multi_key in freq_dict:
			spanish_freqs.append(freq_dict[Multi_key])


		elif Multi_key_indel in freq_dict and Multi_key not in freq_dict:
			spanish_freqs.append(freq_dict[Multi_key_indel])

		else:
			spanish_freqs.append('NaN')

	spanish_freqs = ",".join(spanish_freqs)

	return spanish_freqs



# GETTING GENE INFO

def get_gene_info(gene):
    
    Full_name = None # Default value
    Disease = None # Default value
    
    if gene in gene_dict:
        Full_name = gene_dict[gene][0]
        Disease = gene_dict[gene][1]
    else:
        Full_name = "Not found"
        Disease = "Not found"
        
    return(Full_name, Disease)

def easy_nomenclature(chr, coord, alt):
	
	# Find chr:star-end variable

	alt_split= alt.split(',')
	end_list=[]
	for i in range(len(alt_split)):
		end = int(coord) + (len(alt_split[i]) - 1)
		end_list.append(str(end))

	ends = ".".join(end_list)

	chr_start_end = chr[3:] + ':' + str(coord) + "-" + str(ends)
	return(chr_start_end)

def genomicRegion(vep_consequence):

	priority_list = ["SPLICING", "5UTR", "3UTR", "ncRNA", "regulatory", "UPSTREAM", "DOWNSTREAM", "EXONIC", "INTRONIC", "INTERGENIC", "NA"]
	regions = []
	
	vep_consequence = vep_consequence.split("&")

	for term in vep_consequence: 
		region = region_dict[term]
		regions.append(region)

	value = "&".join(regions)
	my_region = [re.findall('(?i)(' + v + ').*?', value) for v in priority_list if re.findall('(?i)(' + v + ').*?', value)][0][0]

	return my_region



#LORENA


def checkOverlapping(chr, coord, lohs):

	flag=False
	for lohs in lohs:
		loh = lohs.split("_")
		if ("chr"+loh[0]) == chr and int(loh[1]) <= coord and int(loh[2])>=coord:
			flag=True

	return flag


# LORENA

chrDicc = defaultdict(list)

if args.plinkfile!=None:

	with open(args.plinkfile, "r") as plink:
		header=plink.readline()
		for line in plink:
			linesplitok = [x for x in line.split() if x!=""]
			sample = linesplitok[1]
			chr = linesplitok[3]
			start = linesplitok[6]
			end = linesplitok[7]

			chrDicc[sample].append("_".join([chr,start,end]))


output = open(args.output,"w")  
with open(args.input_file, "r") as snp:
	header=snp.readline().strip()
	fields = header.split()
	fields = [s.replace('VEP_', '') for s in fields]
	
	samples = [x.split("_GT")[0] for x in fields if x[-3:] == "_GT"]
	
	if args.plinkfile!=None:

		header= "\t".join(fields[0:4]+["Chr:Start-end"]+[fields[4]]+["gene_full_name"]+fields[6:11]+["Genomic_region"]+fields[11:13]+["CANONICAL"]+fields[14:26]+["NIM_Disease"]+fields[26:70]+["spanish_freq"]+fields[70:-2]+[s + "_VF" for s in samples] + fields[-2:] + [s + "_LOH" for s in samples])
	
	else:
		header= "\t".join(fields[0:4]+["Chr:Start-end"]+[fields[4]]+["gene_full_name"]+fields[6:11]+["Genomic_region"]+fields[11:13]+["CANONICAL"]+fields[14:26]+["NIM_Disease"]+fields[26:70]+["spanish_freq"]+fields[70:-2]+[s + "_VF" for s in samples] + fields[-2:])


	output.write(header + "\n")
	
	for line in snp:
		line = line.strip()
		linesplit = line.split('\t')

		chr = linesplit[0]
		coord = linesplit[1]
		ref = linesplit[2]
		alt = linesplit[3]
		gene = linesplit[4]
		vep_consequence=linesplit[15]


		canonical = linesplit[fields.index('CANONICAL')]
		
		# Fill empty spaces in VEP_CANONICAL
		
		my_canonical=[]
		if canonical in (None, ""):
			my_canonical.append("NO")
		else:
			my_canonical.append("YES")


		# Find corresponding regions
		
#		 dict_key_list=[]
		regions =[]
		region = genomicRegion(vep_consequence)
		regions.append(region)


		my_chr_start_end =[]
		chr_start_end = easy_nomenclature(chr, coord, alt)
		my_chr_start_end.append(chr_start_end)
		
		# Find Variant Frequency
		
		my_vf_list=[]
		myloh=[]
		for sample in samples:
			# vf_dict['vf_lst_%s' % sample] = []
			
			myvf=[]

			ad_field = sample + '_AD'
			dp_field = sample + '_DP'
			vf_field = sample + '_VF'
#		 
			ad_split = linesplit[fields.index(ad_field)].split(',')
		
			for i in range(1, len(ad_split)):
				dp_value = linesplit[fields.index(dp_field)]
				
				if dp_value == "0" or dp_value == "NA" or dp_value == ".":
					myvf.append("")

				elif(ad_split[i] == "NA" or ad_split[i] == "."):
					myvf.append("")

				else:
					vfs = (float(ad_split[i]) / float(dp_value)) *100
					myvf.append(str(vfs))
			
			myvf = ",".join(myvf)
			my_vf_list.append(myvf)

		# Find homozygosity regions
			
			if args.plinkfile!=None:

				if sample in chrDicc:
					x=checkOverlapping(chr, int(coord), chrDicc[sample])
					myloh.append(str(x))
				else:
					myloh.append("False")

		
		
			


		# Get Spanish Frequencies:

		spanish_freq_list =[]
		freq = add_spanish_freq(chr, coord, ref, alt, freq_dict)
		spanish_freq_list.append(freq)

		# Get gene info

		my_full_name=[]
		my_Disease =[]
		Full_name, Disease = get_gene_info(gene)
		my_full_name.append(Full_name)
		my_Disease.append(Disease)



		fline= "\t".join(linesplit[0:4]+my_chr_start_end+[linesplit[4]]+my_full_name+linesplit[6:11]+regions+linesplit[11:13]+my_canonical+linesplit[14:26]+my_Disease+linesplit[26:70]+spanish_freq_list+linesplit[70:-2]+my_vf_list+linesplit[-2:]+myloh)
		output.write(fline + "\n")


