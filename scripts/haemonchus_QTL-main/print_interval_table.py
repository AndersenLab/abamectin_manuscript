#!/usr/bin/env python2.7
import sys

CELE_chr = 'V'
HCON_chr = 'hcontortus_chr5_Celeg_TT_arrow_pilon'

def parse_CELE_interval_file(CELE_interval_file):
	with open(CELE_interval_file, 'r') as CELE_intervals:
		CELE_interval_dict = {} # to store coordinates
		CELE_interval_list = [] # to store order, which is needed to make a rational output file 
		for line in CELE_intervals:
			columns = line.rstrip("\n").split("\t")
			id, start, stop = columns[0], int(columns[1]), int(columns[2])
			CELE_interval_list.append(id)
			CELE_interval_dict[id] = [start, stop]
	return CELE_interval_dict, CELE_interval_list

def parse_HCON_interval_file(HCON_interval_file):
	with open(HCON_interval_file, 'r') as HCON_intervals:
		HCON_interval_dict = {} # to store coordinates
		for line in HCON_intervals:
			columns = line.rstrip("\n").split("\t")
			id, start, stop = columns[0], int(columns[1]), int(columns[2])
			HCON_interval_dict[id] = [start, stop]
	return HCON_interval_dict

def parse_CELE_gff3_file(CELE_gff3_file):
	with open(CELE_gff3_file, 'r') as CELEG_gff3:
		CELEG_gff3_dict = {}
		for line in CELEG_gff3:
			if not line.startswith("#"):
				columns = line.rstrip("\n").split("\t")
				chr, feature, start, stop, col9 = columns[0], columns[2], int(columns[3]), int(columns[4]), columns[8]
				if feature == "mRNA" and chr == CELE_chr:
					transcriptID = ".".join(col9.split(";")[0].split(":")[1].split(".")[0:2])
					geneID = col9.split(";")[1].split(":")[1]
					locusID = '-'
					for field in col9.split(";"):
						if field.startswith("locus="):
							locusID = field.split("=")[1]
					CELEG_gff3_dict[transcriptID] = [chr, start, stop, geneID, locusID]
		return CELEG_gff3_dict

def parse_HCON_gff3_file(HCON_gff3_file):
	with open(HCON_gff3_file, 'r') as HCON_gff3:
		HCON_gff3_dict = {}
		for line in HCON_gff3:
			if not line.startswith("#"):
					columns = line.rstrip("\n").split("\t")
					chr, feature, start, stop, col9 = columns[0], columns[2], int(columns[3]), int(columns[4]), columns[8]
					if feature == "mRNA" and chr == HCON_chr:
						transcriptID = col9.split(";")[0].split(":")[1]
						geneID = col9.split(";")[1].split(":")[1]
						HCON_gff3_dict[transcriptID] = [chr, start, stop, geneID]
		return HCON_gff3_dict

def parse_single_copy_ortho_file(single_copy_ortho_file):
	with open(single_copy_ortho_file, 'r') as single_copy_orthos:
		single_copy_ortho_dict = {}
		for line in single_copy_orthos:
			for seq in line.rstrip("\n").split(" ")[1:]:
				if seq.startswith("HCON_"):
					HCON_seq = seq
				else: 
					CELE_seq = seq
			single_copy_ortho_dict[CELE_seq] = HCON_seq
	return single_copy_ortho_dict

def print_table(CELE_interval_dict, CELE_interval_list, CELEG_gff3_dict, HCON_gff3_dict, single_copy_ortho_dict):
	print "\t".join(["transcriptID", "WormBaseID", "locus", "start", "stop", "VL_am", "VL_lm", "VL_nil", "VC_lm", "VC_nil", "VR_am", "VR_lm", "VR_nil", "HCON_ortho", "HCON_start", "HCON_stop", "interval_1_status", "interval_2_status"])
	for transcriptID, info_list in CELEG_gff3_dict.iteritems(): # for every transcript in the CELE gff3 file
		start, stop, geneID, locusID = info_list[1], info_list[2], info_list[3], info_list[4] # get the start and stop coordinates, geneID and locusID
		interval_presence_list = [] # make a list to store the output for the table
		for intervalID in CELE_interval_list: # for each CELE interval in the ORDERED interval list
			interval_start, interval_stop = CELE_interval_dict[intervalID][0], CELE_interval_dict[intervalID][1] # get the start and stop coordinates	
			if start >= interval_start and stop <= interval_stop: # if the transcript falls within the interval
				interval_presence_list.append("Y") # add Y to the output
			else: 
				interval_presence_list.append("-") # add - to the output
		if "Y" in interval_presence_list: # if the transcript is found within one of the intervals
			try:
				HCON_ortho = single_copy_ortho_dict[transcriptID]
				HCON_start, HCON_stop = HCON_gff3_dict[HCON_ortho][1], HCON_gff3_dict[HCON_ortho][2]
				interval_1_status, interval_2_status = '-', '-'
				if HCON_start >= HCON_interval_dict['interval_1'][0] and HCON_stop <= HCON_interval_dict['interval_1'][1]:
					interval_1_status = "Y"
				if HCON_start >= HCON_interval_dict['interval_2'][0] and HCON_stop <= HCON_interval_dict['interval_2'][1]:
					interval_2_status = "Y"
			except KeyError:
				HCON_ortho, HCON_start, HCON_stop, interval_1_status, interval_2_status = '-', '-', '-', '-', '-'
			print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (transcriptID, geneID, locusID, start, stop, "\t".join(interval_presence_list), HCON_ortho, HCON_start, HCON_stop, interval_1_status, interval_2_status)

if __name__ == "__main__":
	SCRIPT = "print_interval_table.py"
	try:
		CELE_gff3_file = sys.argv[1]
		CELE_interval_file = sys.argv[2]
		HCON_gff3_file = sys.argv[3]
		HCON_interval_file = sys.argv[4]
		single_copy_ortho_file = sys.argv[5]
	except IndexError:
		sys.exit("USAGE: python2.7 %s %s %s %s %s %s" % (SCRIPT, "[CELE_gff3_file]", "[CELE_interval_file]", "[HCON_gff3_file]", "[HCON_interval_file]", "[single_copy_ortho_file]"))
	CELE_interval_dict, CELE_interval_list = parse_CELE_interval_file(CELE_interval_file)
	HCON_interval_dict = parse_HCON_interval_file(HCON_interval_file)
	CELEG_gff3_dict = parse_CELE_gff3_file(CELE_gff3_file)
	HCON_gff3_dict = parse_HCON_gff3_file(HCON_gff3_file)
	single_copy_ortho_dict = parse_single_copy_ortho_file(single_copy_ortho_file)	
	print_table(CELE_interval_dict, CELE_interval_list, CELEG_gff3_dict, HCON_gff3_dict, single_copy_ortho_dict)				
