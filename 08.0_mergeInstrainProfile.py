#!/usr/bin/env python3

"""Merge inStrain-Profile genome outputs


Usage: python3 s18.0_mergeInStrainProfile.py <Sample list> <inStrain-Profile dir>

Args:
	Sample list: File listing all samples (i.e., mapping file)
  inStrain-Profile dir: path to the directory with all inStrain-Profile outputs (ends in "/")
"""

import sys

#Print out the doc string and exit if the number of input parameters is not correct
if len(sys.argv) != 3:
  print(__doc__)
  exit()

sample_list = sys.argv[1]
outdir = sys.argv[2]

samples = []

for sample in open(sample_list, "r"):
  samples.append(sample.strip())

merged_path = outdir + "merged.tsv"
merged_file = open(merged_path, "w")

print("sample\tgenome\tcoverage\tbreadth\tnucl_diversity\tlength\ttrue_scaffolds\tdetected_scaffolds\tcoverage_median\tcoverage_std\tcoverage_SEM\tbreadth_minCov\tbreadth_expected\tnucl_diversity_rarefied\tconANI_reference\tpopANI_reference\tiRep\tiRep_GC_corrected\tlinked_SNV_count\tSNV_distance_mean\tr2_mean\td_prime_mean\tconsensus_divergent_sites\tpopulation_divergent_sites\tSNS_count\tSNV_count\tfiltered_read_pair_count\treads_unfiltered_pairs\treads_mean_PID\treads_unfiltered_reads\tdivergent_site_count", file = merged_file)

for sample in samples:
  output_path = outdir + sample + ".IS/output/" + sample + ".IS_genome_info.tsv"
  for line in open(output_path, "r"):
    if line.startswith("genome"):
      continue
    else:
      line = line.strip()
      line = sample + "\t" + line
      print(line, file = merged_file)

merged_file.close()
