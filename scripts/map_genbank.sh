#!/bin/bash

set -beEu -o pipefail

metadata_csv=$1
fasta_dir=$2

# Run in a temporary directory, remove when done
tmpdir=$(mktemp --directory --tmpdir map_genbank.XXXXXXXX)
cd $tmpdir

# Use NCBI Virus search to download metadata for H5N1 collected in USA since 2023
# Influenza A virus NCBI Taxonomy ID = 11320
tax_id=11320
today=$(date +%F)
ncbi_virus_url='https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?fq=%7B%21tag%3DSeqType_s%7DSeqType_s%3A%28%22Nucleotide%22%29&fq=VirusLineageId_ss%3A%28'$tax_id'%29&fq=%7B%21tag=CollectionDate_dr%7DCollectionDate_dr:%5B2023-01-01T00:00:00.00Z%20TO%20'$today'T23:59:59.00Z%5D&fq=%7B%21edismax%20qf=Serotype_s%7DH5N1&fq=%7B%21tag=Country_s%7DCountry_s:%28%22USA%22%29&q=%2A%3A%2A&cmd=download&dlfmt=csv&fl=genbank_acc%3AAccVer_s%2Cisolate%3AIsolate_s%2Clocation%3ACountryFull_s%2Ccollected%3ACollectionDate_s%2Chost%3AHost_s%2Cbioproject_accession%3ABioProject_s%2Cstrain%3AStrain_s%2Cserotype%3ASerotype_s%2Csegment%3ASegment_s&sort=id+asc'

curl -fSs "$ncbi_virus_url" > ncbi_virus_metadata.csv

cat > ./csvToTab.py <<EOF
import sys,csv
for row in csv.reader(sys.stdin):
    row = [r.strip() for r in row]
    print("\t".join(row))
EOF

# Get a mapping of sample names (minus the "-original" that is not included in GenBank names)
# to SRA run accessions.
python3 ./csvToTab.py < $metadata_csv \
| cut -f 1,31 | tail -n+2 \
| awk -v 'FS=\t' -v 'OFS=\t' '{print $2, $1;}' | sort \
| sed -re 's/-original//; s/-repeat2?//;' \
| sed -re 's/([0-9]{2}-[0-9]{6}-[0-9]{3})-(300|MTM)/\1/g;' \
| sort \
    > sampleToSraRun.tsv

# Extract metadata for GenBank records with those sample names.
python3 ./csvToTab.py < ncbi_virus_metadata.csv > ncbi_virus_metadata.tsv
grep -Fwf <(cut -f 1 sampleToSraRun.tsv) ncbi_virus_metadata.tsv \
    > metadataForSraRuns.tsv

# Extract just the sample name (like 24-005334-001) from the full name and map to GenBank accession,
# full name and (numeric) segment.
cut -f 1,7,9 metadataForSraRuns.tsv \
| perl -wne 'chomp;
    ($gb, $fullName, $segment) = split(/\t/);
    if ($fullName =~ m@/(\d{2}-\d{6}-\d{3})@) {
      $sampleName = $1;
    } else { die "Cant get sample name from $fullName"; }
    print join("\t", $sampleName, $gb, $fullName, $segment) . "\n";' \
| sort \
    > sampleToGbAccNameSegment.tsv

# Join to get SRA run:
join -t$'\t' sampleToSraRun.tsv sampleToGbAccNameSegment.tsv \
    > sampleToSraGb.tsv

# Print output header
echo -e "seg_file\tseg_seq_name\tsra_run\tseg\tgenbank_acc\tgenbank_seg\tgenbank_name"

# Map SRA run to files in fasta_dir
# Some file names have M2, some files have MP; fudge both, then restrict to actual $fasta_dir files
cut -f 2- sampleToSraGb.tsv \
| perl -wne 'chomp;
    ($run, $gb, $fullName, $segNum) = split(/\t/);
    if ($segNum == 1) {
      @segNames = ("PB2");
    } elsif ($segNum == 2) {
      @segNames = ("PB1");
    } elsif ($segNum == 3) {
      @segNames = ("PA");
    } elsif ($segNum == 4) {
      @segNames = ("HA");
    } elsif ($segNum == 5) {
      @segNames = ("NP");
    } elsif ($segNum == 6) {
      @segNames = ("NA");
    } elsif ($segNum == 7) {
      @segNames = ("MP", "M2");
    } elsif ($segNum == 8) {
      @segNames = ("NS");
    } else { die "Unexpected segNum $segNum" }
    foreach $seg (@segNames) {
      print join("\t", "${run}_${seg}_cns.fa",
                 "Consensus_${run}_${seg}_cns_threshold_0.5_quality_20",
                 $run, $seg, $gb, $segNum, $fullName) . "\n";
    }' \
| sort \
| grep -Fwf <(ls -1 $fasta_dir)

cd
rm -r $tmpdir
