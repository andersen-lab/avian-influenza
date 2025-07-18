import argparse
import pandas as pd
import requests
import io
import os
from datetime import date
from Bio import Entrez, SeqIO


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Map GenBank accessions to SRA runs and FASTA files."
    )
    parser.add_argument("metadata_csv", help="Path to the metadata CSV file from SRA.")
    parser.add_argument(
        "fasta_dir", help="Path to the directory containing FASTA files."
    )
    parser.add_argument(
        "output_tsv",
        help="Path to the output TSV file where results will be saved.",
    )
    parser.add_argument(
        "--email",
        required=False,
        default="test@test.com",
        help="Email address for NCBI Entrez API.",
    )
    return parser.parse_args()


def download_ncbi_metadata():
    """Download NCBI virus metadata for H5N1 in the USA since 2023."""
    tax_id = 11320
    today = date.today().strftime("%Y-%m-%d")
    params = {
        "fq": [
            '{!tag=SeqType_s}SeqType_s:("Nucleotide")',
            f"VirusLineageId_ss:({tax_id})",
            f"{{!tag=CollectionDate_dr}}CollectionDate_dr:[2023-01-01T00:00:00.00Z TO {today}T23:59:59.00Z]",
            "{!edismax qf=Serotype_s}H5N1",
            '{!tag=Country_s}Country_s:("USA")',
        ],
        "q": "*:*",
        "cmd": "download",
        "dlfmt": "csv",
        "fl": "genbank_acc:AccVer_s,isolate:Isolate_s,location:CountryFull_s,collected:CollectionDate_s,host:Host_s,bioproject_accession:BioProject_s,strain:Strain_s,serotype:Serotype_s,segment:Segment_s",
        "sort": "id asc",
    }

    base_url = "https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/"
    response = requests.get(base_url, params=params)
    response.raise_for_status()

    ncbi_df = pd.read_csv(io.StringIO(response.text))
    ncbi_df["sample_name_raw"] = ncbi_df["strain"].str.split("/").str[3]
    ncbi_df["sample_name"] = ncbi_df["strain"].str.extract(r"/(\d{2}-\d{6}-\d{3})")
    return ncbi_df


def resolve_samples(df):
    """
    Resolve ambiguous SRA-GenBank mappings by cross-referencing BioSample IDs.
    
    Parameters:
        df (pd.DataFrame): A DataFrame containing SRA metadata with at least the following columns:
            - 'genbank_acc': GenBank accession numbers.
            - 'BioSample': BioSample IDs.
            - 'sra_run': SRA run identifiers.
    
    Returns:
        tuple:
            - correct_matches (pd.DataFrame): A DataFrame of rows where the BioSample ID matches
              between the SRA metadata and GenBank data. Includes resolved mappings.
            - unresolved (pd.DataFrame): A DataFrame of rows that could not be resolved due to
              ambiguous or missing BioSample information.
    
    Resolution Algorithm:
        1. Fetch GenBank records for the given accession numbers in batches.
        2. Extract BioSample IDs from the GenBank records' cross-references.
        3. Merge the extracted BioSample data with the input DataFrame on 'genbank_acc'.
        4. Identify rows where the BioSample ID matches between the input and GenBank data.
        5. Mark rows as unresolved if their 'sra_run' and 'genbank_acc' are not part of the
           resolved matches.
    """
    accessions = df["genbank_acc"].unique()
    if len(accessions) == 0:
        return None, None

    batch_size = 1000
    biosample_data = []

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        try:
            handle = Entrez.efetch(
                db="nucleotide", id=batch.tolist(), rettype="gb", retmode="text"
            )
            for record in SeqIO.parse(handle, "genbank"):
                if hasattr(record, "dbxrefs"):
                    for xref in record.dbxrefs:
                        if xref.startswith("BioSample:"):
                            biosample_id = xref.split(":", 1)[1]
                            biosample_data.append(
                                {
                                    "genbank_acc": record.id,
                                    "genbank_biosample": biosample_id,
                                }
                            )
                            break
            handle.close()
        except Exception as e:
            print(f"An error occurred during Entrez fetch: {e}")

    if not biosample_data:
        print("Could not retrieve BioSample information from GenBank.")
        return None, df

    genbank_biosample_df = pd.DataFrame(biosample_data).drop_duplicates()
    merged_df = pd.merge(df, genbank_biosample_df, on="genbank_acc", how="left")

    # Create a boolean mask for rows where BioSample matches genbank_biosample.
    # This identifies rows that are definitively correct.
    is_correct_row = merged_df["BioSample"] == merged_df["genbank_biosample"]

    # Identify correct matches directly.
    correct_matches = merged_df[is_correct_row].copy()

    # Get unique sra_run and genbank_acc from the correct matches.
    correct_sra_runs = set(correct_matches["sra_run"])
    correct_genbank_accs = set(correct_matches["genbank_acc"])

    # A row is unresolved if its sra_run AND genbank_acc are not in the sets of correct ones.
    # This avoids including rows that are part of a group that has been resolved.
    is_unresolved = ~merged_df["sra_run"].isin(correct_sra_runs) & \
                    ~merged_df["genbank_acc"].isin(correct_genbank_accs)
    
    unresolved = merged_df[is_unresolved]

    return correct_matches, unresolved


def process_sra_metadata(metadata_csv_path):
    """Process SRA metadata to clean up sample names."""
    sra_df = pd.read_csv(metadata_csv_path, usecols=["Run", "Sample Name", "BioSample"])
    sra_df.rename(
        columns={"Run": "sra_run", "Sample Name": "sample_name_raw"}, inplace=True
    )
    sra_df["sample_name"] = sra_df["sample_name_raw"].str.replace(
        r"-original|-repeat2?", "", regex=True
    )
    sra_df["sample_name"] = sra_df["sample_name"].str.replace(
        r"([0-9]{2}-[0-9]{6}-[0-9]{3})-(300|MTM)", r"\1", regex=True
    )
    return sra_df


def create_final_dataframe(sra_df, ncbi_df, fasta_dir):
    """
    Merge SRA and NCBI data, map segments, and filter based on available FASTA files.
    """
    # Merge on raw and manipulated sample names
    exact_match_df = pd.merge(sra_df, ncbi_df, on="sample_name_raw", how="inner")
    sra_fallback = sra_df[~sra_df["sra_run"].isin(exact_match_df["sra_run"])]
    ncbi_fallback = ncbi_df[~ncbi_df["genbank_acc"].isin(exact_match_df["genbank_acc"])]
    fallback_match_df = pd.merge(
        sra_fallback, ncbi_fallback, on="sample_name", how="inner"
    )
    merged_df = pd.concat([exact_match_df, fallback_match_df], ignore_index=True)

    # Identify and resolve faulty SRA assignments
    sra_counts_in_segment = merged_df.groupby(["segment", "sra_run"])[
        "sra_run"
    ].transform("size")
    faulty_sra_assignment = merged_df[sra_counts_in_segment > 1]

    if not faulty_sra_assignment.empty:
        print("Attempting to resolve SRA assignments for segments with multiple SRA runs")
        correct_matches, unresolved = resolve_samples(faulty_sra_assignment)

        print(len(correct_matches), "correct matches found")

        if correct_matches is not None and not correct_matches.empty:
            print('Resolved SRA assignments for segments with multiple SRA runs')
            merged_df = merged_df[
                ~merged_df["sra_run"].isin(faulty_sra_assignment["sra_run"])
            ]
            merged_df = pd.concat([merged_df, correct_matches], ignore_index=True)

        if unresolved is not None and not unresolved.empty:
            print("Number of unresolved SRA assignments:", len(unresolved))
            unresolved.to_csv(
                "unresolved_sra_assignments.tsv", sep="\t", index=False, header=True
            )
    else:
        print("No segments with multiple SRA runs found, proceeding with existing data.")

    # Map segments and filter based on available FASTA files using vectorized operations
    seg_map = {
        1: ["PB2"],
        2: ["PB1"],
        3: ["PA"],
        4: ["HA"],
        5: ["NP"],
        6: ["NA"],
        7: ["MP", "M2"],
        8: ["NS"],
    }
    merged_df["seg_name"] = merged_df["segment"].map(seg_map)
    merged_df.dropna(subset=["seg_name"], inplace=True)

    output_df = merged_df.explode("seg_name")

    output_df["seg_file"] = (
        output_df["sra_run"] + "_" + output_df["seg_name"] + "_cns.fa"
    )

    fasta_files = set(os.listdir(fasta_dir))
    output_df = output_df[output_df["seg_file"].isin(fasta_files)].copy()

    if output_df.empty:
        return pd.DataFrame()

    output_df.rename(columns={"seg_name": "seg"}, inplace=True)
    output_df["seg_seq_name"] = (
        "Consensus_"
        + output_df["sra_run"]
        + "_"
        + output_df["seg"]
        + "_cns_threshold_0.5_quality_20"
    )
    output_df["genbank_seg"] = output_df["segment"]
    output_df["genbank_name"] = output_df["strain"]

    final_cols = [
        "seg_file",
        "seg_seq_name",
        "sra_run",
        "seg",
        "genbank_acc",
        "genbank_seg",
        "genbank_name",
    ]
    output_df = output_df[final_cols]

    output_df.sort_values(by=["seg_file"], inplace=True)
    return output_df


def main():
    """Main function to run the script."""
    args = parse_args()
    ncbi_df = download_ncbi_metadata()
    Entrez.email = args.email
    sra_df = process_sra_metadata(args.metadata_csv)
    output_df = create_final_dataframe(sra_df, ncbi_df, args.fasta_dir)
    if not output_df.empty:
        output_df.to_csv(args.output_tsv, sep="\t", index=False, header=True)
    else:
        print("No matching records found to write to output.")


if __name__ == "__main__":
    main()
