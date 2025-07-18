#!/usr/bin/env python3

import argparse
import sys
import os
import polars as pl
from Bio import Entrez, SeqIO
from dateutil.parser import parse, ParserError
from unidecode import unidecode


# Load GADM data
gadm_data = None
gadm_countries = {}
gadm_states = {}
gadm_counties = {}
gadm_us_state_lookup = {}

def normalize_characters(text):
    """Normalize Unicode characters to ASCII equivalents using unidecode library."""
    if not text:
        return text
    
    return unidecode(text)

def load_gadm_data(gadm_file):
    """Load GADM data from the TSV file.

    Note: 'gadm' refers to the Global Administrative Areas database (https://gadm.org/).
    """
    if gadm_file is None:
        raise ValueError("gadm_file argument must be provided and cannot be None.")
    global gadm_data, gadm_countries, gadm_states, gadm_counties, gadm_us_state_lookup
    
    if gadm_data is not None:
        return
    
    try:
        # Read the TSV file with polars
        try:
            gadm_data = pl.read_csv(gadm_file, separator='\t', encoding='utf-8')
        except Exception as e:
            raise ValueError(f"Error parsing GADM file '{gadm_file}': {e}")
        
        # Process countries
        countries_df = gadm_data.select(['NAME_0', 'GID_0']).filter(
            pl.col('NAME_0').is_not_null() & (pl.col('NAME_0') != "")
        ).unique(subset=['NAME_0'])

        for row in countries_df.iter_rows(named=True):
            gadm_countries[row['NAME_0'].lower()] = {
                'name': row['NAME_0'],
                'gid': row['GID_0']
            }
        
        # Process all states/provinces (not just US)
        all_states_df = gadm_data.filter(
            pl.col('NAME_0').is_not_null() & 
            (pl.col('NAME_0') != "") &
            pl.col('NAME_1').is_not_null() & 
            (pl.col('NAME_1') != "")
        ).select(['NAME_0', 'NAME_1', 'ISO_1']).unique(subset=['NAME_0', 'NAME_1'])

        for row in all_states_df.iter_rows(named=True):
            country_name = row['NAME_0']
            state_name = row['NAME_1']
            iso_1 = row['ISO_1'] or ""
            
            # Create composite key: "country:state"
            state_key = f"{country_name.lower()}:{state_name.lower()}"
            gadm_states[state_key] = {
                'name': state_name,
                'country': country_name,
                'iso': iso_1
            }
            
            # Also create normalized key for character matching
            state_key_normalized = f"{country_name.lower()}:{normalize_characters(state_name.lower())}"
            if state_key_normalized != state_key:
                gadm_states[state_key_normalized] = {
                    'name': state_name,
                    'country': country_name,
                    'iso': iso_1
                }
            
            # For US states, also add to the US-specific lookup
            if country_name == "United States":
                state_code = iso_1.split('-')[-1] if iso_1.startswith('US-') else ""
                gadm_us_state_lookup[state_name.lower()] = state_name
                if state_code:
                    gadm_us_state_lookup[state_code.lower()] = state_name
        
        # Process US counties
        us_counties_df = gadm_data.filter(
            (pl.col('NAME_0') == 'United States') & 
            pl.col('NAME_1').is_not_null() & 
            (pl.col('NAME_1') != "") &
            pl.col('NAME_2').is_not_null() & 
            (pl.col('NAME_2') != "")
        ).select(['NAME_1', 'NAME_2', 'HASC_2', 'ENGTYPE_2']).unique(subset=['NAME_1', 'NAME_2'])

        for row in us_counties_df.iter_rows(named=True):
            county_key = f"{row['NAME_1'].lower()}:{row['NAME_2'].lower()}"
            gadm_counties[county_key] = {
                'name': row['NAME_2'],
                'state': row['NAME_1'],
                'hasc': row['HASC_2'] or "",
                'engtype': row['ENGTYPE_2'] or ""
            }
                
    except Exception as e:
        print(f"Warning: Could not load GADM data: {e}", file=sys.stderr)

def get_country_info(c_name_str):
    """Gets the country name from GADM data."""
    c_name_str_stripped = c_name_str.strip()
    if not c_name_str_stripped:
        return ""
    
    if c_name_str_stripped.upper() == 'USA':
        return "United States"
    
    # Look up in GADM countries
    if c_name_str_stripped.lower() in gadm_countries:
        return gadm_countries[c_name_str_stripped.lower()]['name']
    
    # Try partial matching for common country name variations
    c_name_lower = c_name_str_stripped.lower()
    for country_key, country_data in gadm_countries.items():
        if c_name_lower in country_key or country_key in c_name_lower:
            return country_data['name']
    
    return ""

def get_state_province_name(s_name_str, country_name):
    """Gets the state/province name using GADM data for any country."""
    s_name_str_stripped = s_name_str.strip()
    if not s_name_str_stripped or not country_name:
        return ""

    s_name_lower = s_name_str_stripped.lower()
    s_name_normalized = normalize_characters(s_name_lower)

    if country_name == "United States":
        # Check if input is a state code
        if s_name_lower in gadm_us_state_lookup:
            return gadm_us_state_lookup[s_name_lower]
        
        return None
    
    # For all countries, use the pre-built state lookup
    state_key = f"{country_name.lower()}:{s_name_lower}"
    if state_key in gadm_states:
        return gadm_states[state_key]['name']
    
    state_key_normalized = f"{country_name.lower()}:{s_name_normalized}"
    if state_key_normalized in gadm_states:
        return gadm_states[state_key_normalized]['name']
    
    # Try partial matching within the same country (both original and normalized)
    for key, state_data in gadm_states.items():
        if key.startswith(f"{country_name.lower()}:"):
            state_name_lower = state_data['name'].lower()
            state_name_normalized = normalize_characters(state_name_lower)
            
            # Check if input matches state name (either original or normalized)
            if (s_name_lower in state_name_lower or state_name_lower in s_name_lower or
                s_name_normalized in state_name_normalized or state_name_normalized in s_name_normalized):
                return state_data['name']
    
    return None

def get_city_canonical(city_n_str, state_name):
    """Get canonical city name with proper capitalization."""
    city_n_str_stripped = city_n_str.strip()
    if not city_n_str_stripped:
        return ""

    # Simple capitalization of each word
    result = ' '.join(word.capitalize() for word in city_n_str_stripped.split())
    return result

def infer_county_from_city(city_name, state_name):
    """Infer county from city and state using GADM data."""
    if not city_name or not state_name:
        return ""
    
    # Handle case where input already contains "County", "Borough", or "Parish"
    base_name = city_name
    county_types = ["County", "Borough", "Parish"]
    for county_type in county_types:
        if city_name.endswith(f" {county_type}"):
            base_name = city_name[:-len(f" {county_type}")].strip()
            break
    
    # Check if this base name (or full name) is actually a county
    county_key = f"{state_name.lower()}:{base_name.lower()}"
    if county_key in gadm_counties:
        # Use the pre-loaded county data - return the proper county name + "County"
        county_data = gadm_counties[county_key]
        county_name = county_data['name']
        return f"{county_name} County"
    
    return ""

def normalize_location(input_string):
    """Parses location strings to COUNTRY/STATE/COUNTY/CITY format."""
    if not input_string or not input_string.strip():
        return "////"

    input_stripped = input_string.strip()
    country_final, state_final, county_final, city_final = "", "", "", ""
    
    # case 1: "USA"
    if ":" not in input_string:
        country_final = get_country_info(input_string)
        # If country is not found in GADM data, return ////
        if not country_final:
            return "////"
        result = f"{country_final}/{state_final}/{county_final}/{city_final}"
        return normalize_characters(result)
    
    # split country and rest of components
    parts = input_string.split(':', 1)
    raw_country_input = parts[0].strip()
    country_final = get_country_info(raw_country_input)
    
    # If country is not found in GADM data, return ////
    if not country_final:
        return "////"
    
    # case 2: "Country:State"
    if "," not in input_string or country_final != "United States":
        state_final = get_state_province_name(parts[1].strip(), country_final) or ""
        result = f"{country_final}/{state_final}/{county_final}/{city_final}"
        return normalize_characters(result)
    
    # case 3: "USA:State,County or City" (only for United States)
    other_components_str = parts[1].strip() if len(parts) > 1 else ""
    other_parts = [p.strip() for p in other_components_str.split(',')] if other_components_str else []
    
    if len(other_parts) < 2:
        result = f"{country_final}/{state_final}/{county_final}/{city_final}"
        return normalize_characters(result)
    
    # "Nan County" - just parse country and state
    if len(other_parts) >= 2 and other_parts[1].strip().lower() == "nan county":
        state_final = get_state_province_name(other_parts[0], country_final) or ""
        county_final = ""
        city_final = ""
        result = f"{country_final}/{state_final}/{county_final}/{city_final}"
        return normalize_characters(result)
    
    # US-specific detailed parsing
    state_from_second = get_state_province_name(other_parts[1], country_final)
    state_from_first = get_state_province_name(other_parts[0], country_final)
    
    if state_from_second:
        # Format: "Country: City, State"
        state_final = state_from_second
        city_final = get_city_canonical(other_parts[0], state_final)
        county_final = infer_county_from_city(other_parts[0], state_final)
    elif state_from_first:
        # Format: "Country: State, City" or "Country: State, County"
        state_final = state_from_first
        # Check if second part is a county
        county_result = infer_county_from_city(other_parts[1], state_final)
        if county_result:
            county_final = county_result
        else:
            city_final = get_city_canonical(other_parts[1], state_final)
    else:
        # Fallback: leave empty
        state_final = ""
        city_final = ""

    result = f"{country_final}/{state_final}/{county_final}/{city_final}"
    return normalize_characters(result)

def get_genbank_data(accessions: list) -> list:
    """Get collection dates and geo_loc_name for a batch of GenBank accessions."""

    try:
        handle = Entrez.efetch(db="nucleotide", id=",".join(accessions), rettype="gb", retmode="text")
        records = SeqIO.parse(handle, "genbank")
        
        genbank_data = []
        for record in records:
            data = {"genbank_acc": record.id}
            source_feature = next((f for f in record.features if f.type == "source"), None)
            if source_feature:
                data["Collection_Date_gb"] = source_feature.qualifiers.get("collection_date", [None])[0]
                data["geo_loc_name_gb"] = source_feature.qualifiers.get("geo_loc_name", [None])[0]
            genbank_data.append(data)
            
        handle.close()
        return genbank_data

    except Exception as e:
        raise RuntimeError(f"An unexpected error occurred while fetching GenBank data: {e}")


def populate_fields_ncbi_avian(metadata_df: pl.LazyFrame, genbank_df: pl.LazyFrame) -> pl.DataFrame:
    """Get updated data from NCBI and join with metadata."""

    # Filter for unique HA segment accessions
    ha_accessions_df = (genbank_df
                        .filter(pl.col("seg") == "HA")
                        .group_by("sra_run")
                        .agg(pl.col("genbank_acc").first())
                    )

    accessions = ha_accessions_df.select("genbank_acc").collect().unique().to_series().to_list()

    # Process accessions in batches to query NCBI
    batch_size = 1000
    all_ncbi_data = []
    for i in range(0, len(accessions), batch_size):        
        batch = accessions[i:i+batch_size]
        all_ncbi_data.extend(get_genbank_data(batch))
        print(f"Processed {len(all_ncbi_data)} of {len(accessions)} records so far...")
        time.sleep(0.5)  # Add a delay to avoid overwhelming the NCBI servers

    # Create a DataFrame with the new data from NCBI
    latest_genbank_df = pl.DataFrame(all_ncbi_data).lazy()

    # Join new NCBI data with the HA accessions
    genbank_with_updates = ha_accessions_df.join(latest_genbank_df, on="genbank_acc", how="left")

    print(f"genbank_with_updates: {genbank_with_updates.collect().shape}")

    # Join with main metadata to update fields, coalescing new and old values
    updated_metadata_df = (metadata_df
                           .join(genbank_with_updates, left_on="Run", right_on="sra_run", how="left")
                           .with_columns(
                               geo_loc_name=pl.coalesce("geo_loc_name_gb", "geo_loc_name"),
                               Collection_Date=pl.coalesce("Collection_Date_gb", "Collection_Date")
                           )
                           .drop(["geo_loc_name_gb", "Collection_Date_gb", "genbank_acc"]))

    return updated_metadata_df.collect()


def normalize_date(date_str: str) -> str:
    """Formats a date string into YYYY-MM-DD or YYYY-MM if only month is available or YYYY if only year is available."""
    if not date_str or not date_str.strip():
        return ""

    s_date = str(date_str).strip().lower()
    if not s_date:
        return ""

    try:
        # Parse the date string
        dt_obj = parse(s_date)
        
        # Pre-compute formatted strings to check against original
        day_formats = [dt_obj.strftime(f).lower() for f in ['%d', '%-d', '%j']]
        month_formats = [dt_obj.strftime(f).lower() for f in ['%m', '%-m', '%b', '%B']]
        
        # Check if components are present in the original string
        has_day = any(fmt in s_date for fmt in day_formats)
        has_month = any(fmt in s_date for fmt in month_formats)
        
        if has_day and has_month:
            return dt_obj.strftime("%Y-%m-%d")
        elif has_month:
            return dt_obj.strftime("%Y-%m")
        else:
            return dt_obj.strftime("%Y")

    except (ParserError, ValueError, TypeError):
        return ""


def normalize_metadata(metadata_df: pl.DataFrame) -> pl.DataFrame:
    """Normalize metadata fields."""
    return (metadata_df
            .with_columns(
                pl.col("Collection_Date").map_elements(normalize_date, return_dtype=pl.String),
                pl.col("geo_loc_name").cast(pl.String).str.strip_chars().map_elements(normalize_location, return_dtype=pl.String),
                pl.col("Run").cast(pl.String).str.strip_chars().alias("Run"),
                pl.col("isolation_source").cast(pl.String).str.strip_chars().alias("isolation_source")
            ))


def process_delimited_metadata(metadata_file: str, genbank_file: str, output_file: str) -> bool:
    """Process delimited metadata file, enrich with NCBI data, and save."""
    schema_overrides = {
        'Collection_Date': pl.String,
        'geo_loc_name': pl.String,
        'Run': pl.String,
        'host': pl.String,
        'isolation_source': pl.String
    }
    infer_separator = lambda file: '\t' if file.lower().endswith('.tsv') else ','
    try:
        genbank_df = pl.scan_csv(
            genbank_file, 
            separator=infer_separator(genbank_file),
            schema_overrides=schema_overrides, 
            null_values=''
        )
        metadata_df = pl.scan_csv(
            metadata_file, 
            separator=infer_separator(metadata_file),
            schema_overrides=schema_overrides, 
            null_values=''
        )
        
        updated_df = populate_fields_ncbi_avian(metadata_df, genbank_df)
        updated_df = normalize_metadata(updated_df)
        updated_df.write_csv(output_file, separator='\t', null_value='',)
        
        print(f"Successfully processed and saved updated metadata to {output_file}")
        return True
    except Exception as e:
        print(f"Error processing delimited file: {e}", file=sys.stderr)
        return False


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Normalize metadata for avian influenza sequences.")
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input metadata file')
    parser.add_argument('-g', '--genbank_file', required=True, help='Path to the GenBank metadata file')
    parser.add_argument('-d', '--gadm_file', required=False, default='assets/gadm_pkg_names.tsv', help='Path to the GADM data file')
    parser.add_argument('-e', '--email', required=False, default="test@test.com", help="Email address for NCBI Entrez API.")
    parser.add_argument('-o', '--output_file', required=True, help='Path to save the normalized metadata.')
    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_args()
    load_gadm_data(args.gadm_file)
    Entrez.email = args.email
    process_delimited_metadata(args.input_file, args.genbank_file, args.output_file)


if __name__ == "__main__":
    main()
