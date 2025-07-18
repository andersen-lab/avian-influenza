# Consensus Sequences for U.S. H5N1 Clade 2.3.4.4b

## 📊 Repository Overview

This repository aims to provide consensus sequences, variant calls and depth information for the SRA data associated with BioProjects listed below. The repository checks for new data every 24 hours and updates the consensus sequences, variant calls, depth information,  demixed milk samples and associated metadata accordingly. Additionally, the repository updates mapping of consensus genomes to the respective GenBank sequences by sample name every 24 hours.

### BioProjects

| BioProject | Description |
| --- | --- |
| **[PRJNA1102327](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1102327)** | U.S. H5N1 clade 2.3.4.4b genotype B3.13 immediate releases related to emergence and spread in dairy cattle. |
| **[PRJNA1122849](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1122849)** | U.S. H5N1 clade 2.3.4.4b immediate releases related to peridomestic animals. |
| **[PRJNA1134696](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1134696)** | Viral sequencing of US dairy milk. |
| **[PRJNA1219588](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1219588)** | U.S. H5N1 clade 2.3.4.4b genotype D1.1 immediate releases related to emergence in dairy cattle. |
| **[PRJNA1207547](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1207547)** | U.S. Influenza A Surveillance in Wildlife. |
| **[PRJNA980729](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA980729)** | US HPAI Sequencing. |

### 📈 Processed Data Provided

| **Data Type** | **Description** | **Location** |
|---------------|-----------------|--------------|
| **Consensus Sequences** | Consensus builds from SRA data (min depth: 1, quality: 20, consensus threshold: 50%). | [./fasta/](./fasta) |
| **Variant Calls** | Variant calls results from corresponding SRA data. | [./variants/](./variants) |
| **Depth Information** | Depth information for each SRA dataset. | [./depth/](./depth) |
| **Milk Surveillance** | Estimates of lineage abundance from sources like “Bulk milk tanks” using [Freyja](https://github.com/andersen-lab/Freyja). | [./demixed/](./demixed) |
| **Genotype Information** | Genotyping processed using [GenoFLU](https://github.com/USDA-VS/GenoFLU). | [./metadata/genoflu_results.tsv](./metadata/genoflu_results.tsv) |
| **SRA Metadata** | Metadata provided from SRA sources. | [./metadata/SraRunTable_automated.csv](./metadata/SraRunTable_automated.csv) |
| **Genome to GenBank Mapping** | Links consensus genomes to their respective GenBank sequences by sample name. | [./metadata/genbank_mapping.tsv](./metadata/genbank_mapping.tsv) |
| **Merged SRA and Genbank Metadata** | Metadata provided from both SRA and GenBank with location and date normalization | [./metadata/SraRunTable_automated_normalized.tsv](./metadata/SraRunTable_automated_normalized.tsv) |

For a NextStrain-style formatted version of the genomes and metadata, please see [moncla-lab/avian-flu-USDA-cattle/](https://github.com/moncla-lab/avian-flu-USDA-cattle/).

### Pipeline and Reference Details

The data processing pipeline is available in [andersen-lab/flusra](https://github.com/andersen-lab/flusra).

All data generated from **23rd May 2024** uses the GenBank genome [A/cattle/Texas/24-008749-003/2024(H5N1)](https://www.ncbi.nlm.nih.gov/nuccore/?term=A/cattle/Texas/24-008749-003/2024) as a reference. The reference genome can be found in [./reference/](./reference). Settings include a minimum depth of 5, minimum quality of 20, and a consensus threshold at 75%.

> [!NOTE] 
> Prior to **8th May 2025**, consensus genomes were called with a minimum depth of 1, and a consensus threshold at 50%. Since then we have increased the minimum depth to 5, and the consensus threshold to 75%. This change was made to improve the quality of the consensus genomes. Most SRAs have been reprocessed with these new settings. The SRAs that have not been reprocessed have `cns_threshold_0.5` in the fasta headers.

> [!NOTE]
> Prior to **23rd May 2024**, consensus genomes for 8 segments were generated using `EPI_ISL_19032063` (source: GISAID) as a reference. These were produced with [iVar v1.4.2](https://github.com/andersen-lab/ivar) using the settings (min depth: 1, quality: 20 consensus threshold: 50%).

## 📖 Data Usage

We invite the scientific community to utilize and scrutinize this data to enhance overall quality. For queries or feedback, please contact us.

Please refer to the [NCBI usage policies]( https://www.ncbi.nlm.nih.gov/home/about/policies/) for more details.

---

We gratefully acknowledge the authors, originating and submitting laboratory of the sequences from GISAID's EpiFlu™ Database we used as references for our genome assemblies. The [list](./acknowledgements/gisaid_acknowledge_table_assemby_reference_sequences.xls) is provided in [./acknowledgements](./acknowledgements).
