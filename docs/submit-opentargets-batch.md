Please follow the instructions below in order to create and submit an OpenTargets batch using ClinVar data. Additional diagrams and background explanations can be found in [this presentation](https://docs.google.com/presentation/d/1nai1dvtfow4RkolyITcymXAsQqEwPJ8pUPcgjLDCntM).

# Preparing batch processing folder
Given the year and month the batch is to be released on, run the following command to create the appropriate folders:

```bash
mkdir batch-[yyyy]-[mm] # referred to as 'batch root folder' in next steps
cd batch-[yyyy]-[mm]
mkdir clinvar gene_mapping trait_mapping evidence_strings
```

# Preparing ClinVar data
We can find the ClinVar release associated to an Ensembl release in its [sources page](http://www.ensembl.org/info/genome/variation/species/sources_documentation.html).

Two files need to be downloaded from the ClinVar FTP into the _clinvar_ subfolder:
* Full release XML: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/
* Variant summary: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/

## Transforming ClinVar XML to JSON
The ClinVar XML schema version can be obtained by inspecting the XML file header. If this version is not supported by the ClinVar XML parser, we have to generate the JAXB binding classes to be able to parse the XML. We can do it running the following command:

```bash
python bin/update_clinvar_schema.py -i [path_to_batch_root_folder]/clinvar/ClinVarFullRelease_[yyyy]-[mm].xml.gz -j clinvar-xml-parser/src/main/java
```

A new Java package should have been generated in the directory _clinvar-xml-parser/src/main/java/uk/ac/ebi/eva/clinvar/model_. A Pull Request to merge this code into the main repository should be created, updating the test data to the newer ClinVar version.

Once the new model binding classes have been generated, and the java parser has been built (see [README](https://github.com/EBIvariation/eva-cttv-pipeline/blob/master/README.md#building-java-clinvar-parser)) we can transform the ClinVar XML file running the command below.

```bash
java -jar [path_to_parser_jar_file_containing_dependencies] -i [path_to_batch_root_folder]/clinvar/ClinVarFullRelease_[yyyy]-[mm].xml.gz -o [path_to_batch_root_folder]/clinvar
```

A file named _clinvar.json.gz_ will be created in the output directory.

## Filter Clinvar JSON file
Clinvar JSON file is then filtered, excluding records without an allowed clinical significance.

```bash
python bin/clinvar_jsons/extract_pathogenic_and_likely_pathogenic_variants.py -i [path_to_batch_root_folder]/clinvar/clinvar.json.gz -o [path_to_batch_root_folder]/clinvar/clinvar.filtered.json.gz
```

# Gene and consequence type mappings
A file containing coordinate, allele, variant type and ID information needs to be given to OpenTargets so they can map the ClinVar records to an affected gene and consequence type.

```bash
python bin/gene_mapping/gene_map_coords.py -i [path_to_batch_root_folder]/clinvar/variant_summary_[yyyy]-[mm].txt.gz -o [path_to_batch_root_folder]/gene_mapping/clinvar_[yyyy]-[mm]_coords.tsv.gz
```

The columns in the tsv output file are:
* Chromosome
* Start position
* Stop position
* Reference allele
* Alternate allele
* Strand
* Structural variant type (calculated in script)
* rs ID
* RCV ID (ID for a ClinVar record)
* NCBI gene ID
* sv ID
* Variant type (as specified by ClinVar)

This output file must be uploaded to the [OpenTargets Google Cloud Storage](https://console.cloud.google.com/storage/browser/otar012-eva/). They will fed it into their mapping pipeline, which consists of VEP and some internal custom software. The tsv file returned by OpenTargets has the columns:

* Variant (rs ID, sv ID, coordinate and alleles, or RCV)
* na (unused)
* Ensembl gene ID
* HGNC gene symbol
* Functional consequence
* na (unused)

# Trait mapping pipeline

```bash
python bin/trait_mapping.py -i [path_to_batch_root_folder]/clinvar/clinvar.filtered.json.gz -o [path_to_batch_root_folder]/trait_mapping/automated_trait_mappings.tsv -c [path_to_batch_root_folder]/trait_mapping/traits_requiring_curation.tsv
```

The trait mapping pipeline aims to find a term in EFO corresponding to each trait name from ClinVar. You can find a diagram of the whole workflow [here](https://docs.google.com/presentation/d/1nai1dvtfow4RkolyITcymXAsQqEwPJ8pUPcgjLDCntM/edit#slide=id.g24b2b34015_0_531).

## Querying Zooma
Zooma is first queried using the trait name.

* If there are any high confidence mappings from EFO then they are output for use and no further processing is taken for that trait name.
* If there are lower confidence mappings from EFO then these are output in a separate curation file, and the process stops here for that trait name.
* If there are high confidence mappings not from EFO then their ontology IDs are used as input for querying OxO (see below)

## Querying OxO
[OxO](http://www.ebi.ac.uk/spot/oxo/) is a tool which finds ontology cross references for provided identifiers.

* Any EFO mappings found within a distance of 1 step are output for use, no further processing for this trait name.
* Any EFO mappings within a distance greater than 1 step are output for curation.

## Output
The output consists of 2 files.

### Automated trait mappings
The file _automated_trait_mappings.tsv_ contains the following columns:

* Trait name from ClinVar
* Ontology URI
* Ontology term label

### Manual curation required
The file _traits_requiring_curation.tsv_ contains the entries that require manual curation, including any trait name with no mappings. It contains the following columns:

* Trait name from ClinVar
* Frequency of this trait in ClinVar
* Zooma mappings:
  * URI
  * Ontology label
  * Confidence
  * Datasource/annotator
* OxO mappings:
  * URI
  * Ontology label
  * Confidence
  * Distance

# Manual curation

## Extract information about previous mappings

At this step, mappings produced by the pipeline on the previous iteration (including automated and manual) are 
downloaded to be used to aid the manual curation process.

```bash
# Download the latest eva_clinvar release from FTP
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/eva_clinvar.txt \
  | cut -f4-5 | sort -u > previous_mappings.tsv
```

## Create the final table for manual curation

```bash
python3 bin/trait_mapping/create_table_for_manual_curation.py \
  --traits-for-curation traits_requiring_curation.tsv \
  --previous-mappings previous_mappings.tsv \
  --output table_for_manual_curation.tsv
```

## Sort and export to Google Sheets

Note that the number of columns in the output table is limited to 50, because only a few traits have that many mappings, and in virtually all cases these mappings are not meaningful. However, having a very large table degrades the performance of Google Sheets substantially.  

```bash
cut -f-50 table_for_manual_curation.tsv | sort -t$'\t' -k2,2rn | xclip -selection clipboard
# The paste into Google Sheets
```

Traits with occurence ≥ 10 must have 100% coverage after the manual curation.

To do this manually it helps to open an excel sheet like this one [_traits_requiring_curation_](https://docs.google.com/spreadsheets/d/1ghb7VB6AHMfYKGbIT4L894H9TsPlel7ujWLPUHkPMY4/edit#gid=1707979701) and sort the mappings by frequency. Follow the steps below for assessing the mapping quality of as many traits as possible with frequency greater than 2.
All curated mappings should be stored in a file named _finished_mappings_curation.tsv_.

## Assessing existing mapping quality
Good mappings must be eyeballed to ensure they are actually good. Alternative mappings for medium or low quality mappings can be searched for using OLS. If a mapping cant be found in EFO, look for a mapping to a HP or ORDO trait name. Most HP and ORDO terms will also be in EFO but there are some that arent. These can be imported to EFO using the Webulous submission service. 

The criteria to manually evaluate mapping quality is:

* Exact string for string matches are _good_
* Slight modifications are _good_ e.g. IRAK4 DEFICIENCY -> Immunodeficiency due to interleukin-1  receptor-associated kinase-4 deficiency
* Subtype to parent are _good_ e.g ACHROMATOPSIA 3 -> Achromatopsia
* Parent to subtype are _bad_ e.g. HEMOCHROMATOSIS -> Hemochromatosis type 3
Familial / congenital represented on only one half are _bad_ e.g. Familial renal glycosuria -> Renal glycosuria 
* Susceptibility on only one half is _bad_ e.g Alcohol dependence, susceptibility to -> alcohol dependence
* Early / late onset on only one half is _bad_ e.g. Alzheimer disease, early-onset -> Alzheimer's disease

## Unmapped trait names
Trait names that haven't been automatically mapped against any ontology term can also be searched for using OLS. If a mapping cant be found in EFO, look for a mapping to a HP or ORDO trait name. If these are not already in EFO they should be imported to EFO using the Webulous submission service. 

## Putting together a file with finished mappings after curation
The following mappings must be written to a single file to be used as input for the evidence string generation:

* Mappings generated automatically by the trait mapping pipeline and already considered "finished" (_automated_trait_mappings.tsv_)
* Eyeballed good quality mappings (_finished_mappings_curation.tsv_)
* Manually curated medium and low quality mappings (_finished_mappings_curation.tsv_)
* New mappings for previously unmapped traits (_finished_mappings_curation.tsv_)

The resulting file must be named _trait_names_to_ontology_mappings.tsv_.

## Adding novel trait names to EFO
Traits remaining unmapped or poorly mapped can be submitted to EFO if a suitable parent term is available. Any new terms added will be picked up by Zooma in the next iteration of the trait mapping pipeline.
Novel traits can be submitted to EFO using the [Webulous templates](https://www.ebi.ac.uk/efo/webulous/) **Add EFO disease** and **Import HP term**. Open a new Google spreadsheet and connect with the server using the Webulous Add-on. 
Open a new git issue with EFO to review and import these novel trait names e.g. [https://github.com/EBISPOT/efo/issues/223](https://github.com/EBISPOT/efo/issues/223)

# Generating evidence strings
## Installing pipeline package
Evidence string generation is carried out using the `eva_cttv_pipeline` package in this repository. This package can be installed using `python setup.py install`, but if also developing the pipeline then `python setup.py develop` should be used. For more information, please see the README in the root of the repository.

## Preparing local schema
OpenTargets have a [JSON schema](https://github.com/opentargets/json_schema) used to validate submitted data. Validation of generated evidence strings is carried out during generation. To fetch schema, issue the following commands:

```bash
VERSION=1.6.0  # substitute the version recommended by OT
wget -q https://raw.githubusercontent.com/opentargets/json_schema/$VERSION/opentargets.json
```

The location of `opentargets.json` is not important because it will be passed as a parameter.

## Executing
In order to generate the evidence strings, run the following command:

```bash
python bin/evidence_string_generation.py \
  --out [path_to_batch_root_folder]/evidence_strings/ \
  -e [path_to_batch_root_folder]/trait_mapping/trait_names_to_ontology_mappings.tsv z
  -g [path_to_batch_root_folder]/gene_mapping/clinvar_[yyyy]-[mm]_coords.out \
  -j [path_to_batch_root_folder]/clinvar/clinvar.filtered.json.gz \
  --ot-schema opentargets.json
```

This outputs multiple files, including the file of evidence strings (evidence_strings.json) for submitting to OpenTargets and the file of trait mappings for submitting to Zooma (eva_clinvar.txt).

# Submitting evidence strings
The evidence string file (evidence_strings.json) should be uploaded to the [OpenTargets Google Cloud Storage](https://console.cloud.google.com/storage/browser/otar012-eva/) and be named in the format `cttv012-[dd]-[mm]-[yyyy].json.gz` (e.g. `cttv012-12-06-2017.json.gz`).

More details can be found on [OpenTargets Github wiki](https://github.com/opentargets/data-providers-docs/wiki/Data-(Evidence-Strings)-Submission-Process).

# Submitting feedback to Zooma

## Trait mappings
The file containing trait mappings (eva_clinvar.txt) must be uploaded to the EVA FTP. It will be placed in a folder structure like [this one](ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/2017/07/21/) and shall also be available through [the link to the latest version](ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/).

## ClinVar xrefs
Each ClinVar record is associated with one or more traits. When submitting the data to OpenTargets, the trait needs to be specified as an ontological term present in [the Experimental Factor Ontology (EFO)](http://www.ebi.ac.uk/efo/).

Traits in ClinVar may specify ontological terms for the trait name, but often not to a term present in EFO. Cross-references (xrefs) to ontological terms can help in finding a good EFO term for a trait, whether an xref is in EFO itself or an xref can be used as a starting point for searching for a term in EFO later on.

Within the ClinVar records are ontology xRefs that link the trait name to a controlled vocabulary. We submit these to [Zooma](http://www.ebi.ac.uk/spot/zooma/) under the evidence handle "ClinVar_xRefs".

The mappings from the ClinVar trait name to the specified ontology xrefs are parsed from the ClinVar JSON file into a tsv suitable for submitting to Zooma, excluding any traits which already have mappings in Zooma from trusted data sources (EVA, Open Targets, GWAS, Uniprot). In order to do so, execute the following command:

```bash
python bin/clinvar_jsons/traits_to_zooma_format.py -i [path_to_batch_root_folder]/clinvar/clinvar.filtered.json.gz -o [path_to_output_file]
```

Once the output file is created, the symbolic link in the EVA FTP is updated to point to the new file. The current version of this file on the FTP can be found in ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/clinvar_xrefs.txt
