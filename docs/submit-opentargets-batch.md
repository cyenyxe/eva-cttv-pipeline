# How to submit an OpenTargets batch

Please follow the instructions below in order to create and submit an OpenTargets batch using ClinVar data. Additional
diagrams and background explanations can be found in [this
presentation](https://docs.google.com/presentation/d/1nai1dvtfow4RkolyITcymXAsQqEwPJ8pUPcgjLDCntM).

Before starting the process, follow the [Build instructions](build.md). In particular, you'll need to install Python 3.5
(if you don't have it already), build the Java ClinVar parser, build and install the Python pipeline.



## Step 1. Preparing batch processing folder

The working directory for the processing is `/nfs/production3/eva/opentargets`

Given the year and month the batch is to be released on, run the following command to create the appropriate folders:

```bash
# Referred to as 'batch root folder' in next steps
mkdir batch-[yyyy]-[mm]
cd batch-[yyyy]-[mm]
mkdir clinvar gene_mapping trait_mapping evidence_strings
```



## Step 2. Preparing ClinVar data

Each OpenTargets release is synchronised with a certan Ensembl release version. The specific version is announced in the
e-mail which they send a few weeks before the data submission deadline. Based on Ensembl version, we can find the
ClinVar release associated to an Ensembl release in its [sources
page](http://www.ensembl.org/info/genome/variation/species/sources_documentation.html).

Two files need to be downloaded from the ClinVar FTP into the `clinvar` subfolder:
* Full release XML: `ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/`
* Variant summary: `ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/`

### 2.1. Updating ClinVar schema version (if necessary)
Schema of ClinVar XML files changes from time to time. The schema version can be obtained by inspecting the XML file
header (it should be in the very first line). The current supported version is **1.55**. If the version changes, we have
to regenerate the JAXB binding classes to be able to parse the XML. It can be done using the following command:

```bash
python bin/update_clinvar_schema.py \
  -i [path_to_batch_root_folder]/clinvar/ClinVarFullRelease_[yyyy]-[mm].xml.gz \
  -j clinvar-xml-parser/src/main/java
```

A new Java package should have been generated in the directory
`clinvar-xml-parser/src/main/java/uk/ac/ebi/eva/clinvar/model`.

With each schema version change, test data must be updated as well. See details in [Build
instructions](build.md#regenerating-test-data)

Create a pull request to merge this code into the main repository. It must contain both the updated schema version as
well as the new test data.

After a schema update, you'll also need to rebuild Java parser (see [Build
instructions](build.md#building-java-clinvar-parser))

### 2.2. Converting ClinVar files

Transform the ClinVar XML file:

```bash
java -jar [path_to_parser_jar_file_containing_dependencies] \
  -i [path_to_batch_root_folder]/clinvar/ClinVarFullRelease_[yyyy]-[mm].xml.gz \
  -o [path_to_batch_root_folder]/clinvar
```

A file named `clinvar.json.gz` will be created in the output directory.

### 2.3. Filter Clinvar JSON file
Clinvar JSON file is then filtered, extracting only records with allowed levels of clinical significance:

```bash
python bin/clinvar_jsons/extract_pathogenic_and_likely_pathogenic_variants.py \
  -i [path_to_batch_root_folder]/clinvar/clinvar.json.gz \
  -o [path_to_batch_root_folder]/clinvar/clinvar.filtered.json.gz
```



## Step 3. Gene and consequence type mappings
A file containing coordinate, allele, variant type and ID information needs to be given to OpenTargets so they can map
the ClinVar records to an affected gene and consequence type.

```bash
python bin/gene_mapping/gene_map_coords.py \
  -i [path_to_batch_root_folder]/clinvar/variant_summary_[yyyy]-[mm].txt.gz \
  -o [path_to_batch_root_folder]/gene_mapping/clinvar_[yyyy]-[mm]_coords.tsv.gz
```

The columns in the TSV output file are:
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

This output file must be uploaded to the [OpenTargets Google Cloud
Storage](https://console.cloud.google.com/storage/browser/otar012-eva/). They will fed it into their mapping pipeline,
which consists of VEP and some internal custom software. The tsv file returned by OpenTargets has the columns:

* Variant (rs ID, sv ID, coordinate and alleles, or RCV)
* na (unused)
* Ensembl gene ID
* HGNC gene symbol
* Functional consequence
* na (unused)

After uploading the file to the Cloud Storage bucket, you also need to notify the OpenTargets team via e-mail so that
they can start their mapping pipeline (they won't receive an automatic notification). The reply from OpenTargets will
take some time (usually around a week); however, it doesn't block most of the subsequent steps, which can be done while
waiting for the reply. The only step which _is_ blocked is Step 6, “Evidence string generation”.



## Step 4. Trait mapping pipeline

```bash
python bin/trait_mapping.py \
  -i [path_to_batch_root_folder]/clinvar/clinvar.filtered.json.gz \
  -o [path_to_batch_root_folder]/trait_mapping/automated_trait_mappings.tsv \
  -c [path_to_batch_root_folder]/trait_mapping/traits_requiring_curation.tsv
```

The trait mapping pipeline aims to find a term in EFO corresponding to each trait name from ClinVar. You can find a
diagram of the whole workflow
[here](https://docs.google.com/presentation/d/1nai1dvtfow4RkolyITcymXAsQqEwPJ8pUPcgjLDCntM/edit#slide=id.g24b2b34015_0_531).

### Querying ZOOMA
ZOOMA is first queried using the trait name.

* If there are any high confidence mappings from EFO then they are output for use and no further processing is taken for
  that trait name.
* If there are lower confidence mappings from EFO then these are output in a separate curation file, and the process
  stops here for that trait name.
* If there are high confidence mappings not from EFO then their ontology IDs are used as input for querying OxO (see
  below)

### Querying OxO
[OxO](http://www.ebi.ac.uk/spot/oxo/) is a tool which finds ontology cross references for provided identifiers.

* Any EFO mappings found within a distance of 1 step are output for use, no further processing for this trait name.
* Any EFO mappings within a distance greater than 1 step are output for curation.

### Output
The output consists of 2 files.

#### Automated trait mappings
The file `automated_trait_mappings.tsv` contains the following columns:

* Trait name from ClinVar
* Ontology URI
* Ontology term label

#### Manual curation required
The file `traits_requiring_curation.tsv` contains the entries that require manual curation, including any trait name
with no mappings. It contains the following columns:

* Trait name from ClinVar
* Frequency of this trait in ClinVar
* ZOOMA mappings:
  * URI
  * Ontology label
  * Confidence
  * Datasource/annotator
* OxO mappings:
  * URI
  * Ontology label
  * Confidence
  * Distance



## Step 5. Manual curation

See separate protocol, [Manual curation](manual_curation.md).



## Step 6. Generating evidence strings

### Preparing local schema
OpenTargets have a [JSON schema](https://github.com/opentargets/json_schema) used to validate submitted data. Validation
of generated evidence strings is carried out during generation. To fetch schema, issue the following commands. `VERSION`
needs to be filled with the version number recommended by the OpenTargets in their announcement e-mail.

```bash
VERSION=1.6.0
wget -q https://raw.githubusercontent.com/opentargets/json_schema/$VERSION/opentargets.json
```

File `opentargets.json` should be saved in the batch root directory.

### Generating evidence strings
In order to generate the evidence strings, run the following command. `clinvar_[yyyy]-[mm]_coords.out` file should be
provided by the OpenTargets after processing the file submitted to them on Step 3, “Gene and consequence type mappings”.

```bash
python bin/evidence_string_generation.py \
  --out [path_to_batch_root_folder]/evidence_strings/ \
  -e [path_to_batch_root_folder]/trait_mapping/trait_names_to_ontology_mappings.tsv \
  -g [path_to_batch_root_folder]/gene_mapping/clinvar_[yyyy]-[mm]_coords.out \
  -j [path_to_batch_root_folder]/clinvar/clinvar.filtered.json.gz \
  --ot-schema [path_to_batch_root_folder]/opentargets.json
```

This outputs multiple files, including the file of evidence strings (`evidence_strings.json`) for submitting to
OpenTargets and the file of trait mappings for submitting to ZOOMA (`eva_clinvar.txt`).

After the evidence strings have been generated, summary metrics need to be updated in the Google Sheets
[table](https://docs.google.com/spreadsheets/d/1g_4tHNWP4VIikH7Jb0ui5aNr0PiFgvscZYOe69g191k/).



## Step 7. Submitting evidence strings
The evidence string file (`evidence_strings.json`) should be uploaded to the [OpenTargets Google Cloud
Storage](https://console.cloud.google.com/storage/browser/otar012-eva/) and be named in the format
`cttv012-[dd]-[mm]-[yyyy].json.gz` (e.g. `cttv012-12-06-2017.json.gz`).

More details can be found on [OpenTargets Github
wiki](https://github.com/opentargets/data-providers-docs/wiki/Data-(Evidence-Strings)-Submission-Process).



## Step 8. Submitting feedback to ZOOMA

### Trait mappings
The file containing trait mappings (`eva_clinvar.txt`) must be uploaded to the EVA FTP. It will be placed in a folder
structure like [this one](ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/2017/07/21/) and shall also be available through
[the link to the latest version](ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/).

### ClinVar xrefs
Each ClinVar record is associated with one or more traits. When submitting the data to OpenTargets, the trait needs to
be specified as an ontological term present in [the Experimental Factor Ontology (EFO)](http://www.ebi.ac.uk/efo/).

Traits in ClinVar may specify ontological terms for the trait name, but often not to a term present in EFO.
Cross-references (xrefs) to ontological terms can help in finding a good EFO term for a trait, whether an xref is in EFO
itself or an xref can be used as a starting point for searching for a term in EFO later on.

Within the ClinVar records are ontology xRefs that link the trait name to a controlled vocabulary. We submit these to
[ZOOMA](http://www.ebi.ac.uk/spot/zooma/) under the evidence handle “ClinVar_xRefs”.

The mappings from the ClinVar trait name to the specified ontology xrefs are parsed from the ClinVar JSON file into a
tsv suitable for submitting to ZOOMA, excluding any traits which already have mappings in ZOOMA from trusted data
sources (EVA, Open Targets, GWAS, Uniprot). In order to do so, execute the following command:

```bash
python bin/clinvar_jsons/traits_to_zooma_format.py \
  -i [path_to_batch_root_folder]/clinvar/clinvar.filtered.json.gz \
  -o [path_to_output_file]
```

Once the output file is created, the symbolic link in the EVA FTP is updated to point to the new file. The current
version of this file on the FTP can be found in
`ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/clinvar_xrefs.txt`.



## Step 9. Adding novel trait names to EFO
Traits remaining unmapped or poorly mapped can be submitted to EFO if a suitable parent term is available. Any new terms
added will be picked up by ZOOMA in the next iteration of the trait mapping pipeline.

Novel traits can be submitted to EFO using the [Webulous templates](https://www.ebi.ac.uk/efo/webulous/) **Add EFO
disease** and **Import HP term**. Open a new Google spreadsheet and connect with the server using the Webulous Add-on.

Open a new git issue with EFO to review and import these novel trait names e.g.
[https://github.com/EBISPOT/efo/issues/223](https://github.com/EBISPOT/efo/issues/223)
