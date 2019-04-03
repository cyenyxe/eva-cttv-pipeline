## README ##

[![Build Status](https://travis-ci.com/EBIvariation/eva-cttv-pipeline.svg?branch=master)](https://travis-ci.com/EBIvariation/eva-cttv-pipeline)
[![Coverage Status](https://coveralls.io/repos/github/EBIvariation/eva-cttv-pipeline/badge.svg?branch=master)](https://coveralls.io/github/EBIvariation/eva-cttv-pipeline?branch=master)


Python version required: 3.5

Installing python 3.5 (optional)
-------
In case you have a newer version and need to run the pipeline you can follow the next steps to install python 3.5 without replacing you default python.
1. Download python 3.5 "wget https://www.python.org/ftp/python/3.5.5/Python-3.5.5.tgz"
2. "sudo tar xzf Python-3.5.5.tgz"
3. "cd Python-3.5.5"
4. "sudo apt-get install zlib1g-dev"
5. "sudo ./configure --enable-optimizations --with-zlib=/usr/include"
6. "sudo make altinstall"
7. "python3.5 -V"

Building java Clinvar parser
-------
1. "cd clinvar-xml-parser"
2. "mvn package"
3. Two jar files will be generated in the 'target' directory, one of them including all the dependencies.  

Regenerating test data
-------
All the test does (for the moment) is checking that parsing 10 records from the XML will
(1) not crash and (2) provide 10 records after parsing. So to regenerate test data, we just
have to extract any 10 records (can just be the first 10 records) from the ClinVar XML file:

```bash
CLINVAR_RELEASE="2019-01"  # set the correct one
zcat ClinVarFullRelease_${CLINVAR_RELEASE}.xml.gz \
  | awk 'BEGIN {RS="</ClinVarSet>\n\n"; ORS=RS} {print} NR==10 {exit}' \
  > ClinvarExample.xml
echo "</ReleaseSet>" >> ClinvarExample.xml
gzip -c <ClinvarExample.xml >ClinvarExample.xml.gz
```

Eyeball input & output files to ensure that the ClinVar format has not changed sufficiently
enough to render this snippet invalid. Then put the generated files into
`clinvar-xml-parser/src/test/resources/` directory.

Deploying local OLS installation
-------

During the preparation of 2019_04 release, which had to be synchronized with EFO v3, OLS had to be deployed
locally because the production deployment of OLS on www.ebi.ac.uk/ols only supported EFO v2 at the time.
This can be done using the following command:

```bash
sudo docker run -p 8080:8080 simonjupp/efo3-ols:3.4.0
```

To use the local deployment, uncomment the configuration section at the top of `/eva_cttv_pipeline/trait_mapping/ols.py`
to specify the URL of the local installation. If you have deployed OLS on the different machine than the one you're 
using to run the pipeline, substitute the correct IP address of the machine where the OLS installation is deployed.

Please contact the semantic data integration team at [SPOT](https://www.ebi.ac.uk/about/spot-team) if you have questions
about local OLS installation.


Building python pipeline and (optional) Setting up virtual environment
-------

1. "git clone --recursive git@github.com:EBIvariation/eva-cttv-pipeline.git"
2. "cd eva-cttv-pipeline"
3. [OPTIONAL] "virtualenv -p python3.5 venv"
4. [OPTIONAL] "source venv/bin/activate" ("venv/bin/deactivate" to deactivate virtualenv)
5. pip install -r requirements.txt
6. And then one of:
   * To install: "python3 setup.py install"
   * To install to develop: "python3 setup.py develop"
   * To build a source distribution: "python3 setup.py sdist"

Postprocessing traits for manual curation
-------

### Extract information about previous mappings

At this step, mappings produced by the pipeline on the previous iteration (including automated and manual) are 
downloaded to be used to aid the manual curation process.

```bash
# Download the latest eva_clinvar release from FTP
wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/eva/ClinVar/latest/eva_clinvar.txt \
  | cut -f4-5 | sort -u > previous_mappings.tsv
```

### Create the final table for manual curation

```bash
python3 bin/trait_mapping/create_table_for_manual_curation.py \
  --traits-for-curation traits_requiring_curation.tsv \
  --previous-mappings previous_mappings.tsv \
  --final-table-for-curation table_for_manual_curation.tsv
```

### Sort and export to Google Sheets

Note that the number of columns in the output table is limited to 50, because only a few traits have that many mappings,
and in virtually all cases these mappings are not meaningful. However, having a very large table degrades the 
performance of Google Sheets substantially.  

```bash
cut -f-50 table_for_manual_curation.tsv | sort -t$'\t' -k2,2rn | xclip -selection clipboard
# The paste into Google Sheets
```

Usage
-------

Please see the [GitHub wiki](https://github.com/EBIvariation/eva-cttv-pipeline/wiki/How-to-submit-an-OpenTargets-batch) for usage


Tests
-------

You can run all tests with: `python setup.py test`
