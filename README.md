## README ##

[![Build Status](https://travis-ci.com/EBIvariation/eva-cttv-pipeline.svg?branch=master)](https://travis-ci.com/EBIvariation/eva-cttv-pipeline)
[![Coverage Status](https://coveralls.io/repos/github/EBIvariation/eva-cttv-pipeline/badge.svg?branch=master)](https://coveralls.io/github/EBIvariation/eva-cttv-pipeline?branch=master)

OpenTargets aims to provide evidence on the biological validity of therapeutic targets and provide an initial
assessment of the likely effectiveness of pharmacological intervention on these targets, using genome-scale
experiments and analysis.

For every OpenTargets release, the EVA curates ClinVar records (variants) and submits them to OpenTargets as
“evidence strings”, JSON strings describing:

* Genes the variant affects
* Functional consequence of the variant on the gene
* Traits associated with the variant, e.g. “parkinson disease” or “age-related macular degeneration”
* Other information about the variant and source, such as related publications

Available documentation:
* [Build instructions](docs/build.md)
* [How to submit an OpenTargets batch](docs/submit-opentargets-batch.md)
  + [Manual curation](docs/manual_curation.md)
