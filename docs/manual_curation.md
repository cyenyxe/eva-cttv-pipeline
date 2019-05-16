# Manual trait curation

The goal is for traits with occurence ≥ 10 to have 100% coverage after the manual curation. For the rest of the
traits, curate as many as possible.

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
Note that the number of columns in the output table is limited to 50, because only a few traits have that many
mappings, and in virtually all cases these mappings are not meaningful. However, having a very large table degrades
the performance of Google Sheets substantially.

```bash
cut -f-50 table_for_manual_curation.tsv | sort -t$'\t' -k2,2rn | xclip -selection clipboard
```

Create a Google Sheets table like this one
[_traits_requiring_curation_](https://docs.google.com/spreadsheets/d/1mb_ZAEwlSTLCQYBWsihxvUGWoy-otaKFq8tIxpJVT0U/)
and paste the result to that table.

## Manual curation criteria
Good mappings must be eyeballed to ensure they are actually good. Alternative mappings for medium or low quality
mappings can be searched for using OLS. If a mapping can't be found in EFO, look for a mapping to a HP, ORDO, or
MONDO trait name. Most HP/ORDO/MONDO terms will also be in EFO but some are not. These can be imported to EFO using
the Webulous submission service.

### Criteria to manually evaluate mapping quality
* Exact string for string matches are _good_
* Slight modifications are _good_ e.g. IRAK4 DEFICIENCY → Immunodeficiency due to interleukin-1 receptor-associated
  kinase-4 deficiency
* Subtype to parent are _good_ e.g ACHROMATOPSIA 3 → Achromatopsia
* Parent to subtype are _bad_ e.g. HEMOCHROMATOSIS → Hemochromatosis type 3
Familial / congenital represented on only one half are _bad_ e.g. Familial renal glycosuria → Renal glycosuria
* Susceptibility on only one half is _bad_ e.g Alcohol dependence, susceptibility to → alcohol dependence
* Early / late onset on only one half is _bad_ e.g. Alzheimer disease, early-onset → Alzheimer's disease

### Unmapped trait names
Trait names that haven't been automatically mapped against any ontology term can also be searched for using OLS. If a
mapping can't be found in EFO, look for a mapping to a HP, ORDO, or MONDO trait name. If these are not already in EFO
they should be imported to EFO using the Webulous submission service.

## Curation workflow
Curation should be done by subsequently applying filters to appropriate columns, then making decisions for the traits
in the filtered selection.

* 1\. **There is a previously assigned mapping for this trait.** All of these are the decisions that we made in the
  past, so we trust them (to an extent). Copy and paste previously used mappings into “Mapping to use”. Then review
  them according to the following steps.
  * 1.1. **The previously assigned mapping is in EFO**
    * 1.1.1. **The previously assigned mapping is in EFO and is exact.** Mark as finished immediately. (It's
      extremely unlikely that a better mapping could be found).
    * 1.1.2. **The previously assigned mapping is in EFO and IS NOT exact.** Review the mappings to see if a better
      (more accurate/specific) mapping is available. Then mark as finished.
  * 1.2. **The previously assigned mapping is not contained in EFO.** We need to either find a mapping which is
    already in EFO, or import these terms into EFO.
    * 1.2.1. **The previously used mapping IS NOT contained in EFO and is exact.** These are good candidates to mark
      as finished and them import in EFO afterwards. However, quickly check whether there are non-exact matches which
      are already in EFO are are as good as exact mappings.
      * E. g. if the exact mapping is “erythrocytosis 6, familial” and not in EFO, but there is an inexact mapping
        “familial erythrocytosis 6” which *is* in EFO, we should use the inexact mapping.
      * If a trait does not have any EFO mappings, it's probably safe to mark it as finished (with subsequent import
        to EFO).
    * 1.2.2. **The previously assigned mapping IS NOT contained in EFO and IS NOT exact.** Similarly to 1.2.1,
      attempt to find an acceptable EFO mapping; if not found, use any acceptable mapping (with subsequent import to
      EFO).
* 2\. **There is no previously assigned mappings for the trait, but exact mappings are available.** Because
  letter-to-letter matches are extremely likely to be correct, we can use them after eyeballing for correctness.
  * 2.1. **The exact mapping in the EFO.** Mark as finished immediately.
  * 2.2. **The exact mapping IS NOT in the EFO.** Similarly to 1.2.1, attempt to find an acceptable EFO mapping; if
    not found, use any acceptable mapping (with subsequent import to EFO).
* 3\. **There are no previously used or exact mappings for the trait.** Curate manually as usual.

### Time-saving options
The new manual workflow can be shortened if necessary, while the quality of the results will be _at least as good as
for the old workflow_ (because we're reusing the results of previous curations):
* All subsections 1.\* involve review of mappings previously selected by ourselves. Because we trust them (to an
  extent), this review can be applied not to all mappings, but only to some (selected on a basis of frequency, or
  just randomly sampled/eyeballed).
* If necessary, section 1 can be skipped completely, i. e. copy-paste previous mappings into “Mapping to use” column,
  but skip the review.
* Sections 2.2 and 3 can only be applied to some variants (e. g. based on frequency), depending on the time
  available.

## Exporting curation results
All curated mappings should be stored in a file named `finished_mappings_curation.tsv`.
 
After that, the following mappings must be written to a single file to be used as input for the evidence string
generation:
* Mappings generated automatically by the trait mapping pipeline and already considered “finished”
  (`automated_trait_mappings.tsv`)
* Eyeballed good quality mappings (`finished_mappings_curation.tsv`)
* Manually curated medium and low quality mappings (`finished_mappings_curation.tsv`)
* New mappings for previously unmapped traits (`finished_mappings_curation.tsv`)

The resulting file must be named `trait_names_to_ontology_mappings.tsv`.
