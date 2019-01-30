#!/bin/sh

clinvar_filtered=$1
automated=$2
manual=$3

if [ -z "$clinvar_filtered" ] || [ -z "$automated" ] || [ -z "$manual" ]
then
	echo "Please provide all 3 parameters [clinvar.filtered.json.gz] [automated_trait_mappings.tsv] [traits_requiring_curation.tsv]"
	exit 1
fi

##PROCESS JSON FILE

#uncompress file
gunzip -k ${clinvar_filtered}

#extracts only the value type (Preferred or Alternate) and the trait name
jq .clinvarSet.referenceClinVarAssertion.traitSet.trait[].name[].elementValue clinvar.filtered.json > trait_names.txt

#extracts only the trait names which are preferred
grep -B 1 "Preferred" trait_names.txt > trait_names_preferred.txt

# delete 'see cases', 'not provided' and 'not specified' lines
# delete 'value' attribute
# delete 'type' lines (which are separate from the lines containing the trait names)
# delete lines with only --
# delete trailing ",
sed -i -e '/"see cases"/Id' -e '/"not provided"/Id' -e '/"not specified"/Id' -e 's/"value": "//g' -e '/"type": "Preferred"/d' -e '/^--/d' -e 's/",//g' trait_names_preferred.txt

#remove spaces at the beginning of the line, sort and delete duplicates
awk '{$1=$1;print}' trait_names_preferred.txt | sort | uniq > trait_names_preferred_sorted_uniq.txt


##PROCESS OUTPUT FILES

#delete header for automated
tail -n +2 $automated > automated_no_header.tsv

#concatenate both files
cat automated_no_header.tsv $manual > concatenated.tsv

#leave only first column
cut -f 1 concatenated.tsv > concatenated_terms.tsv

#delete line with the pattern between forward slashes
sed -i -e '/see cases/Id' -e '/not provided/Id' -e '/not specified/Id' concatenated_terms.tsv

#sort and remove any duplicate lines
sort concatenated_terms.tsv | uniq > concatenated_sorted_uniq.tsv


## CHECKS

#check terms in the initial json and not in the output files
fgrep -i -v -f concatenated_sorted_uniq.tsv trait_names_preferred_sorted_uniq.txt > result_terms_in_json_not_in_pipeline_output.txt

#check terms in output files not in initial json (must be zero) or check lines format
fgrep -i -v -f trait_names_preferred_sorted_uniq.txt concatenated_sorted_uniq.tsv > result_terms_in_pipeline_output_not_in_json.txt
