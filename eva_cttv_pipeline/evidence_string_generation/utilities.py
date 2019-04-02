import argparse
import gzip
import sys
import importlib
import importlib._bootstrap
import importlib.util
import os


def open_file(file_path, mode):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


def get_resource_file(package, resource):
    spec = importlib.util.find_spec(package)
    if spec is None:
        return None
    mod = (sys.modules.get(package) or importlib._bootstrap._load(spec))
    if mod is None or not hasattr(mod, '__file__'):
        return None

    parts = resource.split('/')
    parts.insert(0, os.path.dirname(mod.__file__))
    resource_name = os.path.join(*parts)
    return resource_name


class ArgParser:
    """
    Uses argparse module to parse command line arguments.
    Arguments are used in the pipeline, including input file paths, output path, path to files to
    specify EFO urls to either ignore or alter, and clinical significances that will be allowed to
    generate evidence strings.
    """

    def __init__(self, argv):
        usage = """Generates CTTV evidence strings from ClinVar data and trait mappings."""
        parser = argparse.ArgumentParser(usage)

        parser.add_argument("--clinSig", dest="clinical_significance",
                            help="""Optional. String containing a comma-sparated list with the
                            clinical significances that will be allowed to generate
                            evidence-strings. By default all clinical significances will be
                            considered. Possible tags: 'unknown','untested','non-pathogenic',
                            'probable-non-pathogenic','probable-pathogenic','pathogenic',
                            'drug-response','drug response','histocompatibility','other','benign',
                            'protective','not provided','likely benign','confers sensitivity',
                            'uncertain significance','likely pathogenic',
                            'conflicting data from submitters','risk factor','association' """,
                            default="pathogenic,likely pathogenic,protective,association,risk_factor,affects,drug response")
        parser.add_argument("--ignore", dest="ignore_terms_file",
                            help="""Optional. String containing full path to a txt file containing
                            a list of term urls which will be ignored during batch processing """,
                            default=None)
        parser.add_argument("--adapt", dest="adapt_terms_file",
                            help="""Optional. String containing full path to a txt file containing
                            a list of invalid EFO urls which will be adapted to a general valid url
                             during batch processing """,
                            default=None)
        parser.add_argument("--out", dest="out",
                            help="""String containing the name of the file were
                            results will be stored.""",
                            required=True)

        parser.add_argument("-e", "--efoMapFile", dest="efo_mapping_file",
                            help="Path to file with trait name to url mappings", required=True)
        parser.add_argument("-g", "--snp2GeneFile", dest="snp_2_gene_file",
                            help="Path to file with RS id to ensembl gene ID and consequence "
                                 "mappings", required=True)
        parser.add_argument("-j", dest="json_file", help="File containing Clinvar records json "
                                                         "strings in the format of documents in "
                                                         "Cellbase. One record per line.")
        parser.add_argument('--ot-schema', help='OpenTargets schema JSON', required=True)

        args = parser.parse_args(args=argv[1:])

        self.clinical_significance = args.clinical_significance
        self.ignore_terms_file = args.ignore_terms_file
        self.adapt_terms_file = args.adapt_terms_file
        self.out = args.out
        self.efo_mapping_file = args.efo_mapping_file
        self.snp_2_gene_file = args.snp_2_gene_file
        self.json_file = args.json_file
        self.ot_schema = args.ot_schema


def check_dir_exists_create(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
