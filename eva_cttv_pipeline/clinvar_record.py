import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime
import http.client


__author__ = 'Javier Lopez: javild@gmail.com'


def get_rcv_to_rsnsv_mapping(variant_summary_file):

    rcv_to_rs = {}
    rcv_to_nsv = {}

    print('Loading mapping RCV->rs/nsv')
    fdr = open(variant_summary_file, "r")
    fdr.readline()
    for line in fdr:
        parts = line.split('\t')
        rcv_list = parts[8].split(';')
        if parts[6] != '-' and parts[6] != '-1' and parts[6] != '':
            for rcv_id in rcv_list:
                rcv_to_rs[rcv_id] = 'rs' + parts[6]
        if parts[7] != '-' and parts[7] != '-1' and parts[7] != '':
            for rcv_id in rcv_list:
                rcv_to_nsv[rcv_id] = parts[7]
    fdr.close()
    print(' Done.')

    return rcv_to_rs, rcv_to_nsv


class ClinvarRecord(dict):

    cached_symbol_2_ensembl = {}

    score_map = {
        "CLASSIFIED_BY_SINGLE_SUBMITTER": 1,
        "NOT_CLASSIFIED_BY_SUBMITTER": None,
        "CLASSIFIED_BY_MULTIPLE_SUBMITTERS": 2,
        "REVIEWED_BY_EXPERT_PANEL": 3,
        "REVIEWED_BY_PROFESSIONAL_SOCIETY": 4
    }

    def __init__(self, a_dictionary=None):
        if a_dictionary is None:
            dict.__init__(self)
        else:
            dict.__init__(self, a_dictionary)

    @property
    def gene_id(self):
        j = 0
        measure = self['referenceClinVarAssertion']['measureSet']['measure']
        found = False
        while j < len(measure) and not found:
            attribute_set = measure[j]['attributeSet']
            i = 0
            while i < len(attribute_set) and not attribute_set[i]['attribute']['type'].startswith('HGVS'):
                i += 1
            found = (i < len(attribute_set))
            j += 1

        if found:
            return attribute_set[i]['attribute']['value'][:9]
        else:
            return self['referenceClinVarAssertion']['measureSet']['measure'][0]['name'][0]['elementValue']['value']

    @property
    def ensembl_id(self):
        global ensembl_json
        j = 0
        measure = self['referenceClinVarAssertion']['measureSet']['measure']
        while j < len(measure):
            if 'measureRelationship' in measure[j]:
                measureRelationship = measure[j]['measureRelationship']
                l = 0
                while l < len(measureRelationship):
                    if 'symbol' in measureRelationship[l]:
                        symbol = measureRelationship[l]['symbol']
                        i = 0
                        while i < len(symbol):
                            try:
                                return ClinvarRecord.cached_symbol_2_ensembl[symbol[i]['elementValue']['value']]
                            except KeyError:
                                not_solved = True
                                while not_solved:
                                    try:
                                        url = 'http://rest.ensembl.org/lookup/symbol/homo_sapiens/' + symbol[i]['elementValue']['value'] + '?content-type=application/json'
                                        rawreply = urllib.request.urlopen(url).read()
                                        ensembl_json = json.loads(rawreply.decode())
                                        not_solved = False
                                    except urllib.error.HTTPError as e:
                                        if e.code == 400:
                                            print('WARNING: Bad request code returned from ENSEMBL rest.')
                                            print(' ClinVar accession: ' + self.get_acc())
                                            print(' Gene symbol: ' + symbol[i]['elementValue']['value'])
                                            print(' Error: ')
                                            print(e)
                                            print(' Returning None.')
                                            ClinvarRecord.cached_symbol_2_ensembl[
                                                symbol[i]['elementValue']['value']] = None
                                            return None
                                        else:
                                            time.sleep(0.05)
                                    except http.client.BadStatusLine:
                                        time.sleep(3)

                                if len(ensembl_json) > 0:
                                    if ensembl_json['id'] and (ensembl_json['object_type'] == 'Gene'):
                                        ClinvarRecord.cached_symbol_2_ensembl[symbol[i]['elementValue']['value']] = \
                                            ensembl_json['id']
                                        return ensembl_json['id']
                                    else:
                                        print("WARNING at clinvar_record.py: ENSEMBL's REST API returned an unexpected json object")
                                        print(" Queried gene symbol: " + symbol[i]['elementValue']['value'])
                                        print(" ENSEMBL's API response:")
                                        print(ensembl_json)
                                        print("Exiting.")
                                        sys.exit(1)
                                else:
                                    ClinvarRecord.cached_symbol_2_ensembl[symbol[i]['elementValue']['value']] = None
                            i += 1
                    l += 1
            j += 1

        return None

    @property
    def date(self):
        return datetime.fromtimestamp(self['referenceClinVarAssertion']['dateLastUpdated'] / 1000).isoformat()

    @property
    def score(self):
        return self.score_map[self['referenceClinVarAssertion']['clinicalSignificance']['reviewStatus']]

    @property
    def accession(self):
        return self['referenceClinVarAssertion']['clinVarAccession']['acc']

    @property
    def traits(self):
        trait_list = []
        for trait_record in self['referenceClinVarAssertion']['traitSet']['trait']:
            trait_list.append([])
            for name_record in trait_record['name']:
                if (name_record['elementValue'][
                        'type'] == 'Preferred'):  # First trait name in the list will always be the "Preferred" one
                    trait_list[-1] = [name_record['elementValue']['value']] + trait_list[-1]
                else:
                    trait_list[-1].append(name_record['elementValue']['value'])

        return trait_list

    @property
    def trait_pubmed_refs(self):
        pubmed_refs_list = []
        for trait_record in self['referenceClinVarAssertion']['traitSet']['trait']:
            pubmed_refs_list.append([])
            if 'citation' in trait_record:
                for citation_record in trait_record['citation']:
                    if ('id' in citation_record) and citation_record['id'] is not None:
                        if citation_record['id']['source'] == 'PubMed':
                            pubmed_refs_list[-1].append(int(citation_record['id']['value']))

        return pubmed_refs_list

    @property
    def observed_pubmed_refs(self):
        pubmedrefsList = []
        if 'observedIn' in self['referenceClinVarAssertion']:
            for observedInRecord in self['referenceClinVarAssertion']['observedIn']:
                for observedDataRecord in observedInRecord['observedData']:
                    if 'citation' in observedDataRecord:
                        for citationRecord in observedDataRecord['citation']:
                            if ('id' in citationRecord) and citationRecord['id'] is not None:
                                if citationRecord['id']['source'] == 'PubMed':
                                    pubmedrefsList.append(int(citationRecord['id']['value']))
        return pubmedrefsList

    @property
    def measure_set_pubmed_refs(self):
        pubmed_refs_list = []
        for measure_record in self['referenceClinVarAssertion']['measureSet']['measure']:
            if 'citation' in measure_record:
                for ciration_record in measure_record['citation']:
                    if ('id' in ciration_record) and ciration_record['id'] is not None:
                        if ciration_record['id']['source'] == 'PubMed':
                            pubmed_refs_list.append(int(ciration_record['id']['value']))
        return pubmed_refs_list

    @property
    def trait_refs_list(self):
        return [ ['http://europepmc.org/abstract/MED/' + str(ref) for ref in refList]
                 for refList in self.trait_pubmed_refs]

    @property
    def observed_refs_list(self):
        return ['http://europepmc.org/abstract/MED/' + str(ref) for ref in self.observed_pubmed_refs]

    @property
    def measure_set_refs_list(self):
        return ['http://europepmc.org/abstract/MED/' + str(ref) for ref in self.measure_set_pubmed_refs]

    @property
    def hgvs(self):
        hgvs_list = []
        for measure_record in self['referenceClinVarAssertion']['measureSet']['measure']:
            for attribute_set_record in measure_record['attributeSet']:
                if attribute_set_record['attribute']['type'].startswith('HGVS'):
                    hgvs_list.append(attribute_set_record['attribute']['value'])

        return hgvs_list

    @property
    def clinical_significance(self):
        return self['referenceClinVarAssertion']['clinicalSignificance']['description'].lower()

    def get_rs(self, rcv_to_rs):
        try:
            return rcv_to_rs[self.accession]
        except KeyError:
            return None

    def get_nsv(self, rcv_to_nsv):
        try:
            return rcv_to_nsv[self.accession]
        except KeyError:
            return None

    def get_main_consequence_types(self, consequence_type_dict, rcv_to_rs):

        new_rs_id = self.get_rs(rcv_to_rs)
        if new_rs_id is not None and (new_rs_id in consequence_type_dict):
            return consequence_type_dict[new_rs_id]
        else:
            return None

    @property
    def variant_type(self):
        return self['referenceClinVarAssertion']['measureSet']['measure'][0]['type']

    @property
    def allele_origins(self):
        allele_origins = set()
        for clinvar_assetion_document in self['clinVarAssertion']:
            for observed_in_document in clinvar_assetion_document['observedIn']:
                allele_origins.add(observed_in_document['sample']['origin'])

        return list(allele_origins)
