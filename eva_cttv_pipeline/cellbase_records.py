import codecs
import json
import urllib.request

from eva_cttv_pipeline import config


def get_curr_response(skip):
    reader = codecs.getreader("utf-8")
    answer = urllib.request.urlopen('http://' + config.HOST + '/cellbase/webservices/rest/v3/hsapiens/feature/clinical/all?source=clinvar&skip=' + str(skip) + '&limit=' + str(config.BATCH_SIZE))
    curr_response = json.load(reader(answer))['response'][0]
    return curr_response


def get_curr_result_list(skip):
    curr_response = get_curr_response(skip)
    # print(str(curr_response['numTotalResults']) + ' ClinVar records in total.')
    curr_result_list = curr_response['result']
    return curr_result_list


def get_curr_result_lists():
    skip = 0
    while True:
        curr_result_list = get_curr_result_list(skip)
        if len(curr_result_list) == 0:
            break
        skip += config.BATCH_SIZE
        yield curr_result_list


def get_records():
    for curr_result_list in get_curr_result_lists():
        for record in curr_result_list:
            yield record
