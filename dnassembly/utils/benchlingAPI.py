import requests
import json

#baseURL = "https://outpacebiotest.benchling.com/api/" #testURL
baseURL = "https://outpacebio.benchling.com/api/" #realURL

#key = #testkey
key =  #realkey

# These are the backend locations of things in benchlingAPI
#part3a = "lib_iH2yDHof" #Define the locations of each folder for each part
#registry_id = "src_oh4knSj0"

# For Test Only!
newseq = {}
#newseq["folderId"] = "lib_1A6kl8LI" #Test, does this value change for each folder within a project?
newseq["isCircular"] = True
newseq["namingStrategy"] = "NEW_IDS"

#Real Registry
newseq["registryId"] = "src_oh4knSj0" #Real
newseq["schemaId"] = "ts_itQ9daT4" #Real
project = {"value": ["sfso_uEaQhrCf"]} #Real

#Test Registry
#newseq["registryId"] = "src_D2ebrtJZ" #Test
#newseq["schemaId"] = "ts_aPuxhTTQ" #Test
#project = {"value": ["sfso_gVlMamQK"]} #Test

#newseq["fields"] = {"Antibiotic Resistance": resistance, "Project": project} #Need to define each field as a dict

class BadRequestException(Exception):
    def __init__(self, message, rv):
        super(BadRequestException, self).__init__(message)
        self.rv = rv

def getBenchling(path, query):
    request = baseURL+'v2/'+path+query
    r = requests.get(request, auth=(key,"")) #real version

    if r.status_code >= 400:
        raise BadRequestException(
            "Server returned status {}. Response:\n{}".format(
                r.status_code, json.dumps(r.json())
            ),
            r,
        )

    return r.json()

def postPartBenchling(bases, name, LOH):
    request = baseURL+'v2/dna-sequences'

    newseq["bases"] = bases
    newseq["name"] = name

    antibioticId = "sfso_3vpWtjgu" #Real
    #antibioticId =  "sfso_0CQFeeLN" #test
    resistance = {"value": antibioticId}
    newseq["fields"] = {"Antibiotic Resistance": resistance, "Project": project}

    partType = LOH.strip()

    if partType in ['1']:
        newseq["folderId"] = 'lib_ZaMHBULI'

    elif partType in ['2a','2b']:
        newseq["folderId"] = 'lib_QBZvgB3q'

    elif partType in ['3a','3b','3c','3d','3e']:
        newseq["folderId"] = 'lib_tzkMsLGC'

    elif partType in ['4a','4b']:
        newseq["folderId"] = 'lib_Y275V89m'

    elif partType in ['5']:
        newseq["folderId"] = 'lib_RcvOpOZC'

    elif partType in ['6']:
        newseq["folderId"] = 'lib_QOR8IQ4Y'

    elif partType in ['7']:
        newseq["folderId"] = 'lib_QmgULI3v'

    else:
        newseq["folderId"] = 'lib_kCnFLwBS' #Send to the parent Stage 1 folder

    r = requests.post(request, json=newseq, auth=(key,""))

    if r.status_code >= 400:
        raise BadRequestException(
            "Server returned status {}. Response:\n{}".format(
                r.status_code, json.dumps(r.json())
            ),
            r,
        )

    return r.json()

def searchSeqBenchling(bases):
    request = baseURL+'v2-experimental/dna-sequences:search-bases'
    query = {"bases": bases}
    r = requests.post(request, json=query, auth=(key,""))

    if r.status_code >= 400:
        raise BadRequestException(
            "Server returned status {}. Response:\n{}".format(
                r.status_code, json.dumps(r.json())
            ),
            r,
        )

    templates = r.json()
    return templates['dnaSequences']

#def main():
    """
    newseq = {}
    newseq["bases"] =
    newseq["folderId"] =
    newseq["isCircular"] = True
    newseq["name"] =
    newseq["registryId"] = registry_id
    newseq["namingStrategy"] = "NEW_IDS"

    BenchlingAPI("dna-sequences")
    """

#if __name__ == "__main__":
    #main()
