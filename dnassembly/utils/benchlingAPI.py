import requests
import json
import os
from .http_status import HTTPStatus
from typing import List, Tuple

#BASE_URL = "https://outpacebiotest.benchling.com/api/" #testURL
#BASE_URL = "https://outpacebio.benchling.com/api/" #realURL
BASE_URL = os.environ.get('BENCHLING_URL')

API_VERSION = os.environ.get('BENCHLING_API_VERSION')

#key = #testkey
key =  os.environ.get('BENCHLING_API_KEY') #realkey

#Define the common parameters for new DNA plasmid sequences
newSeq = {}
newSeq["isCircular"] = True
newSeq["namingStrategy"] = os.environ.get('SEQ_NAMING_STRAT')

#Real Registry
newSeq["registryId"] = os.environ.get('SEQ_REGISTRY_ID') #Outpace Registry
newSeq["schemaId"] = os.environ.get('SEQ_SCHEMA_ID') #ID code for Plasmid
project = {"value": [os.environ.get('BENCHLING_PROJ_ID')]} #MoClo Project
CARB = os.environ.get('CARB') #Real
KAN = os.environ.get('KAN') #Real

#Test Registry
#newSeq["folderId"] = "lib_1A6kl8LI" #Test, does this value change for each folder within a project?
#newSeq["registryId"] = "src_D2ebrtJZ" #Test
#newSeq["schemaId"] = "ts_aPuxhTTQ" #Test
#project = {"value": ["sfso_gVlMamQK"]} #Test

#Define the common parameters for new DNA Parts
newPart = {}
newPart["isCircular"] = False
newPart["namingStrategy"] = os.environ.get('PART_NAMING_STRAT')

#Real Registry
newPart["registryId"] = os.environ.get('PART_REGISTRY_ID') #Outpace Registry
newPart["schemaId"] = os.environ.get('PART_SCHEMA_ID') #ID code for DNA Part

#TODO - populate with test registry
#newPart["registryId"] = "" #Outpace Registry
#newPart["schemaId"] = "" #ID code for DNA Part

class BadRequestException(Exception):
    #TODO - figure out rv type
    def __init__(self, message: str, rv):
        super(BadRequestException, self).__init__(message)
        self.rv = rv

def handle_response(res: requests.Response) -> requests.Response:
    """
    Catches response codes indicating a fault and raises a
    BadRequestException, otherwise no action is taken.

    Args:
        res (requests.Response): Response from a HTTP/HTTPS request

    Raises:
        BadRequestException: [description]

    Returns:
        requests.Response: Response from executed HTTP/HTTPS request
    """
    if res.status_code >= HTTPStatus.BAD_REQUEST:
        raise BadRequestException(
            "Server returned status {}. Response:\n{}".format(
                res.status_code, json.dumps(res.json())
            ),
            res,
        )

    return res

def getBenchling(api_endpoint: str, query: str):
    """
    Sends a GET request to the Benchling tenant

    Args:
        api_endpoint (str): endpoint to query
        query (str): query parameters to add to the specified endpoint

    Returns:
        [a]: Results of the query
    """
    request = f"{BASE_URL}{API_VERSION}/{api_endpoint}{query}"
    res = requests.get(request, auth=(key,""))

    handle_response(res)

    return res.json()

def postSeqBenchling(squence_bases: str, name: str, assembly_type: str, assembledFrom: str = '', assembledFromID: str = ''):
    request = f"{BASE_URL}{API_VERSION}/dna-sequences"

    newSeq["bases"] = squence_bases
    newSeq["name"] = name

    #Test
    #antibioticId =  "sfso_0CQFeeLN" #test
    #resistance = {"value": antibioticId}
    #newSeq["fields"] = {"Antibiotic Resistance": resistance, "Project": project}

    #Real
    if assembly_type == 'cassette': #Set the antibiotic to Carb for Stage 2
        resistance = { "value": CARB }
    else: #Set the antibiotic to Kan for both Stage 1 and Stage 3
        resistance = { "value": KAN }

    MoClo_Assembled_From = {"value": assembledFrom}
    MoClo_Assembled_From_seqID = {"value": assembledFromID}

    #Set the fields
    newSeq["fields"] = {
        "Antibiotic Resistance": resistance,
        "Project": project,
        "MoClo_Assembled_From": MoClo_Assembled_From,
        "MoClo_Assembled_From_seqID": MoClo_Assembled_From_seqID
    }

    if assembly_type == 'cassette': #Set location to the parent Stage 2 folder
        newSeq["folderId"] = 'lib_lIpZ86uz'

    elif assembly_type == 'MC': #Set location to the parent Stage 3 folder
        newSeq["folderId"] = 'lib_hWRdTDLG'

    else: #Set location to somewhere in Stage 1 folder
        partType = assembly_type.strip()

        if partType in ['1']:
            newSeq["folderId"] = 'lib_ZaMHBULI'

        elif partType in ['2a','2b']:
            newSeq["folderId"] = 'lib_QBZvgB3q'

        elif partType in ['3a','3b','3c','3d','3e']:
            newSeq["folderId"] = 'lib_tzkMsLGC'

        elif partType in ['4a','4b']:
            newSeq["folderId"] = 'lib_Y275V89m'

        elif partType in ['5']:
            newSeq["folderId"] = 'lib_RcvOpOZC'

        elif partType in ['6']:
            newSeq["folderId"] = 'lib_QOR8IQ4Y'

        elif partType in ['7']:
            newSeq["folderId"] = 'lib_QmgULI3v'

        else:
            newSeq["folderId"] = 'lib_kCnFLwBS' #Send to the parent Stage 1 folder

    res = requests.post(request, json=newSeq, auth=(key,""))

    handle_response(res)

    return res.json()

def postPartBenchling(squence_bases: str, name: str, part_type: str):
    """
    Sends a POST request to the Benchling tenant to create a new
    part

    Args:
        squence_bases (str): [description]
        name (str): [description]
        part_type (str): [description]

    Returns:
        [type]: [description]
    """

    request = f"{BASE_URL}{API_VERSION}/dna-sequences"

    newPart["bases"] = squence_bases
    newPart["name"] = name

    #TODO - have folderId as environment variable
    newPart["folderId"] = 'lib_R5pTYwpd' #Sends part to DNA Part Folder

    if part_type.strip() in ['1']:
        DNA_type = {"value": 'sfso_ZIYN45Gl'}

    elif part_type.strip() in ['2a','2b']:
        DNA_type = {"value": 'sfso_xvdSPEIX'}

    elif part_type.strip() in ['3a','3b','3c','3d','3e']:
        DNA_type = {"value": 'sfso_a5vEtqm5'}

    elif part_type.strip() in ['4a','4b']:
        DNA_type = {"value": "sfso_U3jQSIjd"}

    elif part_type.strip() in ['5']:
        DNA_type = {"value": 'sfso_gyH8kxAF'}

    elif part_type.strip() in ['6','7']:
        DNA_type = {"value": 'sfso_eDKuXLe5'}

    newPart["fields"] = {
        "Project": project,
        "DNA Type": DNA_type
    }

    res = requests.post(request, json=newPart, auth=(key,""))

    handle_response(res)

    return res.json()

def searchSeqBenchling(squence_bases: str) -> List[dict]:
    """
    Searches Benchling tenant for DNA sequence entitites that contain the 
    specified DNA sequence bases.

    Args:
        squence_bases (str): Bases of a DNA sequence to search for.

    Returns:
        List[dict]: A filtered list of DNA sequences that contain the provided bases.
    """
    request = f"{BASE_URL}{API_VERSION}-beta/dna-sequences:search-bases"
    
    #TODO - have registryId and schemaId as environment variable
    query = {
        "bases": squence_bases,
        "registryId": "src_oh4knSj0", #real
        "schemaId": "ts_itQ9daT4"
    }

    #Test
    #"registryId": "src_D2ebrtJZ",
    #"schemaId": "ts_aPuxhTTQ"}

    res = requests.post(request, json=query, auth=(key,""))

    handle_response(res)

    templates = res.json()
    return templates['dnaSequences']

# Does this return none?
def annotatePartBenchling(sequence_ids: List[str]):
    """
    Sends a list of DNA sequence IDs to Benchling tenant, which will initiate 
    an annotation of DNA sequence parts.

    Args:
        sequence_ids (List[str]): list of sequence of IDs to annotate
    """
    request = f"{BASE_URL}{API_VERSION}/dna-sequences:autofill-parts"

    annotateList = { "dnaSequenceIds": sequence_ids }
    res = requests.post(request, json=annotateList, auth=(key,""))

    handle_response(res)
    return