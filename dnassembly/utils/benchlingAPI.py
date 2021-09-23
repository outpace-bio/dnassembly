import requests
import json
import os
from .http_status import HTTPStatus
from typing import List, Tuple

# Benchling tenant environmental variables
BASE_URL = os.environ.get('BENCHLING_URL')
API_VERSION = os.environ.get('BENCHLING_API_VERSION')
API_KEY = os.environ.get('BENCHLING_API_KEY')
PROJECT_ID = os.environ.get('BENCHLING_PROJ_ID')
SEQ_REGISTRY_ID = os.environ.get('SEQ_REGISTRY_ID')
SEQ_SCHEMA_ID = os.environ.get('SEQ_SCHEMA_ID')
SEQ_NAMING_STRAT = os.environ.get('SEQ_NAMING_STRAT')
SEQ_PARENT_DIR_ID = os.environ.get('PART_PARENT_DIR_ID') #lib_kCnFLwBS
PART_1_SEQ_DIR_ID = os.environ.get('PART_1_SEQ_DIR_ID') #lib_ZaMHBULI
PART_2_SEQ_DIR_ID = os.environ.get('PART_2_SEQ_DIR_ID') #lib_QBZvgB3q
PART_3_SEQ_DIR_ID = os.environ.get('PART_3_SEQ_DIR_ID') #lib_tzkMsLGC
PART_4_SEQ_DIR_ID = os.environ.get('PART_4_SEQ_DIR_ID') #lib_Y275V89m
PART_5_SEQ_DIR_ID = os.environ.get('PART_5_SEQ_DIR_ID') #lib_RcvOpOZC
PART_6_SEQ_DIR_ID = os.environ.get('PART_6_SEQ_DIR_ID') #lib_QOR8IQ4Y
PART_7_SEQ_DIR_ID = os.environ.get('PART_7_SEQ_DIR_ID') #lib_QmgULI3v


PART_REGISTRY_ID = os.environ.get('PART_REGISTRY_ID')
PART_SCHEMA_ID = os.environ.get('PART_SCHEMA_ID')
PART_NAMING_STRAT = os.environ.get('PART_NAMING_STRAT')
DNA_PART_FOLDER_ID = os.environ.get('DNA_PART_FOLDER_ID')
DNA_PART_1_ID = os.environ.get('DNA_PART_1_ID') #sfso_ZIYN45Gl
DNA_PART_2_ID = os.environ.get('DNA_PART_2_ID') #sfso_xvdSPEIX
DNA_PART_3_ID = os.environ.get('DNA_PART_3_ID') #sfso_a5vEtqm5
DNA_PART_4_ID = os.environ.get('DNA_PART_4_ID') #sfso_U3jQSIjd
DNA_PART_5_ID = os.environ.get('DNA_PART_5_ID') #sfso_gyH8kxAF
DNA_PART_6_AND_7_ID = os.environ.get('DNA_PART_6_AND_7_ID') #sfso_eDKuXLe5

CARB = os.environ.get('CARB')
KAN = os.environ.get('KAN')

#Test Registry
#newSeq["folderId"] = "lib_1A6kl8LI" #Test, does this value change for each folder within a project?
#newSeq["registryId"] = "src_D2ebrtJZ" #Test
#newSeq["schemaId"] = "ts_aPuxhTTQ" #Test
#project = {"value": ["sfso_gVlMamQK"]} #Test

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

def get_dna_type_id(parts: str) -> str:
    part_type = parts.strip()

    if part_type in ['1']:
        return DNA_PART_1_ID
    elif part_type in ['2a','2b']:
        return DNA_PART_2_ID
    elif part_type in ['3a','3b','3c','3d','3e']:
        return DNA_PART_3_ID
    elif part_type in ['4a','4b']:
        return DNA_PART_4_ID
    elif part_type in ['5']:
        return DNA_PART_5_ID
    elif part_type in ['6', '7']:
        return DNA_PART_6_AND_7_ID

#TODO - change assembly type to enum class and share class between files
def get_sequence_folder_id(assembly_type: str) -> str:
    if assembly_type == 'cassette': #Set location to the parent Stage 2 folder
        return 'lib_lIpZ86uz'
    elif assembly_type == 'MC': #Set location to the parent Stage 3 folder
        return 'lib_hWRdTDLG'
    else: #Set location to somewhere in Stage 1 folder
        part_type = assembly_type.strip()

        if part_type in ['1']:
            return PART_1_SEQ_DIR_ID
        elif part_type in ['2a','2b']:
            return PART_2_SEQ_DIR_ID
        elif part_type in ['3a','3b','3c','3d','3e']:
            return PART_3_SEQ_DIR_ID
        elif part_type in ['4a','4b']:
            return PART_4_SEQ_DIR_ID
        elif part_type in ['5']:
            return PART_5_SEQ_DIR_ID
        elif part_type in ['6']:
            return PART_6_SEQ_DIR_ID
        elif part_type in ['7']:
            return PART_7_SEQ_DIR_ID
        else:
            return SEQ_PARENT_DIR_ID #Send to the parent Stage 1 folder

def get_resistance(assembly_type: str) -> str:
    #Test
    #antibioticId =  "sfso_0CQFeeLN" #test
    #resistance = {"value": antibioticId}
    #newSeq["fields"] = {"Antibiotic Resistance": resistance, "Project": project}

    #Real
    if assembly_type == 'cassette': #Set the antibiotic to Carb for Stage 2
        return CARB
    else: #Set the antibiotic to Kan for both Stage 1 and Stage 3
        return KAN

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
    res = requests.get(request, auth=(API_KEY,""))

    handle_response(res)

    return res.json()

def postSeqBenchling(squence_bases: str, name: str, assembly_type: str, assembled_from: str = '', assembled_from_id: str = ''):
    request = f"{BASE_URL}{API_VERSION}/dna-sequences"

    sequence_folder_id = get_sequence_folder_id(assembly_type)
    resistance = get_resistance(assembly_type)

    newSeq = {
        "isCircular": True,
        "namingStrategy": SEQ_NAMING_STRAT,
        "registryId": SEQ_REGISTRY_ID,
        "schemaId": SEQ_SCHEMA_ID,
        "folderId": sequence_folder_id,
        "bases": squence_bases,
        "name": name,
        "fields": {
            "Antibiotic Resistance": { "value": resistance },
            "Project": project,
            "MoClo_Assembled_From": { "value": assembled_from},
            "MoClo_Assembled_From_seqID": {"value": assembled_from_id}
        }
    }

    res = requests.post(request, json=newSeq, auth=(API_KEY,""))

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
    
    DNA_type_id = get_dna_type_id(part_type)

    newPart = {
        "isCircular": False,
        "namingStrategy": PART_NAMING_STRAT,
        "registryId": PART_REGISTRY_ID,
        "schemaId": PART_SCHEMA_ID,
        "folderId": DNA_PART_FOLDER_ID,
        "bases": squence_bases,
        "name": name,
        "fields": {
            "Project": { 
                "value": [
                    PROJECT_ID
                ]
            },
            "DNA Type": {"value": DNA_type_id}
        }
    }   

    res = requests.post(request, json=newPart, auth=(API_KEY,""))

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
    
    query = {
        "bases": squence_bases,
        "registryId": SEQ_REGISTRY_ID,
        "schemaId": SEQ_SCHEMA_ID
    }

    res = requests.post(request, json=query, auth=(API_KEY,""))

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
    res = requests.post(request, json=annotateList, auth=(API_KEY,""))

    handle_response(res)
    return