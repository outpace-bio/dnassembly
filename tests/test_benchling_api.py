from dnassembly.utils.benchlingAPI import handle_response, get_dna_type_id, get_sequence_folder_id,\
    get_resistance, getBenchling, postSeqBenchling, postPartBenchling, searchSeqBenchling, annotatePartBenchling
from typing import List, Tuple
import pytest

@pytest.fixture(scope='class')
def parts_list():
    return [
        ()
    ]

def test_get_dna_type_id(parts_list: List[Tuple[str, str]]) -> None:

    for parts in parts_list:

        input_parts = parts[0]
        expected_output = parts[1]

        dna_type_id = get_dna_type_id(input_parts)

        assert dna_type_id == expected_output