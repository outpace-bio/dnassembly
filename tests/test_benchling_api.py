from dnassembly.utils.benchlingAPI import handle_response, get_dna_type_id,\
    DNA_PART_1_ID, DNA_PART_2_ID, DNA_PART_3_ID, DNA_PART_4_ID, DNA_PART_5_ID,\
    DNA_PART_6_AND_7_ID, get_sequence_folder_id, get_resistance
from typing import List, Tuple
import pytest

@pytest.fixture(scope='class')
def parts_list():
    return [
        ("1", DNA_PART_1_ID),
        ("2b 2c", DNA_PART_2_ID),
        ("3c", DNA_PART_3_ID),
        ("4b", DNA_PART_4_ID),
        ("5", DNA_PART_5_ID),
        ("7", DNA_PART_6_AND_7_ID),
        ("6", DNA_PART_6_AND_7_ID),
        ("", None)
    ]

def test_get_dna_type_id(parts_list: List[Tuple[str, str]]) -> None:

    for parts in parts_list:

        input_parts = parts[0]
        expected_output = parts[1]

        actual_output = get_dna_type_id(input_parts)

        assert actual_output == expected_output

@pytest.fixture(scope='class')
def assembly_types():
    return [
        ("cassette", ),
        ("MC", ),
        ("", ),
        ("1", ),
        ("2b 2c", ),
        ("3c", ),
        ("4b", ),
        ("5", ),
        ("7", ),
        ("6", ),
    ]

def test_get_sequence_folder_id(assembly_types: List[Tuple[str, str]]) -> None:

    for assembly_type in assembly_types:

        input_assembly_type = assembly_type[0]
        expected_output = assembly_type[1]

        actual_output = get_sequence_folder_id(input_assembly_type)

        assert actual_output == expected_output


def test_get_resistance() -> None:

    for assembly_type in assembly_types:

        input_assembly_type = assembly_type[0]
        expected_output = assembly_type[1]

        actual_output = get_sequence_folder_id(input_assembly_type)

        assert actual_output == expected_output