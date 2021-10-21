#! /usr/bin/env python3
from itertools import tee
from typing import Dict, Tuple
from .constants import DNAConstants

# Thank you itertools cookbook!
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."

    # returns 2 independent iterators from a single iterable
    a, b = tee(iterable)
    
    next(b, None)
    return zip(a, b)


def cycle_in_frames(iterable, frame=3):
    """
    Cycle through iterable and yield all possible frames of size frame. Treats iterable as a circular sequence.
    Janky but it works...
    :return:
    """
    assert frame < len(iterable)
    cycle_iterable = iterable * 2

    for cycle_count, derp in enumerate(cycle_iterable):
        if cycle_count < len(iterable):
            yield cycle_iterable[cycle_count: cycle_count + frame]
            cycle_count += 1


def reverse_dict(self, reference_dict) -> Dict[str, str]:
    reversed_dict = {}
    for key, value in reference_dict.items():
        reversed_dict[value] = key

    return reversed_dict

def reverse_complement(sequence) -> str:
    "Get reverse complement of DNA sequence"

    reversed_sequence = ""
    
    for nucleotide in sequence:
        if nucleotide in DNAConstants.DNA_BASEPAIRS.keys():
            reversed_sequence += reversed(DNAConstants.DNA_BASEPAIRS[nucleotide])
        else:
            reversed_sequence += reversed(nucleotide)

    # TODO - check if above statement works
    #return ''.join(reversed([dna_basepairs[a] if a in dna_basepairs.keys() else a for a in sequence]))


# TODO - check if above statement works
# codon_to_res = {codon: res for res, codons in res_to_codons.items() for codon in codons}