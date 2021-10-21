from abc import ABC
from typing import Dict, List
from .helper_classes import Const
from .utils import reverse_dict

class RestrictionEnzymes(Const):

    BBSI_ANNOTATION_F = {
        'TCAA': '1',
        'ACGA': '2a',
        'AACG': '2b',
        'CACC': '3a',
        'TTCT': '3b',
        'AGCA': '3c',
        'AGGC': '3d',
        'TCCA': '3e',
        'ATCC': '4a',
        'GTAA': '4b',
        'GCTA': '5',
        'ACTC': '6',
        'ATTG': '7'
    }

    BBSI_ANNOTATION_R = {
        'ACGA': '1',
        'AACG': '2a',
        'CACC': '2b',
        'TTCT': '3a',
        'AGCA': '3b',
        'AGGC': '3c',
        'TCCA': '3d',
        'ATCC': '3e',
        'GTAA': '4a',
        'GCTA': '4b',
        'ACTC': '5',
        'ATTG': '6',
        'TCAA': '7'
    }

    BSMBI_ANNOTATION_F = {
        'ACTC': 'LS-1',
        'TTCT': '1-2',
        'TCCA': '2-3',
        'ATCC': '3-RE',
        'ACAG': 'RE-LS'
    }

    BSMBI_ANNOTATION_R = {
        'TTCT': 'LS-1',
        'TCCA': '1-2',
        'ATCC': '2-3',
        'ACAG': '3-RE',
        'ACTC': 'RE-LS'
    }


class RestrictionEnzymeBase(ABC):

    def __init__(self, forward: Dict[str, str], reverse: Dict[str, str]) -> None:
        self.__forward = forward
        self.__reverse = reverse
        self.__overhangs_5prime = reverse_dict(forward)
        self.__overhangs_3prime = reverse_dict(reverse)

    @property
    def forward_parts(self) -> List[str]:
        return list(self.__forward.values())

    @property
    def reverse_parts(self) -> List[str]:
        return list(self.__reverse.values())

    #TODO - should there be a return indicating an issue other than empty str?
    def get_5prime_overhang(self, part: str) -> str:
        return self.__overhangs_5prime.get(part, "")

    #TODO - should there be a return indicating an issue other than empty str?
    def get_3prime_overhang(self, part: str) -> str:
        return self.__overhangs_3prime.get(part, "")

class BsmBI(RestrictionEnzymeBase):
    
    def __init__(self) -> None:
        super().__init__(RestrictionEnzymes.BSMBI_ANNOTATION_F, 
            RestrictionEnzymes.BSMBI_ANNOTATION_R)


class BbsI(RestrictionEnzymeBase):
    
    def __init__(self) -> None:
         super().__init__(RestrictionEnzymes.BBSI_ANNOTATION_F, 
            RestrictionEnzymes.BBSI_ANNOTATION_R)


class AminoAcids():
    pass


class DNAConstants(Const):
    
    DNA_BASEPAIRS = {
        'A': 'T', 
        'T': 'A', 
        'C': 'G', 
        'G': 'C', 
        'a': 't', 
        't': 'a', 
        'g': 'c', 
        'c': 'g'
    }

    # E. coli codon table: http://www.sci.sdsu.edu/~smaloy/MicrobialGenetics/topics/in-vitro-genetics/codon-usage.html
    # Listed in order of frequency
    RES_TO_CODONS = {
        'A': ['GCG', 'GCC', 'GCA', 'GCT'],
        'C': ['TGC', 'TGT'],
        'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'],
        'F': ['TTT', 'TTC'],
        'G': ['GGC', 'GGT', 'GGG', 'GGA'],
        'H': ['CAT', 'CAC'],
        'I': ['ATT', 'ATC', 'ATA'],
        'K': ['AAA', 'AAG'],
        'L': ['CTG', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA'],
        'M': ['ATG'],
        'N': ['AAC', 'AAT'],
        'P': ['CCG', 'CCA', 'CCT', 'CCC'],
        'Q': ['CAG', 'CAA'],
        'R': ['CGT', 'CGC', 'CGG', 'CGA', 'AGA', 'AGG'],
        'S': ['AGC', 'TCT', 'TCC', 'TCG', 'AGT', 'TCA'],
        'T': ['ACC', 'ACA', 'ACG', 'ACT'],
        'V': ['GTG', 'GTT', 'GTC', 'GTA'],
        'W': ['TGG'],
        'Y': ['TAT', 'TAC'],
        '*': ['TAA', 'TGA', 'TAG'],
    }
    

    def __init__(self) -> None:
        super().__init__()

        # DOES THIS STILL WORK AS INTENDED?
        CODON_TO_RES = self.build_codon_to_res()


    def build_codon_to_res(stlf) -> Dict[str, str]:
        codon_to_res = {}

        for res, codons in DNAConstants.RES_TO_CODONS.items():
            for codon in codons:
                codon_to_res[codon] = res
        
        return codon_to_res