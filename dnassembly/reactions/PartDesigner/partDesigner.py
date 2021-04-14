"""
GGpart.py

Modified from code originally written by Will DeLoache.
"""

import sys, os
import pdb
from Bio import SearchIO, SeqIO
#from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Restriction import BbsI, BsaI, Restriction, FormattedSeq

from .GGfrag import *
from .MergeSegments import *
from .OHfinder import *
from .AssemblyInstructions import *
from .removeRS import *
from .findTemplates import *

#from CodonOptimize.codonOptimize import codonOptimize
#from makeBLASTdb import makeBLASTdb

# BASE_PATH is the absolute path of ../../.. relative to this script location
#BASE_PATH = reduce (lambda l,r: l + os.path.sep + r, os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep )[:-3] )
#sys.path.append( os.path.join( BASE_PATH, "main" ) )
#sys.path.append('../../../main')
#from siteinfo import *

#Define static variables
gBlockMaxSize = 3000
PCAMaxSize = 800
gBlockMinSize = 125
oligoAssemblySize = 100
annealingLength = 20
oligoTM = 48

#Overhangs of the entry vector
entry_ohs = ('TTCT', 'AACG')

#Sequence between entry vector OH and part ohs
entry_extensions = ('GAAGACTC', 'GAGT')

#part assembly overhangs
partOH_5 = {'1': 'TCAA',
'2a': 'ACGA',
'2b': 'AACG',
'3a': 'CACC',
'3b': 'TTCT',
'3c': 'AGCA',
'3d': 'AGGC',
'3e': 'TCCA',
'4a': 'ATCC',
'4b': 'GTAA',
'5': 'GCTA',
'6': 'ACTC',
'7': 'ATTG'
}

partOH_3 = {'1': 'ACGA',
'2a': 'AACG',
'2b': 'CACC',
'3a': 'TTCT',
'3b': 'AGCA',
'3c': 'AGGC',
'3d': 'TCCA',
'3e': 'ATCC',
'4a': 'GTAA',
'4b': 'GCTA',
'5': 'ACTC',
'6': 'ATTG',
'7': 'TCAA'
}

#Sequence between part and part overhangs
partextension_5 = {'1': 'GAATTCGCATCTAGA',
'2a': '',
'2b': '', # This is where the standardized 5' UTR and 3' UTR can go
'3a': 'GC',
'3b': 'GG',
'3c': '',
'3d': 'TC',
'3e': 'TCCA',
'4a': 'GG', # This is where the standardized 5' UTR and 3' UTR can go
'4b': 'TA',
'5': '',
'6': '',
'7': ''
}

partextension_3 = {'1': '',
'2a': '',
'2b': '',
'3a': '',
'3b': '',
'3c': 'GT',
'3d': '',
'3e': 'GT',
'4a': '',
'4b': 'TAA',
'5': 'ACTAGTGCACTGCAG',
'6': '',
'7': ''
}

multigene_ohs = ["CTGA", "CCAA", "GATG", "GTTC", "GGTA", "AAGT", "AGCA"] # - NEED TO FIX

#Sequence to add to ends of a PCR in order to digest with BsmBI for entry into Stage 1 vector
tails = {"BsmBI":("GCATCGTCTCA", "TGAGACGGCAT"),
        "BbsI":("GCATGGTCTCA", "TGAGACCGCAT")}

entryVectors = {"OPL4772":"gcatcaaatgaaactgcaatttattcatatcaggattatcaataccatatttttgaaaaagccgtttctgtaatgaaggagaaaactcaccgaggcagttccataggatggcaagatcctggtatcggtctgcgattccgactcgtccaacatcaatacaacctattaatttcccctcgtcaaaaataaggttatcaagtgagaaatcaccatgagtgacgactgaatccggtgagaatggcaaaagtttatgcatttctttccagacttgttcaacaggccagccattacgctcgtcatcaaaatcactcgcatcaaccaaaccgttattcattcgtgattgcgcctgagccagacgaaatacgcgatcgctgttaaaaggacaattacaaacaggaatcgaatgcaaccggcgcaggaacactgccagcgcatcaacaatattttcacctgaatcaggatattcttctaatacctggaatgctgtttttccggggatcgcagtggtgagtaaccatgcatcatcaggagtacggataaaatgcttgatggtcggaagaggcataaattccgtcagccagtttagtctgaccatctcatctgtaacatcattggcaacgctacctttgccatgtttcagaaacaactctggcgcatcgggcttcccatacaagcgatagattgtcgcacctgattgcccgacattatcgcgagcccatttatacccatataaatcagcatccatgttggaatttaatcgcggcctcgacgtttcccgttgaatatggctcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccagatacctgaaacaaaacccatcgtacggccaaggaagtctccaataactgtgatccaccacaagcgccagggttttcccagtcacgacgttgtaaaacgacggccagtcatgcataatccgcacgcatctggaataaggaagtgccattccgcctgacctttctagagacgttctagagcacagctaacaccacgtcgtccctatctgctgccctaggtctatgagtggttgctggataactttacgggcatgcataaggctcgtataatatattcagggagaccacaacggtttccctctacaaataattttgtttaacttttactagatcacacaggaaagtactagatgcgtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgatgttaatgggcacaaattttctgtcagtggagagggtgaaggtgatgcaacatacggaaaacttacccttaaatttatttgcactactggaaaactacctgttccatggccaacacttgtcactactttcggttatggtgttcaatgctttgcgagatacccagatcatatgaaacagcatgactttttcaagagtgccatgcccgaaggttatgtacaggaaagaactatatttttcaaagatgacgggaactacaagacacgtgctgaagtcaagtttgaaggtgatacccttgttaatagaatcgagttaaaaggtattgattttaaagaagatggaaacattcttggacacaaattggaatacaactataactcacacaatgtatacatcatggcagacaaacaaaagaatggaatcaaagttaacttcaaaattagacacaacattgaagatggaagcgttcaactagcagaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattacctgtccacacaatctgccctttcgaaagatcccaacgaaaagagagaccacatggtccttcttgagtttgtaacagctgctgggattacacatggcatggatgaactatacaaatagtagtactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttatacgtctcgaacgaggctaggtggaggctcagtgatgataagtctgcgatggtggatgcatgtgtcatggtcatagctgtttcctgtgtgaaattgttatccgctcagagggcacaatcctattccgcgctatccgacaatctccaagacattaggtggagttcagttcggcgtatggcatatgtcgctggaaagaacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctctattcaacaaagccgccgtcccgtcaagtcagcgtaaatgggtagggggcttcaaatcgtcctcgtgataccaattcggagcctgcttttttgtacaaacttgttgataatggcaattcaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttagaaaaactcatcga"}

entryVectorDigests = {"OPL4772":"aggctaggtggaggctcagtgatgataagtctgcgatggtggatgcatgtgtcatggtcatagctgtttcctgtgtgaaattgttatccgctcagagggcacaatcctattccgcgctatccgacaatctccaagacattaggtggagttcagttcggcgtatggcatatgtcgctggaaagaacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctctattcaacaaagccgccgtcccgtcaagtcagcgtaaatgggtagggggcttcaaatcgtcctcgtgataccaattcggagcctgcttttttgtacaaacttgttgataatggcaattcaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttagaaaaactcatcgagcatcaaatgaaactgcaatttattcatatcaggattatcaataccatatttttgaaaaagccgtttctgtaatgaaggagaaaactcaccgaggcagttccataggatggcaagatcctggtatcggtctgcgattccgactcgtccaacatcaatacaacctattaatttcccctcgtcaaaaataaggttatcaagtgagaaatcaccatgagtgacgactgaatccggtgagaatggcaaaagtttatgcatttctttccagacttgttcaacaggccagccattacgctcgtcatcaaaatcactcgcatcaaccaaaccgttattcattcgtgattgcgcctgagccagacgaaatacgcgatcgctgttaaaaggacaattacaaacaggaatcgaatgcaaccggcgcaggaacactgccagcgcatcaacaatattttcacctgaatcaggatattcttctaatacctggaatgctgtttttccggggatcgcagtggtgagtaaccatgcatcatcaggagtacggataaaatgcttgatggtcggaagaggcataaattccgtcagccagtttagtctgaccatctcatctgtaacatcattggcaacgctacctttgccatgtttcagaaacaactctggcgcatcgggcttcccatacaagcgatagattgtcgcacctgattgcccgacattatcgcgagcccatttatacccatataaatcagcatccatgttggaatttaatcgcggcctcgacgtttcccgttgaatatggctcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccagatacctgaaacaaaacccatcgtacggccaaggaagtctccaataactgtgatccaccacaagcgccagggttttcccagtcacgacgttgtaaaacgacggccagtcatgcataatccgcacgcatctggaataaggaagtgccattccgcctgacct"}

class GGpart():
    def __init__(self, partName, leftPartType, rightPartType, seq, method="None", plasmidName="", inputSeqType="DNA", rsToRemove=['BbsI','BsmBI'], removeHomology=False,
                skipPartPlasmid=False, fiveprime="", threeprime="", maxPrimerLength=60, possibleTemplates={}, databaseTemplates=[], addOHs=[]):
        self.partName = partName
        self.leftPartType = leftPartType
        self.rightPartType = rightPartType
        self.inputSeq = seq
        self.method = method
        self.plasmidName = plasmidName

        self.inputSeqType = inputSeqType
        self.rsToRemove = rsToRemove
        self.removeHomology = removeHomology

        self.maxPrimerLength = maxPrimerLength
        self.skipPartPlasmid = skipPartPlasmid

        self.possibleTemplates = possibleTemplates
        for each in self.possibleTemplates:
            self.possibleTemplates[each]["seq"] = self.possibleTemplates[each]["seq"].replace("\r","").replace("\n","").replace(" ","")

        self.additionalOHs = addOHs

        self.databaseTemplates = databaseTemplates
        self.warnings = []
        self.errors = []

        ## These will all be initialized later unless the an error is thrown in which case they will be left blank
        self.vectorName = ""
        self.vectorSeq = ""
        self.vectorDigest = ""
        self.maxOverhangConstraints = 1
        self.plasmidSeq = ""
        self.enzyme = ""
        self.subclone = ""
        self.partSeq = ""
        self.GGfrags = []
        self.fragments = []

        if self.checkInput_beforeOptimization():
            #Generate partSeq by removing RS and codon optimizing
            self.partSeq = ""
            self.partSeq = self.partSeq.lower() #change all bases to lowercase
            self.removeInternalRS()

            if self.checkInput_afterOptimization():
                #Initialize GG fragment by making it a single piece
                self._initGGfrag(fiveprime, threeprime)

                #Find the best overhangs
                self.pickGGfrags()

                #Combine adjacent Oligo Assemblies if possible
                if self.method in ["None", "Oligo Assembly"]:
                    self.mergeOligoAssemblies()
                try:
                    self.vectorDigest = entryVectorDigests[self.vectorName]
                    self.vectorSeq = entryVectors[self.vectorName]
                except:
                    self.vectorDigest = ""
                    self.vectorSeq = ""
                self.plasmidSeq = self.getPlasmidSeq()

                #List of AssemblyInstruction objects
                self.fragments = self.getAssemblyInstructions()


#       ################ Other functions #################
    # Check to make sure input is in the form of DNA letters
    def _verify_alphabet(sequence):
        letters = sequence.alphabet.letters
        if not letters:
            raise ValueError("Alphabet does not define letters.")
        for letter in sequence:
            if letter not in letters:
                return False
        return True

    # Check to ensure the input sequence is DNA and does not break part 1 and 5 construction
    def checkInput_beforeOptimization(self):
        #seq = Seq(self.inputSeq.upper())
        #if not (_verify_alphabet(seq)):
        #    self.errors.append("The DNA sequence you entered in invalid.")
        if self.leftPartType in ["1", "5"] and "BsmBI" in self.rsToRemove:
            self.warnings.append("You chose to remove BsmBI sites from a type connector part. Doing so will break multigene assembly.")
        if self.rightPartType in ["1", "5"] and "BsmBI" in self.rsToRemove:
            self.warnings.append("You chose to remove BsmBI sites from a type connector part. Doing so will break multigene assembly.")
        if len(self.errors) != 0:
            return False
        else:
            return True

    # Use custom script to remove defined RS. By default these are defined as BbsI and BsmBI
    def removeInternalRS(self): #goto removeRS.py
        self.partSeq = removeRS(self.inputSeq, self.rsToRemove)

    def checkInput_afterOptimization(self):
        f_seq = FormattedSeq(Seq(self.partSeq), True)

        #Check to make sure BbsI/BsmBI sites aren't present
        if BbsI.search(f_seq) and BsmBI.search(f_seq):
            #comment the line below to let through parts with BbsI and BsmBI
            #be careful, though, these assemblies may be problematic
            #self.errors.append("Your part contains both BbsI and BsmBI sites and cannot be assembled using golden gate.")
            pass
        elif self.leftPartType in ['1', '2a', '2b', '3a', '3b', '3c', '3d', '3e', '4a', '4b', '5', '6', '7'] and BbsI.search(f_seq):
            self.errors.append("Your part contains a BbsI site which must be removed prior to assembly.")

        elif self.rightPartType in ['1', '2a', '2b', '3a', '3b', '3c', '3d', '3e', '4a', '4b', '5', '6', '7'] and BbsI.search(f_seq):
            self.errors.append("Your part contains a BbsI site which must be removed prior to assembly.")

        elif self.leftPartType in ['2a', '2b', '3a', '3b', '3c', '3d', '3e', '4a', '4b', '6', '7'] and BsmBI.search(f_seq):
            self.errors.append("Your part contains a BsmBI site which must be removed prior to assembly.")

        elif self.rightPartType in ['2a', '2b', '3a', '3b', '3c', '3d', '3e', '4a', '4b', '6', '7'] and BsmBI.search(f_seq):
            self.errors.append("Your part contains a BsmBI site which must be removed prior to assembly.")

        #Check to make sure connector parts have BsmBI sites
        bsmBIFor = self.partSeq.upper().count(BsmBI.site)
        bsmBIRev = self.partSeq.upper().count(str(Seq(BsmBI.site).reverse_complement()))
        if str(self.leftPartType) == "1" or str(self.rightPartType) == "1":
            if bsmBIFor + bsmBIRev != 1:
                self.warnings.append("Your type 1 part should have exactly 1 BsmBI site. Consider modifying for multigene assembly.")
        if str(self.leftPartType) == "5" or str(self.rightPartType) == "5":
            if bsmBIFor + bsmBIRev != 1:
                self.warnings.append("Your type 5 part should have exactly 1 BsmBI site. Consider modifying for multigene assembly.")

        #Warn if ORF has a start codon
        if self.leftPartType in ["3a"]:
            if str(self.partSeq.upper())[:3] != "ATG":
                self.warnings.append("Your part is missing a start codon.")
        #Warn if ORF has a stop codon
        if str(self.rightPartType) in ['3a', '3b', '3c', '3d', '3e', '4a']:
            if len(self.partSeq) % 3 != 0:
                self.warnings.append("Your part appears to be out of frame (length is not a multiple of 3). If this is a coding sequence, check to make sure it is correct.")
            if Seq(self.partSeq).translate().find("*") > -1:
                self.warnings.append("Your part has a stop codon. If it is not removed, the part cannot be used for making N-terminal fusions.")

        if len(self.errors) != 0:
            return False
        else:
            return True


    def _initGGfrag(self, fiveprime, threeprime):
        #fiveprime and threeprime are variables that will get passed in the case of custom part design
        #only set to true if you have to gel purify the fragments
        self.subclone = False

        frag_seq = self.partSeq
        #These are the strings that will get appended onto the ends of the fratment
        fiveprimeOH = ""
        threeprimeOH = ""
        fiveprimeExt = ""
        threeprimeExt = ""

        #Fix this section for custom part design
        if self.leftPartType == "Custom" or self.rightPartType == "Custom":
            part_oh = (fiveprime[:4].upper(), threeprime[-4:].upper())
            part_extension = (fiveprime[4:], threeprime[:-4])
            #Custom parts will go into the part entry vector unless skipPartPlasmid is selected
            entry_oh = entry_ohs
            entry_extension = entry_extensions
            #could possibly change this default if I want to allow custom mV parts to be made
            self.enzyme = "BsmBI"
        else:
            self.enzyme = "BsmBI"
            fiveprimeOH = entry_ohs[0].upper()
            threeprimeOH = entry_ohs[1].upper()
            fiveprimeExt = entry_extensions[0] + partOH_5[self.leftPartType] + partextension_5[self.leftPartType]
            threeprimeExt = partextension_3[self.rightPartType] + partOH_3[self.rightPartType] + entry_extensions[1]

            #Set vector
            self.subclone = False
            self.vectorName = "OPL4772" #Remember to hardcode in this vector

            newFrag = GGfrag(fiveprimeOH, fiveprimeExt, frag_seq, threeprimeExt, threeprimeOH)
            self.GGfrags = [newFrag]

    def pickGGfrags(self):
        allowableSize = 0
        forced_method = ""
        if self.method == "gBlocks":
            allowableSize = gBlockMaxSize - 24
        elif self.method == "Oligo Assembly":
            allowableSize = gBlockMinSize - 3
        elif self.method == "PCA":
            allowableSize = PCAMaxSize - 24

        #Divide the initialized GGfrag into chunks based on assembly method
        if self.method in ["gBlocks", "Oligo Assembly", "PCA"]:
            self.GGfrags[0].forced_method = self.method
            self.GGfrags = self.GGfrags[0].divideBySize(allowableSize)
        else:
            ##### Run findTemplates to find PCR templates in database####################
            '''Skip this functionality for now - in the future we may want to add the ability to use the API to search Benchling
            dbpath = makeBLASTdb(self.possibleTemplates, self.databaseTemplates)

            results = findTemplates(self.GGfrags[0].seq, dbpath) #results is in the form of 0=
            self.GGfrags[0].seq = results[0] #duplicates sequence? Not sure what this is for
            con = MySQLdb.connect(DATABASE_HOST, DATABASE_USER, DATABASE_PASSWD, DATABASE_NAME);
            cur = con.cursor(MySQLdb.cursors.DictCursor)
            for templateID in results[1].keys():
                try:
                    cur.execute("SELECT * FROM Plasmids WHERE (ID=%s)" % templateID)
                    plasmid = cur.fetchone()
                    seq = plasmid["Sequence"]
                    if plasmid["Linear"] == 1:
                        self.possibleTemplates[plasmid["Name"]] = {"seq":plasmid["Sequence"], "Linear":True}
                    else:
                        self.possibleTemplates[plasmid["Name"]] = {"seq":plasmid["Sequence"], "Linear":False}
                except:
                    if templateID in ["S288C_Genomic_DNA", "DH10B_Genomic_DNA", "GS115_Genomic_DNA"]:
                        self.possibleTemplates[templateID] = {"seq":results[1][templateID], "Linear":True}
            '''
            ################### End Find Templates from Database ######################
            self.GGfrags = self.GGfrags[0].divideByCaseChange() #seems like the uppercase acts as a divider? I think the only uppercase should be from the part entry OH
            self._mergeFragments() #
            ### Split sequence without template into gBlock sized chunks###
            if self.method == "None":
                tempFrags = []
                for index, each in enumerate(self.GGfrags):
                    #if not findTemplate(each.seq, self.possibleTemplates):
                    if each.getLength() > 3100:
                        each.forced_method = "PCA"
                        tempFrags += each.divideBySize(PCAMaxSize - 24)
                    elif each.getLength() < gBlockMinSize:
                        each.forced_method = "Oligo Assembly"
                        tempFrags += [each]
                    #If length is just over gBlock length, make the assembly include one oligo assembly to save on a gBlock
                    elif 0 < each.getLength() % (gBlockMaxSize - 24) < gBlockMinSize - 4 * each.getLength()/(gBlockMaxSize-24):
                        oligoAssemSize = max(20, each.getLength() % (gBlockMaxSize - 24) + 4 * each.getLength()/(gBlockMaxSize-24))
                        if index == 0:
                            oligoFrag = GGfrag(each.getFragSeq()[:4], "", each.getFragSeq()[4:oligoAssemSize], "", "", forced_method="Oligo Assembly")
                        else:
                            oligoFrag = GGfrag("", "", each.getFragSeq()[:oligoAssemSize], "", "", forced_method="Oligo Assembly")
                        if index == len(self.GGfrags)-1:
                            gBlockFrag = GGfrag("", "", each.getFragSeq()[oligoAssemSize:-4], "", each.getFragSeq()[-4:], forced_method="gBlocks")
                        else:
                            gBlockFrag = GGfrag("", "", each.getFragSeq()[oligoAssemSize:], "", "", forced_method="gBlocks")
                        tempFrags += [oligoFrag]
                        tempFrags += gBlockFrag.divideBySize(gBlockMaxSize - 24)
                    else:
                        each.forced_method = "gBlocks"
                        tempFrags += each.divideBySize(gBlockMaxSize - 24)
                    #else:
                        #tempFrags.append(each)
                self.GGfrags = tempFrags
        self.pickGGfragOHs()
        if self.maxOverhangConstraints > 12:
            self.warnings.append("Due to its complexity, this assembly uses overhangs that are somewhat similar and may be prone to misannealing.")
        for each in self.GGfrags:
            each.assignTails(tails[self.enzyme], self.enzyme)

    #combines items from self.GGfrags where possible based on maxPrimer length
    #works by converting the GGfrags into segments of form [oh+ext, seq, ext+oh]
    #calls mergeSegments to combine them and then turns them back into GGfrags.
    #If I have time, this should be written since it is stupid to convert into segments
    def _mergeFragments(self):
        #Maximum length of primer
        maxPrimerLength = self.maxPrimerLength - 11 #11 = length of type IIs tail
        if len(self.GGfrags) <= 1: #no need to merge
            return
        #a segment is just a ist of 3 sequences, the middle piece is the template and the flanking pieces are extra sequence
        segments = []
        for i in range(len(self.GGfrags)):
            newSeg = [self.GGfrags[i].fiveprimeOH + self.GGfrags[i].fiveprimeExt, self.GGfrags[i].seq, self.GGfrags[i].threeprimeExt + self.GGfrags[i].threeprimeOH]
            segments.append(newSeg)
        mergedSegs = mergeSegments(segments, maxPrimerLength, annealingLength)
        #convert segments into GGfrags
        newFrags = []
        for i in range(len(mergedSegs)):
            fiveprimeOH = ""
            fiveprimeExt = mergedSegs[i][0]
            seq = mergedSegs[i][1]
            threeprimeExt = mergedSegs[i][2]
            threeprimeOH = ""
            if i == 0:
                fiveprimeOH = self.GGfrags[i].fiveprimeOH
                fiveprimeExt = mergedSegs[i][0][4:]
            if i == len(mergedSegs) - 1:
                threeprimeOH = self.GGfrags[-1].threeprimeOH
                threeprimeExt = mergedSegs[i][2][:-4]
            newFrag = GGfrag(fiveprimeOH, fiveprimeExt, seq, threeprimeExt, threeprimeOH)
            newFrags.append(newFrag)
        self.GGfrags = newFrags

    #Find the best overhangs
    def pickGGfragOHs(self):
        allowableSize = 0
        if self.method == "gBlocks":
            allowableSize = gBlockMaxSize - 22
        elif self.method == "Oligo Assembly":
            allowableSize = gBlockMinSize

        #Find possible overhangs at each junction
        possOHs = []
        if self.method in ["gBlocks", "Oligo Assembly"]:
            possOHs = findPossOH_byFragLength(self.GGfrags, allowableSize, annealingLength)
        else:
            possOHs = findPossOH_byPrimerLength(self.GGfrags, self.maxPrimerLength - 11, annealingLength, gBlockMaxSize, self.enzyme)
        if possOHs == False:
            self.errors.append("No valid overhangs could be found. Consider increasing the maximum primer length or specifying a different assembly strategy.")
            return False

        #Compile all additional overhangs to avoid conflicts with
        #AHN: I think this is only necessary if you are trying to skip the part assembly step? May comment out
        additionalOHs = list(self.additionalOHs)

        '''
        if self.skipPartPlasmid and self.GGtype != "Circular":
            if self.GGtype in ["mV", "mVa", "mVb", "V", "Va", "Vb"]:
                additionalOHs += list(multigene_ohs)
            else:
                for each in part_ohs:
                    if part_ohs[each][0] not in additionalOHs:
                        additionalOHs.append(part_ohs[each][0])
            if self.GGfrags[0].fiveprimeOH not in additionalOHs:
                additionalOHs.append(self.GGfrags[0].fiveprimeOH)
            if self.GGfrags[-1].threeprimeOH not in additionalOHs:
                additionalOHs.append(self.GGfrags[-1].threeprimeOH)
        '''

        #Add existing overhangs at beginning and end to the list
        existingOHs = []
        existingOHs.append((self.GGfrags[0].fiveprimeOH, 0))
        existingOHs.append((self.GGfrags[-1].threeprimeOH, 0))

        #Find the best set of overhangs
        constraints = 1
        bestOHs = False
        while not bestOHs:
            if constraints > 16:
                self.errors.append("No valid overhangs could be found")
                return False
            bestOHs = findBestOverhangs(existingOHs, possOHs, additionalOHs, constraints)
            constraints += 1
        #reorder the overhangs so that the 3' OH is last
        bestOHs.append(bestOHs.pop(1))
        if constraints - 1 > self.maxOverhangConstraints:
            self.maxOverhangConstraints = constraints - 1
        self.GGfrags = reindexGGfrags(self.GGfrags, bestOHs)

    def mergeOligoAssemblies(self):
        fragsMerged = True
        while fragsMerged:
            fragsMerged = False
            for i in range(len(self.GGfrags)-1):
                currFrag = self.GGfrags[i]
                nextFrag = self.GGfrags[i+1]
                if (currFrag.getLength() + nextFrag.getLength()) < oligoAssemblySize:
                    fiveprimeOH = currFrag.fiveprimeOH
                    fiveprimeExt = currFrag.fiveprimeExt + currFrag.seq + currFrag.threeprimeExt + currFrag.threeprimeOH + nextFrag.fiveprimeExt
                    seq = nextFrag.seq
                    threeprimeExt = nextFrag.threeprimeExt
                    threeprimeOH = nextFrag.threeprimeOH
                    mergedFrag = GGfrag(fiveprimeOH, fiveprimeExt, seq, threeprimeExt, threeprimeOH)
                    self.GGfrags.pop(i)
                    self.GGfrags[i] = mergedFrag
                    fragsMerged = True
                    break

    def getPlasmidSeq(self):
        s = ""
        for index, each in enumerate(self.GGfrags):
            s += each.fiveprimeOH + each.fiveprimeExt + each.seq + each.threeprimeExt
            if index == len(self.GGfrags) - 1:
                s += each.threeprimeOH
        if self.vectorName == "None":
            #temporary solution until we support skipPartPlasmid better
            if self.skipPartPlasmid:
                s = tails[self.enzyme][0] + s + tails[self.enzyme][1]
        elif self.vectorName in entryVectorDigests.keys():
            s += entryVectorDigests[self.vectorName]
        else:
            vector = "OPL4772"

            ''' #Just specify the one vector
            if self.GGtype in ["1","234","2","2a","2b","3","3a","3ab","3ac","3ad","3ae","3b","3bc","3bd","3be","3c","3cd","3ce","3d","3de","3e","4","4a","4b","5","6","7"]:
                vector = "OPL4772"
            '''

            s = entry_ohs[0] + entry_extensions + s
            s += entry_extensions[1] + entry_ohs[1]
            s += entryVectorDigests[vector]
            self.vectorDigest = self.GGfrags[-1].threeprimeOH + entry_extensions[1] + entry_ohs[1] + entryVectorDigests[vector] + entry_ohs[0] + entry_extensions[0] + self.GGfrags[0].fiveprimeOH
            self.vectorSeq = "NNNNNNNNNN" + self.vectorDigest #mike changed to self.vectorDigest
        return s

    def getAssemblyInstructions(self):
        assemIns = []
        results = self.getPrimersAndMethods()
        primers = results[0]
        methods = results[1]

##### Add unique cutters surrounding each PCA fragment ###########
        if self.method == "PCA":
            tails[self.enzyme] = (tails[self.enzyme][0][:4] + "GAATTCA" + tails[self.enzyme][0][4:], tails[self.enzyme][1][:-4] + "AGGATCC" + tails[self.enzyme][1][-4:])
##################################################################

        for i in range(len(self.GGfrags)):
            assemIns.append(AssemblyInstructions(methods[i], primers[i], self.GGfrags[i], self.possibleTemplates))
        return assemIns

    def getPrimersAndMethods(self):
        primers = []
        methods = []
        for each in self.GGfrags:
            if self.method == "gBlocks":
                primers.append([])
                methods.append("gBlocks")
            elif self.method == "Oligo Assembly":
                primers.append(each.getOligoAssemPrimers(self.maxPrimerLength))
                methods.append("Oligo Assembly")
            #elif self.method == "PCA":
            #	primers.append(each.getPCAprimers(self.maxPrimerLength))
            #	methods.append("PCA")
            #Check to make sure PCR templates are present
            elif self.method == "PCR":
                primers.append(each.getPCRprimers(self.maxPrimerLength, annealingLength, oligoTM=oligoTM))
                methods.append("PCR")
                if not findTemplate(each.seq, self.possibleTemplates):
                    message = "PCR templates could not be found for this assembly. If you have a template, add it to the templates text box. If you don't have a template, consider using an assembly method other than PCR."
                    if message not in self.errors:
                        self.errors.append(message)
            elif self.method == "None":
                if not self.possibleTemplates: #Change - if no template, then
                    if each.getLength() <= gBlockMinSize:
                        primers.append(each.getOligoAssemPrimers(self.maxPrimerLength))
                        methods.append("Oligo Assembly")
                    #elif each.getLength() + len(tails[self.enzyme][0]) + len(tails[self.enzyme][1]) >  gBlockMaxSize:
                        #primers.append(each.getPCAprimers(self.maxPrimerLength))
                        #methods.append("PCA")
                    else:
                        primers.append([])
                        methods.append("gBlocks")
                else:
                    if each.getLength() <= oligoAssemblySize:
                        primers.append(each.getOligoAssemPrimers(self.maxPrimerLength))
                        methods.append("Oligo Assembly")
                    else:
                        primers.append(each.getPCRprimers(self.maxPrimerLength, annealingLength, oligoTM=oligoTM))
                        methods.append("PCR")
        return (primers, methods)
