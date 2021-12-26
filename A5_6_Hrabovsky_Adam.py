# link to GitHub repistory: https://github.com/hraboada/Assignment-5-6.git
# for the sequence alignment, I used a tool called MUSCLE
# I provided the .exe file in the submission folder
# in order for the code to work properly, the .exe file for the MUSCLE needs to be in the same repository as the main file (in your local computer)

import Bio
import pandas as pd
from Bio import Entrez, Phylo
from Bio.Seq import UndefinedSequenceError
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

Entrez.email = "be21x008@technikum-wien.at"

# At first, I had to define a new class which will contain all the requires information from database
# I also need to ensure that there won't be the same organism again, but with a different transcript from Genebank
# Hence the equality and hash functions need to be included
class Org:
    def __init__(self, acc_number, organism, sequence, description, percentageGC, instability, aromaticity, isoelectricPoint):
        self.acc_number = acc_number
        self.organism = organism
        self.sequence = sequence
        self.description = description
        self.percentageGC = percentageGC
        self.instability = instability
        self.aromaticity = aromaticity
        self.isoelectricPoint  = isoelectricPoint

    def __eq__(self, other):
        return self.organism == other.organism
    def __hash__(self):
        # hash(custom_object)
        return hash((self.organism))

# In some of the records in Genbank, there is no sequence included
# So I need to make sure, that such records will not be downloaded with this function
# is_defined = False
def defined(sequence):
   is_defined = True
   try:
        bytes(sequence)
   except UndefinedSequenceError:
        is_defined = False
   return is_defined

# Searching for the results from Genbank
# I got this example from the cookbook tutorial
# http://biopython.org/DIST/docs/tutorial/Tutorial.html - my source with
number_of_records = 0
handle = Entrez.egquery(term="ZNF141")
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    if row["DbName"]=="nuccore":
        number_of_records = int(row["Count"])
        print("Number of records in Genbank:",number_of_records)
handle = Entrez.esearch(db="nucleotide", term="ZNF141", retmax = number_of_records, idtype="acc")
record = Entrez.read(handle)
handle.close()

# here I actually download the records I searched for in the previous step
idlist = ",".join(record["IdList"])
handle = Entrez.efetch(db="nucleotide", id=idlist, retmode="xml")
records = Entrez.read(handle)

# at this point, I am getting the required information from previously downloaded records
# and I am loading them into a set
setOrganizmov = set()
i=0
while i<len(records):
    acc_number = records[i]["GBSeq_primary-accession"]
    organism = records[i]["GBSeq_organism"]
    strOrg = str(organism)
    strOrg = strOrg.replace(" ", "_")

    description = records[i]["GBSeq_definition"]
    handle = Entrez.efetch(db="nucleotide", id=acc_number, rettype="gb", retmode="text")
    seq_mon = SeqIO.read(handle, "genbank")
    sequence = seq_mon.seq
    org = Org(acc_number, strOrg, sequence, description, "", "", "", "")

    # calling the function to check if in the record is also included the seuence
    # if not, it won't be added to the set
    jeDef = defined(sequence)
    if jeDef == True:
        setOrganizmov.add(org)
    i=i+1

# here I randomly choose 5 organisms from the set and load them into a list
# I need to ensure that the sequence isn't too long (because in some cases the sequence was just too long for the alignment
# and the the program terminated with an error
# and also I want 5 organisms - not viruses nor synthetic sequences - I am excluding those as well
list5Org = []
for o in setOrganizmov:
    cosi = str(o.organism)
    if len(o.sequence) > 30000:
        continue
    elif cosi.__contains__("virus") or cosi.__contains__("synthetic"):
        continue
    else:
        list5Org.append(o)
    if len(list5Org) == 5:
        break

# here I create FASTA files (file name is the organism)
# files contain the name of the organism, its accession number and its sequence
# and I load them into a list
list_of_files = []
for o in list5Org:
    o.percentageGC = GC(o.sequence)
    name = str(o.organism)
    name_of_the_file = name+".fasta"
    list_of_files.append(name_of_the_file)
    idNumber = o.organism + "-" + o.acc_number
    seq = SeqRecord(seq=o.sequence, id=idNumber, name=o.organism, description='', dbxrefs=[] )
    # writting into the files
    SeqIO.write(seq,name_of_the_file,"fasta")

# protein analysis
# I am loading the sequences from the previously created files,
# translationg them to proteins
# and analyzing them - calculating the required info:
# aromaticity, insatbility, isoelectric point and secondary structure
# https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html - my source for protein analysis
print("\nProtein structure analysis:")
for file_name in list_of_files:
    org = None
    for o in list5Org:
        subor = file_name.replace(".fasta", "")
        if subor == o.organism:
            org = o
    nucleotide = SeqIO.read(file_name, "fasta")

    # finding number of letters to trim
    trim_char = len(nucleotide) % 3

    # trim sequence to translate it
    if trim_char > 0:
        trim_seq = nucleotide[:-trim_char]
    else:
        trim_seq = nucleotide

    protein = trim_seq.translate()
    sequence_of_protein = protein.seq
    sequence_of_protein = str(sequence_of_protein)
    protein_analysis = ProteinAnalysis(sequence_of_protein)
    org.aromaticity = protein_analysis.aromaticity()
    org.isoelectricPoint = protein_analysis.isoelectric_point()
    if sequence_of_protein.__contains__("*"):
        sequence_of_protein = sequence_of_protein.replace("*", "")
    if sequence_of_protein.__contains__("X"):
        sequence_of_protein = sequence_of_protein.replace("X", "")
    protein_analysis = ProteinAnalysis(sequence_of_protein)
    org.instability = protein_analysis.instability_index()
    sequence_of_protein = protein.seq
    sequence_of_protein = str(sequence_of_protein)
    protein_analysis = ProteinAnalysis(sequence_of_protein)
    sec_struc = protein_analysis.secondary_structure_fraction()
    print(org.organism)
    print("Fractions of amino acids that tend to be helix:", sec_struc[0])
    print("Fractions of amino acids that tend to be turn:", sec_struc[1])
    print("Fractions of amino acids that tend to be sheet:", sec_struc[2])

# performing the sequence alignment and creating a phylogenetic tree from it
# https://www.youtube.com/watch?v=wBdz3vFQ4Ks&t=922s
# I was following this tutorial on youtube and for the alignment they use a tool called MUSCLE
# I included an .exe file in my submission folder
# it needs to be in the same repository as the main file for it to work
# https://biopython.org/docs/1.75/api/Bio.Align.Applications.html

input_sequences = "input_sequences.fasta"
output_alignment = "output_alignment.fasta"

# creating one file containing all the sequences
listSeqIO = []
for name in list_of_files:
    listSeqIO.append(SeqIO.read(name, "fasta"))
SeqIO.write(listSeqIO, input_sequences, "fasta")

# so the input file for the MUSCLE tool is the file containing all of the sequences
# and it returns another file with the alignments performed

muscle_exe = r"muscle3.8.31_i86win32.exe"
def sequence_alignment (fasta):
    muscle_cline = MuscleCommandline(muscle_exe, input=fasta, out=output_alignment)
    muscle_cline()
sequence_alignment(input_sequences) # performing the sequence alignment

# opening and reading the alignment file
with open(output_alignment, "r") as aln:
    alignment = AlignIO.read(aln, "fasta")
# opening and initiating the distance calculator using the identity model
from Bio.Phylo.TreeConstruction import DistanceCalculator
calculator = DistanceCalculator("identity")
# performing the distance calculation
distance_matrix = calculator.get_distance(alignment)
# opening and initiating the tree consturctor
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
constructor = DistanceTreeConstructor(calculator)
# building the actual tree
organism_tree = constructor.build_tree(alignment)
organism_tree.rooted = True

# here I created lists and I loaded them with all the required information
# so that I could build a dataframe table to properly show my results
organism = []
accNumber = []
title = []
dlzka = []
percentage = []
instability = []
aroma = []
isoElPoint = []

for org in list5Org:
    organism.append(org.organism)
    accNumber.append(org.acc_number)
    title.append(org.description)
    length = len(org.sequence)
    length = int(length)
    dlzka.append(length)
    percentage.append(org.percentageGC)
    instability.append(org.instability)
    aroma.append(org.aromaticity)
    isoElPoint.append(org.isoelectricPoint)

result = pd.DataFrame(list(zip(organism, accNumber, title, dlzka, percentage, instability, aroma, isoElPoint)),
          columns = ["Organism", "Accession number", "Title", "Length of the sequence", "GC percentage", "Instability", "Aromaticity", "Isoelectric point"])
pd.set_option("display.max_rows", None, "display.max_columns", None)
print("\nTable with all results:")
print(result)

print("\nDistance matrix for the phylogenetic tree:")
print(distance_matrix)
print("\nPhylogenetic tree:")
Phylo.draw_ascii(organism_tree)