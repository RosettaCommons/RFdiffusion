import torch
import numpy as np

num2aa=[
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL',
    'UNK','MAS',
    ]

# Mapping 3 letter AA to 1 letter AA (e.g. ALA to A)
one_letter = ["A", "R", "N", "D", "C", \
             "Q", "E", "G", "H", "I", \
             "L", "K", "M", "F", "P", \
             "S", "T", "W", "Y", "V", "?", "-"]

aa2num= {x:i for i,x in enumerate(num2aa)}

aa_321 = {a:b for a,b in zip(num2aa,one_letter)}
aa_123 = {val:key for key,val in aa_321.items()}


# create single letter code string from parsed integer sequence
def seq2chars(seq):
    out = ''.join([aa_321[num2aa[a]] for a in seq])
    return out

# full sc atom representation (Nx14)
aa2long=[
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # ala
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2",  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD "," HE ","1HH1","2HH1","1HH2","2HH2"), # arg
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HD2","2HD2",  None,  None,  None,  None,  None,  None,  None), # asn
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None,  None), # asp
    (" N  "," CA "," C  "," O  "," CB "," SG ",  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ",  None,  None,  None,  None,  None,  None,  None,  None), # cys
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE2","2HE2",  None,  None,  None,  None,  None), # gln
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ",  None,  None,  None,  None,  None,  None,  None), # glu
    (" N  "," CA "," C  "," O  ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  ","1HA ","2HA ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None), # gly
    (" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2",  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HD2"," HE1"," HE2",  None,  None,  None,  None,  None,  None), # his
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1",  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG2","2HG2","3HG2","1HG1","2HG1","1HD1","2HD1","3HD1",  None,  None), # ile
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ","1HD1","2HD1","3HD1","1HD2","2HD2","3HD2",  None,  None), # leu
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ ",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ","1HE ","2HE ","1HZ ","2HZ ","3HZ "), # lys
    (" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE ",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE ","2HE ","3HE ",  None,  None,  None,  None), # met
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ ",  None,  None,  None," H  "," HA ","1HB ","2HB "," HD1"," HD2"," HE1"," HE2"," HZ ",  None,  None,  None,  None), # phe
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD ",  None,  None,  None,  None,  None,  None,  None," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ",  None,  None,  None,  None,  None,  None), # pro
    (" N  "," CA "," C  "," O  "," CB "," OG ",  None,  None,  None,  None,  None,  None,  None,  None," H  "," HG "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None), # ser
    (" N  "," CA "," C  "," O  "," CB "," OG1"," CG2",  None,  None,  None,  None,  None,  None,  None," H  "," HG1"," HA "," HB ","1HG2","2HG2","3HG2",  None,  None,  None,  None,  None,  None), # thr
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2"," H  "," HA ","1HB ","2HB "," HD1"," HE1"," HZ2"," HH2"," HZ3"," HE3",  None,  None,  None), # trp
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH ",  None,  None," H  "," HA ","1HB ","2HB "," HD1"," HE1"," HE2"," HD2"," HH ",  None,  None,  None,  None), # tyr
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2",  None,  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG1","2HG1","3HG1","1HG2","2HG2","3HG2",  None,  None,  None,  None), # val
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # unk
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # mask
]

# build the "alternate" sc mapping
aa2longalt=[
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # ala
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2",  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD "," HE ","1HH1","2HH1","1HH2","2HH2"), # arg
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HD2","2HD2",  None,  None,  None,  None,  None,  None,  None), # asn
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD2"," OD1",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None,  None), # asp
    (" N  "," CA "," C  "," O  "," CB "," SG ",  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ",  None,  None,  None,  None,  None,  None,  None,  None), # cys
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE2","2HE2",  None,  None,  None,  None,  None), # gln
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE2"," OE1",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ",  None,  None,  None,  None,  None,  None,  None), # glu
    (" N  "," CA "," C  "," O  ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  ","1HA ","2HA ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None), # gly
    (" N  "," CA "," C  "," O  "," CB "," CG "," NE2"," CD2"," CE1"," ND1",  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HD2"," HE1"," HE2",  None,  None,  None,  None,  None,  None), # his
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1",  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG2","2HG2","3HG2","1HG1","2HG1","1HD1","2HD1","3HD1",  None,  None), # ile
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ","1HD1","2HD1","3HD1","1HD2","2HD2","3HD2",  None,  None), # leu
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ ",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ","1HE ","2HE ","1HZ ","2HZ ","3HZ "), # lys
    (" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE ",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE ","2HE ","3HE ",  None,  None,  None,  None), # met
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD2"," CD1"," CE2"," CE1"," CZ ",  None,  None,  None," H  "," HD2"," HE2"," HZ "," HE1"," HD1"," HA ","1HB ","2HB ",  None,  None,  None,  None), # phe
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD ",  None,  None,  None,  None,  None,  None,  None," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ",  None,  None,  None,  None,  None,  None), # pro
    (" N  "," CA "," C  "," O  "," CB "," OG ",  None,  None,  None,  None,  None,  None,  None,  None," H  "," HG "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None), # ser
    (" N  "," CA "," C  "," O  "," CB "," OG1"," CG2",  None,  None,  None,  None,  None,  None,  None," H  "," HG1"," HA "," HB ","1HG2","2HG2","3HG2",  None,  None,  None,  None,  None,  None), # thr
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2"," H  "," HA ","1HB ","2HB "," HD1"," HE1"," HZ2"," HH2"," HZ3"," HE3",  None,  None,  None), # trp
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD2"," CD1"," CE2"," CE1"," CZ "," OH ",  None,  None," H  "," HA ","1HB ","2HB "," HD2"," HE2"," HE1"," HD1"," HH ",  None,  None,  None,  None), # tyr
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2",  None,  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG1","2HG1","3HG1","1HG2","2HG2","3HG2",  None,  None,  None,  None), # val
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # unk
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # mask
]

aabonds=[
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB ","1HB "),(" CB ","2HB "),(" CB ","3HB ")) , # ala
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," CD "),(" CG ","1HG "),(" CG ","2HG "),(" CD "," NE "),(" CD ","1HD "),(" CD ","2HD "),(" NE "," CZ "),(" NE "," HE "),(" CZ "," NH1"),(" CZ "," NH2"),(" NH1","1HH1"),(" NH1","2HH1"),(" NH2","1HH2"),(" NH2","2HH2")) , # arg
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," OD1"),(" CG "," ND2"),(" ND2","1HD2"),(" ND2","2HD2")) , # asn
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," OD1"),(" CG "," OD2")) , # asp
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," SG "),(" CB ","1HB "),(" CB ","2HB "),(" SG "," HG ")) , # cys
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," CD "),(" CG ","1HG "),(" CG ","2HG "),(" CD "," OE1"),(" CD "," NE2"),(" NE2","1HE2"),(" NE2","2HE2")) , # gln
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," CD "),(" CG ","1HG "),(" CG ","2HG "),(" CD "," OE1"),(" CD "," OE2")) , # glu
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA ","1HA "),(" CA ","2HA "),(" C  "," O  ")) , # gly
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," ND1"),(" CG "," CD2"),(" ND1"," CE1"),(" CD2"," NE2"),(" CD2"," HD2"),(" CE1"," NE2"),(" CE1"," HE1"),(" NE2"," HE2")) , # his
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG1"),(" CB "," CG2"),(" CB "," HB "),(" CG1"," CD1"),(" CG1","1HG1"),(" CG1","2HG1"),(" CG2","1HG2"),(" CG2","2HG2"),(" CG2","3HG2"),(" CD1","1HD1"),(" CD1","2HD1"),(" CD1","3HD1")) , # ile
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," CD1"),(" CG "," CD2"),(" CG "," HG "),(" CD1","1HD1"),(" CD1","2HD1"),(" CD1","3HD1"),(" CD2","1HD2"),(" CD2","2HD2"),(" CD2","3HD2")) , # leu
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," CD "),(" CG ","1HG "),(" CG ","2HG "),(" CD "," CE "),(" CD ","1HD "),(" CD ","2HD "),(" CE "," NZ "),(" CE ","1HE "),(" CE ","2HE "),(" NZ ","1HZ "),(" NZ ","2HZ "),(" NZ ","3HZ ")) , # lys
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," SD "),(" CG ","1HG "),(" CG ","2HG "),(" SD "," CE "),(" CE ","1HE "),(" CE ","2HE "),(" CE ","3HE ")) , # met
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," CD1"),(" CG "," CD2"),(" CD1"," CE1"),(" CD1"," HD1"),(" CD2"," CE2"),(" CD2"," HD2"),(" CE1"," CZ "),(" CE1"," HE1"),(" CE2"," CZ "),(" CE2"," HE2"),(" CZ "," HZ ")) , # phe
    ((" N  "," CA "),(" N  "," CD "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," CD "),(" CG ","1HG "),(" CG ","2HG "),(" CD ","1HD "),(" CD ","2HD ")) , # pro
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," OG "),(" CB ","1HB "),(" CB ","2HB "),(" OG "," HG ")) , # ser
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," OG1"),(" CB "," CG2"),(" CB "," HB "),(" OG1"," HG1"),(" CG2","1HG2"),(" CG2","2HG2"),(" CG2","3HG2")) , # thr
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," CD1"),(" CG "," CD2"),(" CD1"," NE1"),(" CD1"," HD1"),(" CD2"," CE2"),(" CD2"," CE3"),(" NE1"," CE2"),(" NE1"," HE1"),(" CE2"," CZ2"),(" CE3"," CZ3"),(" CE3"," HE3"),(" CZ2"," CH2"),(" CZ2"," HZ2"),(" CZ3"," CH2"),(" CZ3"," HZ3"),(" CH2"," HH2")) , # trp
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG "),(" CB ","1HB "),(" CB ","2HB "),(" CG "," CD1"),(" CG "," CD2"),(" CD1"," CE1"),(" CD1"," HD1"),(" CD2"," CE2"),(" CD2"," HD2"),(" CE1"," CZ "),(" CE1"," HE1"),(" CE2"," CZ "),(" CE2"," HE2"),(" CZ "," OH "),(" OH "," HH ")) , # tyr
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB "," CG1"),(" CB "," CG2"),(" CB "," HB "),(" CG1","1HG1"),(" CG1","2HG1"),(" CG1","3HG1"),(" CG2","1HG2"),(" CG2","2HG2"),(" CG2","3HG2")), # val
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB ","1HB "),(" CB ","2HB "),(" CB ","3HB ")) , # unk
    ((" N  "," CA "),(" N  "," H  "),(" CA "," C  "),(" CA "," CB "),(" CA "," HA "),(" C  "," O  "),(" CB ","1HB "),(" CB ","2HB "),(" CB ","3HB ")) , # mask
]

aa2type = [
    ("Nbb", "CAbb","CObb","OCbb","CH3",   None,  None,  None,  None,  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hapo",  None,  None,  None,  None,  None,  None,  None,  None), # ala
    ("Nbb", "CAbb","CObb","OCbb","CH2", "CH2", "CH2", "NtrR","aroC","Narg","Narg",  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hpol","Hpol","Hpol","Hpol","Hpol"), # arg
    ("Nbb", "CAbb","CObb","OCbb","CH2", "CNH2","ONH2","NH2O",  None,  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hpol","Hpol",  None,  None,  None,  None,  None,  None,  None), # asn
    ("Nbb", "CAbb","CObb","OCbb","CH2", "COO", "OOC", "OOC",   None,  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo",  None,  None,  None,  None,  None,  None,  None,  None,  None), # asp
    ("Nbb", "CAbb","CObb","OCbb","CH2", "SH1",   None,  None,  None,  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","HS",    None,  None,  None,  None,  None,  None,  None,  None), # cys
    ("Nbb", "CAbb","CObb","OCbb","CH2", "CH2", "CNH2","ONH2","NH2O",  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hapo","Hapo","Hpol","Hpol",  None,  None,  None,  None,  None), # gln
    ("Nbb", "CAbb","CObb","OCbb","CH2", "CH2", "COO", "OOC", "OOC",   None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hapo","Hapo",  None,  None,  None,  None,  None,  None,  None), # glu
    ("Nbb", "CAbb","CObb","OCbb",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None), # gly
    ("Nbb", "CAbb","CObb","OCbb","CH2", "CH0", "Nhis","aroC","aroC","Ntrp",  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hpol","Hapo","Hapo",  None,  None,  None,  None,  None,  None), # his
    ("Nbb", "CAbb","CObb","OCbb","CH1", "CH2", "CH3", "CH3",   None,  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo",  None,  None), # ile
    ("Nbb", "CAbb","CObb","OCbb","CH2", "CH1", "CH3", "CH3",   None,  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo",  None,  None), # leu
    ("Nbb", "CAbb","CObb","OCbb","CH2", "CH2", "CH2", "CH2", "Nlys",  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hpol","Hpol","Hpol"), # lys
    ("Nbb", "CAbb","CObb","OCbb","CH2", "CH2", "S",   "CH3",   None,  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo",  None,  None,  None,  None), # met
    ("Nbb", "CAbb","CObb","OCbb","CH2", "CH0", "aroC","aroC","aroC","aroC","aroC",  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Haro","Haro","Haro","Haro","Haro",  None,  None,  None,  None), # phe
    ("Npro","CAbb","CObb","OCbb","CH2", "CH2", "CH2",   None,  None,  None,  None,  None,  None,  None,"Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo",  None,  None,  None,  None,  None,  None), # pro
    ("Nbb", "CAbb","CObb","OCbb","CH2", "OH",    None,  None,  None,  None,  None,  None,  None,  None,"HNbb","Hpol","Hapo","Hapo","Hapo",  None,  None,  None,  None,  None,  None,  None,  None), # ser
    ("Nbb", "CAbb","CObb","OCbb","CH1", "OH",  "CH3",   None,  None,  None,  None,  None,  None,  None,"HNbb","Hpol","Hapo","Hapo","Hapo","Hapo","Hapo",  None,  None,  None,  None,  None,  None), # thr
    ("Nbb", "CAbb","CObb","OCbb","CH2", "CH0", "aroC","CH0", "Ntrp","CH0", "aroC","aroC","aroC","aroC","HNbb","Haro","Hapo","Hapo","Hapo","Hpol","Haro","Haro","Haro","Haro",  None,  None,  None), # trp
    ("Nbb", "CAbb","CObb","OCbb","CH2", "CH0", "aroC","aroC","aroC","aroC","CH0", "OHY",   None,  None,"HNbb","Haro","Haro","Haro","Haro","Hapo","Hapo","Hapo","Hpol",  None,  None,  None,  None), # tyr
    ("Nbb", "CAbb","CObb","OCbb","CH1", "CH3", "CH3",   None,  None,  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo","Hapo",  None,  None,  None,  None), # val
    ("Nbb", "CAbb","CObb","OCbb","CH3",   None,  None,  None,  None,  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hapo",  None,  None,  None,  None,  None,  None,  None,  None), # unk
    ("Nbb", "CAbb","CObb","OCbb","CH3",   None,  None,  None,  None,  None,  None,  None,  None,  None,"HNbb","Hapo","Hapo","Hapo","Hapo",  None,  None,  None,  None,  None,  None,  None,  None), # mask
]

# tip atom
aa2tip = [
        " CB ", # ala
        " CZ ", # arg
        " ND2", # asn
        " CG ", # asp
        " SG ", # cys
        " NE2", # gln
        " CD ", # glu
        " CA ", # gly
        " NE2", # his
        " CD1", # ile
        " CG ", # leu
        " NZ ", # lys
        " SD ", # met
        " CZ ", # phe
        " CG ", # pro
        " OG ", # ser
        " OG1", # thr
        " CH2", # trp
        " OH ", # tyr
        " CB ", # val
        " CB ", # unknown (gap etc)
        " CB " # masked
        ]


torsions=[
    [ None, None, None, None ],  # ala
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," CD "], [" CB "," CG "," CD "," NE "], [" CG "," CD "," NE "," CZ "] ],  # arg
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," OD1"], None, None ],  # asn
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," OD1"], None, None ],  # asp
    [ [" N  "," CA "," CB "," SG "], [" CA "," CB "," SG "," HG "], None, None ],  # cys
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," CD "], [" CB "," CG "," CD "," OE1"], None ],  # gln
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," CD "], [" CB "," CG "," CD "," OE1"], None ],  # glu
    [ None, None, None, None ],  # gly
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," ND1"], [" CD2"," CE1"," HE1"," NE2"], None ],  # his (protonation handled as a pseudo-torsion)
    [ [" N  "," CA "," CB "," CG1"], [" CA "," CB "," CG1"," CD1"], None, None ],  # ile
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," CD1"], None, None ],  # leu
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," CD "], [" CB "," CG "," CD "," CE "], [" CG "," CD "," CE "," NZ "] ],  # lys
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," SD "], [" CB "," CG "," SD "," CE "], None ],  # met
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," CD1"], None, None ],  # phe
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," CD "], [" CB "," CG "," CD ","1HD "], None ],  # pro
    [ [" N  "," CA "," CB "," OG "], [" CA "," CB "," OG "," HG "], None, None ],  # ser
    [ [" N  "," CA "," CB "," OG1"], [" CA "," CB "," OG1"," HG1"], None, None ],  # thr
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," CD1"], None, None ],  # trp
    [ [" N  "," CA "," CB "," CG "], [" CA "," CB "," CG "," CD1"], [" CE1"," CZ "," OH "," HH "], None ],  # tyr
    [ [" N  "," CA "," CB "," CG1"], None, None, None ],  # val
    [ None, None, None, None ],  # unk
    [ None, None, None, None ],  # mask
]

# ideal N, CA, C initial coordinates
init_N = torch.tensor([-0.5272, 1.3593, 0.000]).float()
init_CA = torch.zeros_like(init_N)
init_C = torch.tensor([1.5233, 0.000, 0.000]).float()
INIT_CRDS = torch.full((27, 3), np.nan)
INIT_CRDS[:3] = torch.stack((init_N, init_CA, init_C), dim=0) # (3,3)

norm_N = init_N / (torch.norm(init_N, dim=-1, keepdim=True) + 1e-5)
norm_C = init_C / (torch.norm(init_C, dim=-1, keepdim=True) + 1e-5)
cos_ideal_NCAC = torch.sum(norm_N*norm_C, dim=-1) # cosine of ideal N-CA-C bond angle

#fd Rosetta ideal coords
#fd   - uses same "frame-building" as AF2
ideal_coords = [
    [ # 0 ala
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3341, -0.4928,  0.9132)],
        [' CB ', 8, (-0.5289,-0.7734,-1.1991)],
        ['1HB ', 8, (-0.1265, -1.7863, -1.1851)],
        ['2HB ', 8, (-1.6173, -0.8147, -1.1541)],
        ['3HB ', 8, (-0.2229, -0.2744, -2.1172)],
    ],
    [ # 1 arg
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3467, -0.5055,  0.9018)],
        [' CB ', 8, (-0.5042,-0.7698,-1.2118)],
        ['1HB ', 4, ( 0.3635, -0.5318,  0.8781)],
        ['2HB ', 4, ( 0.3639, -0.5323, -0.8789)],
        [' CG ', 4, (0.6396,1.3794, 0.000)],
        ['1HG ', 5, (0.3639, -0.5139,  0.8900)],
        ['2HG ', 5, (0.3641, -0.5140, -0.8903)],
        [' CD ', 5, (0.5492,1.3801, 0.000)],
        ['1HD ', 6, (0.3637, -0.5135,  0.8895)],
        ['2HD ', 6, (0.3636, -0.5134, -0.8893)],
        [' NE ', 6, (0.5423,1.3491, 0.000)],
        [' NH1', 7, (0.2012,2.2965, 0.000)],
        [' NH2', 7, (2.0824,1.0030, 0.000)],
        [' CZ ', 7, (0.7650,1.1090, 0.000)],
        [' HE ', 7, (0.4701,-0.8955, 0.000)],
        ['1HH1', 7, (-0.8059,2.3776, 0.000)],
        ['1HH2', 7, (2.5160,0.0898, 0.000)],
        ['2HH1', 7, (0.7745,3.1277, 0.000)],
        ['2HH2', 7, (2.6554,1.8336, 0.000)],
    ],
    [ # 2 asn
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3233, -0.4967,  0.9162)],
        [' CB ', 8, (-0.5341,-0.7799,-1.1874)],
        ['1HB ', 4, ( 0.3641, -0.5327,  0.8795)],
        ['2HB ', 4, ( 0.3639, -0.5323, -0.8789)],
        [' CG ', 4, (0.5778,1.3881, 0.000)],
        [' ND2', 5, (0.5839,-1.1711, 0.000)],
        [' OD1', 5, (0.6331,1.0620, 0.000)],
        ['1HD2', 5, (1.5825, -1.2322, 0.000)],
        ['2HD2', 5, (0.0323, -2.0046, 0.000)],
    ],
    [ # 3 asp
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3233, -0.4967,  0.9162)],
        [' CB ', 8, (-0.5162,-0.7757,-1.2144)],
        ['1HB ', 4, ( 0.3639, -0.5324,  0.8791)],
        ['2HB ', 4, ( 0.3640, -0.5325, -0.8792)],
        [' CG ', 4, (0.5926,1.4028, 0.000)],
        [' OD1', 5, (0.5746,1.0629, 0.000)],
        [' OD2', 5, (0.5738,-1.0627, 0.000)],
    ],
    [ # 4 cys
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3481, -0.5059,  0.9006)],
        [' CB ', 8, (-0.5046,-0.7727,-1.2189)],
        ['1HB ', 4, ( 0.3639, -0.5324,  0.8791)],
        ['2HB ', 4, ( 0.3638, -0.5322, -0.8787)],
        [' SG ', 4, (0.7386,1.6511, 0.000)],
        [' HG ', 5, (0.1387,1.3221, 0.000)],
    ],
    [ # 5 gln
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3363, -0.5013,  0.9074)],
        [' CB ', 8, (-0.5226,-0.7776,-1.2109)],
        ['1HB ', 4, ( 0.3638, -0.5323,  0.8789)],
        ['2HB ', 4, ( 0.3638, -0.5322, -0.8788)],
        [' CG ', 4, (0.6225,1.3857, 0.000)],
        ['1HG ', 5, ( 0.3531, -0.5156,  0.8931)],
        ['2HG ', 5, ( 0.3531, -0.5156, -0.8931)],
        [' CD ', 5, (0.5788,1.4021, 0.000)],
        [' NE2', 6, (0.5908,-1.1895, 0.000)],
        [' OE1', 6, (0.6347,1.0584, 0.000)],
        ['1HE2', 6, (1.5825, -1.2525, 0.000)],
        ['2HE2', 6, (0.0380, -2.0229, 0.000)],
    ],
    [ # 6 glu
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3363, -0.5013,  0.9074)],
        [' CB ', 8, (-0.5197,-0.7737,-1.2137)],
        ['1HB ', 4, ( 0.3638, -0.5323,  0.8789)],
        ['2HB ', 4, ( 0.3638, -0.5322, -0.8788)],
        [' CG ', 4, (0.6287,1.3862, 0.000)],
        ['1HG ', 5, ( 0.3531, -0.5156,  0.8931)],
        ['2HG ', 5, ( 0.3531, -0.5156, -0.8931)],
        [' CD ', 5, (0.5850,1.3849, 0.000)],
        [' OE1', 6, (0.5752,1.0618, 0.000)],
        [' OE2', 6, (0.5741,-1.0635, 0.000)],
    ],
    [ # 7 gly
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        ['1HA ', 0, ( -0.3676, -0.5329,  0.8771)],
        ['2HA ', 0, ( -0.3674, -0.5325, -0.8765)],
    ],
    [ # 8 his
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3299, -0.5180,  0.9001)],
        [' CB ', 8, (-0.5163,-0.7809,-1.2129)],
        ['1HB ', 4, ( 0.3640, -0.5325,  0.8793)],
        ['2HB ', 4, ( 0.3637, -0.5321, -0.8786)],
        [' CG ', 4, (0.6016,1.3710, 0.000)],
        [' CD2', 5, (0.8918,-1.0184, 0.000)],
        [' CE1', 5, (2.0299,0.8564, 0.000)],
        [' HE1', 5, (2.8542, 1.5693,  0.000)],
        [' HD2', 5, ( 0.6584, -2.0835, 0.000) ],
        [' ND1', 6, (-1.8631, -1.0722,  0.000)],
        [' NE2', 6, (-1.8625,  1.0707, 0.000)],
        [' HE2', 6, (-1.5439,  2.0292, 0.000)],
    ],
    [ # 9 ile
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3405, -0.5028,  0.9044)],
        [' CB ', 8, (-0.5140,-0.7885,-1.2184)],
        [' HB ', 4, (0.3637, -0.4714,  0.9125)],
        [' CG1', 4, (0.5339,1.4348,0.000)],
        [' CG2', 4, (0.5319,-0.7693,-1.1994)],
        ['1HG2', 4, (1.6215, -0.7588, -1.1842)],
        ['2HG2', 4, (0.1785, -1.7986, -1.1569)],
        ['3HG2', 4, (0.1773, -0.3016, -2.1180)],
        [' CD1', 5, (0.6106,1.3829, 0.000)],
        ['1HG1', 5, (0.3637, -0.5338,  0.8774)],
        ['2HG1', 5, (0.3640, -0.5322, -0.8793)],
        ['1HD1', 5, (1.6978,  1.3006, 0.000)],
        ['2HD1', 5, (0.2873,  1.9236, -0.8902)],
        ['3HD1', 5, (0.2888, 1.9224, 0.8896)],
    ],
    [ # 10 leu
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.525, -0.000, -0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3435, -0.5040,  0.9027)],
        [' CB ', 8, (-0.5175,-0.7692,-1.2220)],
        ['1HB ', 4, ( 0.3473, -0.5346,  0.8827)],
        ['2HB ', 4, ( 0.3476, -0.5351, -0.8836)],
        [' CG ', 4, (0.6652,1.3823, 0.000)],
        [' CD1', 5, (0.5083,1.4353, 0.000)],
        [' CD2', 5, (0.5079,-0.7600,1.2163)],
        [' HG ', 5, (0.3640, -0.4825, -0.9075)],
        ['1HD1', 5, (1.5984,  1.4353, 0.000)],
        ['2HD1', 5, (0.1462,  1.9496, -0.8903)],
        ['3HD1', 5, (0.1459, 1.9494, 0.8895)],
        ['1HD2', 5, (1.5983, -0.7606,  1.2158)],
        ['2HD2', 5, (0.1456, -0.2774,  2.1243)],
        ['3HD2', 5, (0.1444, -1.7871,  1.1815)],
    ],
    [ # 11 lys
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3335, -0.5005,  0.9097)],
        ['1HB ', 4, ( 0.3640, -0.5324,  0.8791)],
        ['2HB ', 4, ( 0.3639, -0.5324, -0.8790)],
        [' CB ', 8, (-0.5259,-0.7785,-1.2069)],
        ['1HG ', 5, (0.3641, -0.5229,  0.8852)],
        ['2HG ', 5, (0.3637, -0.5227, -0.8841)],
        [' CG ', 4, (0.6291,1.3869, 0.000)],
        [' CD ', 5, (0.5526,1.4174, 0.000)],
        ['1HD ', 6, (0.3641, -0.5239,  0.8848)],
        ['2HD ', 6, (0.3638, -0.5219, -0.8850)],
        [' CE ', 6, (0.5544,1.4170, 0.000)],
        [' NZ ', 7, (0.5566,1.3801, 0.000)],
        ['1HE ', 7, (0.4199, -0.4638,  0.9482)],
        ['2HE ', 7, (0.4202, -0.4631, -0.8172)],
        ['1HZ ', 7, (1.6223, 1.3980, 0.0658)],
        ['2HZ ', 7, (0.2970,  1.9326, -0.7584)],
        ['3HZ ', 7, (0.2981, 1.9319, 0.8909)],
    ],
    [ # 12 met
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3303, -0.4990,  0.9108)],
        ['1HB ', 4, ( 0.3635, -0.5318,  0.8781)],
        ['2HB ', 4, ( 0.3641, -0.5326, -0.8795)],
        [' CB ', 8, (-0.5331,-0.7727,-1.2048)],
        ['1HG ', 5, (0.3637, -0.5256,  0.8823)],
        ['2HG ', 5, (0.3638, -0.5249, -0.8831)],
        [' CG ', 4, (0.6298,1.3858,0.000)],
        [' SD ', 5, (0.6953,1.6645,0.000)],
        [' CE ', 6, (0.3383,1.7581,0.000)],
        ['1HE ', 6, (1.7054,  2.0532, -0.0063)],
        ['2HE ', 6, (0.1906,  2.3099, -0.9072)],
        ['3HE ', 6, (0.1917, 2.3792, 0.8720)],
    ],
    [ # 13 phe
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3303, -0.4990,  0.9108)],
        ['1HB ', 4, ( 0.3635, -0.5318,  0.8781)],
        ['2HB ', 4, ( 0.3641, -0.5326, -0.8795)],
        [' CB ', 8, (-0.5150,-0.7729,-1.2156)],
        [' CG ', 4, (0.6060,1.3746, 0.000)],
        [' CD1', 5, (0.7078,1.1928, 0.000)],
        [' CD2', 5, (0.7084,-1.1920, 0.000)],
        [' CE1', 5, (2.0900,1.1940, 0.000)],
        [' CE2', 5, (2.0897,-1.1939, 0.000)],
        [' CZ ', 5, (2.7809, 0.000, 0.000)],
        [' HD1', 5, (0.1613, 2.1362, 0.000)],
        [' HD2', 5, (0.1621, -2.1360, 0.000)],
        [' HE1', 5, (2.6335,  2.1384, 0.000)],
        [' HE2', 5, (2.6344, -2.1378, 0.000)],
        [' HZ ', 5, (3.8700, 0.000, 0.000)],
    ],
    [ # 14 pro
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' HA ', 0, (-0.3868, -0.5380,  0.8781)],
        ['1HB ', 4, ( 0.3762, -0.5355,  0.8842)],
        ['2HB ', 4, ( 0.3762, -0.5355, -0.8842)],
        [' CB ', 8, (-0.5649,-0.5888,-1.2966)],
        [' CG ', 4, (0.3657,1.4451,0.0000)],
        [' CD ', 5, (0.3744,1.4582, 0.0)],
        ['1HG ', 5, (0.3798, -0.5348,  0.8830)],
        ['2HG ', 5, (0.3798, -0.5348, -0.8830)],
        ['1HD ', 6, (0.3798, -0.5348,  0.8830)],
        ['2HD ', 6, (0.3798, -0.5348, -0.8830)],
    ],
    [ # 15 ser
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3425, -0.5041,  0.9048)],
        ['1HB ', 4, ( 0.3637, -0.5321,  0.8786)],
        ['2HB ', 4, ( 0.3636, -0.5319, -0.8782)],
        [' CB ', 8, (-0.5146,-0.7595,-1.2073)],
        [' OG ', 4, (0.5021,1.3081, 0.000)],
        [' HG ', 5, (0.2647, 0.9230, 0.000)],
    ],
    [ # 16 thr
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3364, -0.5015,  0.9078)],
        [' HB ', 4, ( 0.3638, -0.5006,  0.8971)],
        ['1HG2', 4, ( 1.6231, -0.7142, -1.2097)],
        ['2HG2', 4, ( 0.1792, -1.7546, -1.2237)],
        ['3HG2', 4, ( 0.1808, -0.2222, -2.1269)],
        [' CB ', 8, (-0.5172,-0.7952,-1.2130)],
        [' CG2', 4, (0.5334,-0.7239,-1.2267)],
        [' OG1', 4, (0.4804,1.3506,0.000)],
        [' HG1', 5, (0.3194,  0.9056, 0.000)],
    ],
    [ # 17 trp
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3436, -0.5042,  0.9031)],
        ['1HB ', 4, ( 0.3639, -0.5323,  0.8790)],
        ['2HB ', 4, ( 0.3638, -0.5322, -0.8787)],
        [' CB ', 8, (-0.5136,-0.7712,-1.2173)],
        [' CG ', 4, (0.5984,1.3741, 0.000)],
        [' CD1', 5, (0.8151,1.0921, 0.000)],
        [' CD2', 5, (0.8753,-1.1538, 0.000)],
        [' CE2', 5, (2.1865,-0.6707, 0.000)],
        [' CE3', 5, (0.6541,-2.5366, 0.000)],
        [' NE1', 5, (2.1309,0.7003, 0.000)],
        [' CH2', 5, (3.0315,-2.8930, 0.000)],
        [' CZ2', 5, (3.2813,-1.5205, 0.000)],
        [' CZ3', 5, (1.7521,-3.3888, 0.000)],
        [' HD1', 5, (0.4722, 2.1252,  0.000)],
        [' HE1', 5, ( 2.9291,  1.3191,  0.000)],
        [' HE3', 5, (-0.3597, -2.9356,  0.000)],
        [' HZ2', 5, (4.3053, -1.1462,  0.000)],
        [' HZ3', 5, ( 1.5712, -4.4640,  0.000)],
        [' HH2', 5, ( 3.8700, -3.5898,  0.000)],
    ],
    [ # 18 tyr
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3305, -0.4992,  0.9112)],
        ['1HB ', 4, ( 0.3642, -0.5327,  0.8797)],
        ['2HB ', 4, ( 0.3637, -0.5321, -0.8785)],
        [' CB ', 8, (-0.5305,-0.7799,-1.2051)],
        [' CG ', 4, (0.6104,1.3840, 0.000)],
        [' CD1', 5, (0.6936,1.2013, 0.000)],
        [' CD2', 5, (0.6934,-1.2011, 0.000)],
        [' CE1', 5, (2.0751,1.2013, 0.000)],
        [' CE2', 5, (2.0748,-1.2011, 0.000)],
        [' OH ', 5, (4.1408, 0.000, 0.000)],
        [' CZ ', 5, (2.7648, 0.000, 0.000)],
        [' HD1', 5, (0.1485, 2.1455,  0.000)],
        [' HD2', 5, (0.1484, -2.1451,  0.000)],
        [' HE1', 5, (2.6200, 2.1450,  0.000)],
        [' HE2', 5, (2.6199, -2.1453,  0.000)],
        [' HH ', 6, (0.3190, 0.9057,  0.000)],
    ],
    [ # 19 val
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3497, -0.5068,  0.9002)],
        [' CB ', 8, (-0.5105,-0.7712,-1.2317)],
        [' CG1', 4, (0.5326,1.4252, 0.000)],
        [' CG2', 4, (0.5177,-0.7693,1.2057)],
        [' HB ', 4, (0.3541, -0.4754, -0.9148)],
        ['1HG1', 4, (1.6228,  1.4063,  0.000)],
        ['2HG1', 4, (0.1790,  1.9457, -0.8898)],
        ['3HG1', 4, (0.1798, 1.9453, 0.8903)],
        ['1HG2', 4, (1.6073, -0.7659,  1.1989)],
        ['2HG2', 4, (0.1586, -0.2971,  2.1203)],
        ['3HG2', 4, (0.1582, -1.7976,  1.1631)],
    ],
    [ # 20 unk
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3341, -0.4928,  0.9132)],
        [' CB ', 8, (-0.5289,-0.7734,-1.1991)],
        ['1HB ', 8, (-0.1265, -1.7863, -1.1851)],
        ['2HB ', 8, (-1.6173, -0.8147, -1.1541)],
        ['3HB ', 8, (-0.2229, -0.2744, -2.1172)],
    ],
    [ # 21 mask
        [' N  ', 0, (-0.5272, 1.3593, 0.000)],
        [' CA ', 0, (0.000, 0.000, 0.000)],
        [' C  ', 0, (1.5233, 0.000, 0.000)],
        [' O  ', 3, (0.6303, 1.0574, 0.000)],
        [' H  ', 2, (0.4920,-0.8821,  0.0000)],
        [' HA ', 0, (-0.3341, -0.4928,  0.9132)],
        [' CB ', 8, (-0.5289,-0.7734,-1.1991)],
        ['1HB ', 8, (-0.1265, -1.7863, -1.1851)],
        ['2HB ', 8, (-1.6173, -0.8147, -1.1541)],
        ['3HB ', 8, (-0.2229, -0.2744, -2.1172)],
    ],
]
