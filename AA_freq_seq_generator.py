import sys, random

# BioE 131 Assignment 4, Michael Vredenburgh --- From command line -h for help
# Analyzes amino acid frequency in a gene, and displays those frequencies, as
# well as generate new sequences with similar frequencies from that of the input sequence.

class InputFASA:
    fileName = ""
    myRawFile = [""]
    myGenes = []
    def __init__(self, myFileName):
        self.fileName = myFileName
        self.myRawFile = InputFASA.extractSequence(myFileName)
        self.parseGene()
    def extractSequence(myFileName):
        # Parses the FASTA file, and creates a string for the sequence.
        infile = open(myFileName)  # Creates a variable to store the opened file to reference later
        myRawFile= infile.readlines()
        return myRawFile
    def parseGene(self):
        # Creates an array of Gene objects called myGenes.
        n=0
        HeaderList = []
        for line in self.myRawFile:
            if line[0][:1]=='>':
                HeaderList.append(n)
            n+=1
        HeaderList.append(len(self.myRawFile)+1)
        n=0
        while n < len(HeaderList)-1:
            self.myGenes.append(Gene(self.myRawFile[HeaderList[n]:HeaderList[n+1]-1]))
            n+=1
                
class Gene:
    header = ""
    mySequence = ""
    length=0
    A=0
    T=0
    G=0
    C=0
    myNewSequence = ""
    def __init__(self, FASTAchunk):
        # FASTAchunk is the header line plus the genetic infomation.
        self.header = FASTAchunk[0]
        self.makeSequence(FASTAchunk)
        self.geneParams()
        self.makeNewSequence()
    def makeSequence(self, FASTAchunk):
        # Makes a long string out of the FASTAarray. Also removes all instances of "\n".
        FASTAgene=FASTAchunk[1:]
        for line in FASTAgene:
            line = line.replace("\n","")
            self.mySequence += line
    def addLength(self):
        myLength = len(self.mySequence)
        self.header = self.header.replace("\n","")
        self.header += " %d bp" % myLength
    def geneParams(self):
        self.length = len(self.mySequence)
        for basepair in self.mySequence:
            self.length+=1
            if basepair == "A":
                self.A+=1
            if basepair == "T":
                self.T+=1
            if basepair == "G":
                self.G+=1
            if basepair == "C":
                self.C+=1
    def makeNewSequence(self):
        self.myNewSequence = newSeq(self.length, self.A, self.T, self.G, self.C).newSeq

class newSeq:
    newSeq = "TAC"
    length= 0
    A = 0
    T = 0
    G = 0
    C = 0
    def __init__(self, length, A, T, G, C):
        self.length = length
        self.A=A
        self.T=T
        self.G=G
        self.C=C
        self.makeSeq()
        self.toString()
    def makeSeq(self):
        numOfCodons = int(self.length/3) # Rounds off to insure multiple of 3
        for i in range(0,numOfCodons):
            self.newSeq += self.genCodon()
        self.newSeq+="TAG"
        #print(self.newSeq)
    def genCodon(self):
        BP=[]
        for i in range(0,self.A):
            BP.append("A")
        for i in range(0,self.T):
            BP.append("T")
        for i in range(0,self.G):
            BP.append("G")
        for i in range(0,self.C):
            BP.append("C")
        codon = random.sample(BP,3)
        codon = codon[0] + codon[1] + codon[2]
        if codon == "TAA": 
            return self.genCodon()
        if codon == "TAG":
            return self.genCodon()
        if codon == "TGA":
            return self.genCodon()
        return codon
    def toString(self):
        print("> Length:%s A:%s T:%s G:%s C:%s" % (self.length, self.A, self.T, self.G, self.C))
        copySeq = self.newSeq
        while True:
            if len(copySeq) <= 80:
                print(copySeq[:80])
                copySeq = copySeq[80:]
            else:
                print(copySeq)
                break
class loadParam:
    length= 0
    A = 0
    T = 0
    G = 0
    C = 0
    def __init__(self, myFileName):
        infile = open(myFileName)
        myRawFile= infile.readlines()
        myRawFileArray = []
        for lines in myRawFile:
            myRawFileArray.append(lines)
        self.length = int(myRawFileArray[0][2:])
        self.A = int(myRawFileArray[1][2:])
        self.T = int(myRawFileArray[2][2:])
        self.G = int(myRawFileArray[3][2:])
        self.C = int(myRawFileArray[4][2:])
        
def inputCMD(argv1,argv2,argv3,argv4):
    if argv1 == "--calc":
        param=InputFASA(argv2)
        newSeq(param.myGenes[0].length, param.myGenes[0].A, param.myGenes[0].T, param.myGenes[0].G, param.myGenes[0].C)
        if argv3 == "--save":
            file = open(argv4,"w")
            paramString = """L:%s
A:%s
T:%s
G:%s
C:%s
""" % (param.myGenes[0].length, param.myGenes[0].A, param.myGenes[0].T, param.myGenes[0].G, param.myGenes[0].C)
            file.write(paramString)
            file.close()
            return param.myGenes[0].myNewSequence
    if argv1 == "--load":
        param=loadParam(argv2)
        mySeq = newSeq(param.length, param.A, param.T, param.G, param.C)
        if argv3 == "--save":
            file = open(argv4,"w")
            paramString = """L:%s
A:%s
T:%s
G:%s
C:%s
""" % (param.length, param.A, param.T, param.G, param.C)
            file.write(paramString)
            file.close()
            return mySeq.newSeq
def aminoAcid(sequence):
    myDict = {"ATT":"I", "ATC":"I", "ATA":"I", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L", "TTA":"L", "TTG":"L", "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V", "TTT":"F", "TTC":"F", "ATG":"M", "TGT":"C", "TGC":"C", "GCT":"A",
              "GCC":"A", "GCA":"A", "GCG":"A", "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G", "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P", "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T", "TCT":"S", "TCC":"S", "TCA":"S", "TCA":"S",
              "TCG":"S", "AGT":"S", "AGC":"S", "TAT":"Y", "TAC":"Y", "TGG":"W", "CAA":"Q", "CAG":"Q", "AAT":"N", "AAC":"N", "CAT":"H", "CAC":"H", "GAA":"E", "GAG":"E", "GAT":"D", "GAC":"D", "AAA":"K", "AAG":"K", "CGT":"R",
              "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R"}
    length = len(sequence)/3
    I=0
    L=0
    V=0
    F=0
    M=0
    C=0
    A=0
    G=0
    P=0
    T=0
    S=0
    Y=0
    W=0
    Q=0
    N=0
    H=0
    E=0
    D=0
    K=0
    R=0
    for i in range(0,len(sequence)-6):
        if i % 3 == 0:
            if myDict[sequence[3+i:6+i]] == "I":
                I+=1
            if myDict[sequence[3+i:6+i]] == "L":
                L+=1
            if myDict[sequence[3+i:6+i]] == "V":
                V+=1
            if myDict[sequence[3+i:6+i]] == "F":
                F+=1
            if myDict[sequence[3+i:6+i]] == "M":
                M+=1
            if myDict[sequence[3+i:6+i]] == "C":
                C+=1
            if myDict[sequence[3+i:6+i]] == "A":
                A+=1
            if myDict[sequence[3+i:6+i]] == "G":
                G+=1
            if myDict[sequence[3+i:6+i]] == "P":
                P+=1
            if myDict[sequence[3+i:6+i]] == "T":
                T+=1
            if myDict[sequence[3+i:6+i]] == "S":
                S+=1
            if myDict[sequence[3+i:6+i]] == "Y":
                Y+=1
            if myDict[sequence[3+i:6+i]] == "W":
                W+=1
            if myDict[sequence[3+i:6+i]] == "Q":
                Q+=1
            if myDict[sequence[3+i:6+i]] == "N":
                N+=1
            if myDict[sequence[3+i:6+i]] == "H":
                H+=1
            if myDict[sequence[3+i:6+i]] == "E":
                E+=1
            if myDict[sequence[3+i:6+i]] == "D":
                D+=1
            if myDict[sequence[3+i:6+i]] == "K":
                K+=1
            if myDict[sequence[3+i:6+i]] == "R":
                R+=1
    AminoAcidString="""
The following is the frequency of each amino acid, represented by their one letter code.
I:%s
L:%s
V:%s
F:%s
M:%s
C:%s
A:%s
G:%s
P:%s
T:%s
S:%s
Y:%s
W:%s
Q:%s
N:%s
H:%s
E:%s
D:%s
K:%s
R:%s""" % ( I/length, L/length, V/length, F/length, M/length, C/length, A/length, G/length, P/length, T/length, S/length, Y/length, W/length, Q/length, N/length, H/length, E/length, D/length, K/length, R/length)
    print(AminoAcidString)
if len(sys.argv) == 6:
    if sys.argv[5] == "--more":
        mySeq = inputCMD(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
        aminoAcid(mySeq)
    else:
        raise("The 5th option is reserved only for --more, to display amino acid distribution.")
if len(sys.argv) == 5:      
    inputCMD(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
if len(sys.argv) == 4:
    raise("You must give a name to your saved configuration.")
if len(sys.argv) == 3:
    inputCMD(sys.argv[1],sys.argv[2],"","") 
if len(sys.argv) == 2:
    if sys.argv[1] == "-h":
        helpInfo = """HELP INFORMATION
(Script by Michael Vredenburgh for BioE131 Assignment 4)
This program generates semi-random sequences of DNA that begin and end with the appropriate start and stop codons. The frequency of base pair expression can be determined in two ways. First, the –calc can be used to calculate the expression from a FASTA DNA file. The second method is to load in a parameters file. The specifics of the parameter file will be covered at the end of the help menu. Statistics about the amino acid distribution can also be displayed.
The following three lines are example command line inputs.
$ python3 Assignment4.py –calc myglobin.txt –save myparams
$ python3 Assignment4.py –load myparams
$ python3 Assignment4.py –load myparams –save myparams2
The previous inputs will not display any information about the amino acid content. This feature was added to try to meet the extra credit requirement. To see the distribution there must be all four options used, plus a fifth –more option. To do that, enter something similar to the following example.
$ python3 Assignment4.py –calc myglobin.txt –save myparams –more
The parameter files are simple .txt files that hold only the following information.
**
L:#
A:#
T:#
G:#
C:#
**
L being length, and the other letters representing the number of respected base pairs.
"""
        print(helpInfo)
    else:
        raise("Not enough information provided. Try entering like this example. Assignment4.py --calc seqfile.fasa --save myparams. For additional help try the help menue, -h.")
if len(sys.argv) == 1:
    raise("Not enough information provided. Try entering like this example. Assignment4.py --calc seqfile.fasa --save myparams. For additional help try the help menue, -h.")
