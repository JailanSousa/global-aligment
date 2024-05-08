import re
import numpy as np

from visualize_align import wrap_text

class GlobalAligment():

    """Perform a global amino acid sequence aligment 
       based on needleman-wunsch algorithm.
    """
    
    def __init__(self, query, subject, gap=-1):

        self.query = query
        self.subject = subject
        self.seq1 = ""
        self.seq2 = ""
        self.m = 0
        self.n = 0
        self.gap = gap

    def read_fastas(self):
        """Read fasta files and return their sequences."""

        fastas = [self.query, self.subject]

        sequences = []
        for fasta in fastas:
            with open(fasta) as obj:
                seqs = obj.readlines()
                seq = "".join(seqs[1:]).replace('\n', '')
                sequences.append(seq)
        
        self.seq1 = sequences[0]
        self.seq2 = sequences[1]
        self.m = len(sequences[0])
        self.n = len(sequences[1])

        return sequences


    def score_system(*aminiacid, file='blosum62.txt'):

        """Load blosum62 matriz and return amino acid substitution score."""

        with open(file) as obj:
            content = obj.read()

            blosum = content.split('\n')
            amino_acids = [x for x in blosum[0] if x.isalpha()] # get the 20 amino acids.

            # get scores
            matriz = []
            for i in blosum:
                values = re.findall(r'-?\d+(?:\.\d+)?', i) # get only numbers.
                if values:
                    matriz.append(values)

            # Get score by accesing amino acids positions.
            row, column = 0, 0
            for aa in enumerate(amino_acids):
                if aa[1] == aminiacid[1]:
                    row = aa[0]
                if aa[1] == aminiacid[2]:
                    column = aa[0]

            return float(matriz[row][column])

    def matrix_filling(self):

        """Return two matrizes: the first one is a dynamic programming matriz
           and the second one is the backtracking pointers information.
        """
        seq1 = self.read_fastas()[0]
        seq2 = self.read_fastas()[1]

        # Initializing dynamic programming matriz.
        matrix = np.zeros((self.m + 1, self.n + 1))
        traceback = []

        # backtracking pointers
        for i in range(self.m+1):
            temp = []
            for j in range(self.n+1):
                temp.append(0)
            traceback.append(temp)

        # Initializing the first rows and columns of both matrizes.
        for i in range(self.m + 1):
            matrix[i, 0] = i * self.gap # multiplayng the firts colunm by gap penalty
            traceback[i][0] = 'up' # fulling the first colunm with up.

        for j in range(self.n + 1):
            matrix[0, j] = j * self.gap # multiplayng the firts row by gap penalty
            traceback[0][j]= 'left' # fulling the first colunm with left.
        
        traceback[0][0] = 'Done'

        # filling matrizes.
        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                # getting amino acid substitution score form blosum62 matriz.
                score = self.score_system(seq1[i-1], seq2[j-1])
                # filling dynamic programming matriz.
                sij = [matrix[i,j-1]+self.gap, matrix[i-1,j]+self.gap, matrix[i-1,j-1] + score]
                matrix[i,j] = max(sij)

                # getting backtracking pointers.
                if max(sij) == sij[2]:
                    traceback[i][j] = 'diag'
                elif max(sij) == sij[1]:
                    traceback[i][j] = 'up'
                else:
                    traceback[i][j] = 'left'
            
        return [matrix, traceback]

    def backtracking(self):

        """Return a list with both sequences aligned."""

        # backtracking pointers
        traceback = self.matrix_filling()[1]

        i = self.m
        j = self.n
        aligment = [[],[]]


        while (i > 0 or j > 0):
            # put both amino acids in aligment.
            if traceback[i][j] == 'diag':
                aligment[0].append(self.seq1[i-1])
                aligment[1].append(self.seq2[j-1])
                i = i-1
                j = j-1
            
            elif traceback[i][j] == 'left':
                aligment[0].append('-') # put a gap in the first sequence.
                aligment[1].append(self.seq2[j-1])
                j = j-1
            
            elif traceback[i][j] == 'up':
                aligment[0].append(self.seq1[i-1])
                aligment[1].append('-') # put a gap in the second sequence.
                i = i-1
            elif traceback[i][j] == 'Done':
                break
        
        return aligment
    
    def show_matriz(self):

        """Show the dynamic programming table."""
        
        matriz_filling = self.matrix_filling()

        score_matriz = matriz_filling[0].tolist()
        traceback = matriz_filling[1]

        # Inicializating a new matriz.
        matriz = []
        for i in range(self.m + 1):
            temp = []
            for j in range(self.n + 1):
                temp.append(0)
            matriz.append(temp)

        # filling the new matriz.
        for i in range(self.m + 1):
            for j in range(self.n + 1):
                
                # getting scores and backtracking pointers.
                if traceback[i][j] == 'diag':
                    matriz[i][j] = f"{score_matriz[i][j]}{u'\u2196'}" 
                                                        # Unicode northwest arrow caracter.
                elif traceback[i][j] == 'left':
                    matriz[i][j] = f"{score_matriz[i][j]}{u'\u2190'}"
                                                        # Unicode left arrow caracter.
                elif traceback[i][j] == 'up':
                    matriz[i][j] = f"{score_matriz[i][j]}{u'\u2191'}"
                                                        # Unicode up arrow caracter.
                else:
                     matriz[i][j] = '*'

        # add the first sequence in the first column.
        for i in range(1, self.m + 1):
            matriz[i][0] = self.seq1[i-1]

        # add the first sequence in the first row.
        for j in range(1, self.n + 1):
            matriz[0][j] = self.seq2[j-1]

        # Showing dynamic programming table.
        for i in matriz:
            for j in i:
                print(j, end='\t')
            print('\n')
    
    def show_aligment(self):

        """Show aligment, score aligment and identity."""

        backtracking = self.backtracking()
        seq1 = ''.join(reversed(backtracking[0]))
        seq2 = ''.join(reversed(backtracking[1]))

        alig_score = 0
        caracters = ""
        for i, j in zip(seq1, seq2):
            if i == j:
                alig_score += self.score_system(i,j)
                caracters += '|' # match caracter
            else:
                if i == '-' or j == '-':
                    alig_score += self.gap
                    caracters += " " # when gap
                else:
                    alig_score += self.score_system(i,j)
                    caracters += '*' # mismach
        
        # Calculating identity.
        matchs = caracters.count('|')
        lenght = self.m

        identity = 100 * (matchs / lenght)

        # Showing aligment.
        algn1 = wrap_text(seq1)
        char = wrap_text(caracters)
        alig2 = wrap_text(seq2)

        for i,j,k in zip(algn1,char,alig2):
            print(f"seq1 {i}")
            print(f"     {j}")
            print(f"seq2 {k}\n")

        # Showing aligment metrics.
        print(f"\nNeedleman-Wunsh aligment score: {alig_score}")
        print(f"Indentity: {identity:.1f}%")
    
    def run(self):
        """Show results."""
        print(f"\n{20*'---'}Dinamic programming table.{20*'---'}")
        self.show_matriz()
        print(20*'----')

        print("Sequence Aligment")
        self.show_aligment()


    
w = "seq1.fasta"
v = "seq2.fasta"

alig = GlobalAligment(w, v)
alig.run()