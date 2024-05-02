import re

def score(*aminiacid, file='blosum62.txt'):

    with open(file) as obj:
        content = obj.read()

        blosum = content.split('\n')
        amino_acids = [x for x in blosum[0] if x.isalpha()]

        matriz = []
        for i in blosum:
            values = re.findall(r'-?\d+(?:\.\d+)?', i)
            if values:
                matriz.append(values)
        row, colunm = 0,0

        print(aminiacid)
        for aa in enumerate(amino_acids):
            if aa[1] == aminiacid[0]:
                row = aa[0]
            if aa[1] == aminiacid[1]:
                colunm = aa[0]
    
        return float(matriz[row][colunm])
    
score('D', 'D')