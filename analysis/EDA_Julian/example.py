# This is a sample Python script.

# Press Umschalt+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from CELLO_multiplier import compute_embedding,classify,cellcompare
import pandas as pd
import json
import matplotlib.pyplot as plt

normalize = True
cell_label = "osteoblast"
n_LVs = 10

def init(cell_label, n_LVs):

    B = pd.read_csv("B.csv").transpose()

    with open("bulk_labels.json", 'r') as f:
        cell_types_ids = json.load(f)

    Y = pd.read_csv('Ewing_NT_cell_lines.csv', index_col=0)
    mat4 = pd.read_csv('mat4-CellO_train-ewing_genes.csv', index_col=0)
    Z = pd.read_csv('Z.csv')

    b, Y, mat4 = compute_embedding(Y, Z, mat4)
    Y_samples = Y.index.values
    mat4_samples = mat4.index.values
    cell_types = classify(B, b, Y_samples, mat4_samples, cell_types_ids, normalize)
    print(cell_types)

    for cell_type in cell_types[:10].index:
        LVs = cellcompare(cell_types_ids, cell_type, B, n_LVs, b)
        print(LVs)
        fig = LVs.plot.barh( figsize=(40, 40), fontsize=10).get_figure()
        fig.savefig(str(cell_type)+'.pdf')
if __name__ == '__main__':
    init(cell_label, n_LVs)

