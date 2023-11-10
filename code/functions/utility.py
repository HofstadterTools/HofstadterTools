"""Functions for input/output and miscellaneous use."""

import os
import csv
import numpy as np


def read_t_from_file():

    dir = "configuration/" if os.path.isdir('configuration') else ""
    with open(dir+"hopping_input.txt", 'r') as csvfile:
        data = csv.reader(csvfile, delimiter='\t')
        NN, t = [], []
        for i, row in enumerate(data):
            NN.append(int(row[0]))
            t.append(float(row[1]))

    t_list = np.zeros(max(NN))
    for i, val in enumerate(NN):
        t_list[val-1] = t[i]

    return t_list

