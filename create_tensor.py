from os.path import join as jn
import numpy as np
import argparse
import matplotlib.pyplot as plt
from Tensor import MDTensor


parser = argparse.ArgumentParser(description='Create tensor from directories')
parser.add_argument('-d',  type=str, nargs='+',
                    help='List of directories to create the tensors from')
parser.add_argument('-o', type=str, nargs='+',
                    help='List of output names')
args = parser.parse_args()

assert args.d == args.o, 'Number of folder and output names must be the same'

for output, folder in zip(args.o, args.d):
    ts = MDTensor(folder)
    ts.save_tensor(output)






