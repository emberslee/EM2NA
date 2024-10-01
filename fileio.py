import os
def readlines(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    return lines
