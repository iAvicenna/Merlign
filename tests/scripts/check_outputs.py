#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 11:30:29 2023

@author: avicenna
"""

import sys
import pandas as pd

if __name__ == "__main__":

  #if len(sys.argv) != 4:
  #  raise ValueError("Provide two arguments, file1 and file2")

  #file1 = sys.argv[1]
  #file2 = sys.argv[2]

  file1 = "../test_files/simulated_ign_output.csv"
  file2 = "../test_files/simulated_ign_output_test.csv"

  output1 = pd.read_csv(file1, index_col=0, header=0)
  output2 = pd.read_csv(file2, index_col=0, header=0)

  breakpoint()
