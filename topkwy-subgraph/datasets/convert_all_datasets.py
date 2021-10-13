import math
import numpy as np
import sys
import os
import re
import time


datasets = ["AIDS","BZR","COX2","DHFR","Mutagenicity","DD","ENZYMES","MUTAG","NCI1","NCI109","Tox21_AHR"]


for dataset in datasets:

    cmd = "python convert_dataset.py "+str(dataset)+"/"+str(dataset)
    os.system(cmd)
