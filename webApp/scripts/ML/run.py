#!/usr/bin/python3

##Header
__author__ = "Gurdeep Singh"
__credits__ = ["Gurdeep Singh"]
__license__ = "GPL"
__email__ = "gurdeep330@gmail.com"

import os, sys
from sklearn.metrics import r2_score, roc_curve
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pylab
import make_dic
import make_fm
import logreg_ml
import prediction

## Load all mutation data
dic, protein = make_dic.main()
## Make a training matrix
X, y, X_test, y_test, features, new_test, exclude, sfm, scaler = make_fm.main(dic, protein)
## Make the predictor
clf = logreg_ml.main(X, y, X_test, y_test, features, dic, new_test, exclude)
sys.exit()
## Save the predictions
prediction.main(dic, clf, sfm, scaler)
