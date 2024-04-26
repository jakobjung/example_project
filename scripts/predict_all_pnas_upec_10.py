# script that uses the trained model in pickle format (./data/models/rf_optimized_model.sav)
# to predict the PNAS class of all UPEC genomes in the dataset

import pandas as pd
import numpy as np
import pickle
from sklearn.ensemble import RandomForestRegressor
