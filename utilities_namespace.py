import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
from data_frames import MyDataFrame as DataFrame
from helpers.gui.table import display_table
from helpers.gui.source import embed_source_styling, show_source

from IPython.display import HTML
from functools import reduce
from pandas import read_table, read_csv, concat, Series
import numpy as np
import seaborn as sns
from matplotlib import pyplot
from typing import Iterable
from copy import deepcopy, copy
from mpl_toolkits.mplot3d import Axes3D
from itertools import chain
from helpers.gui import Namespace, NeatNamespace, HorizontalNamespace
from helpers.r import *
from tqdm.auto import tqdm
from statistics import mean
from types import SimpleNamespace

MyDataFrame = DataFrame

pd.options.display.max_rows = 10
pd.options.display.max_columns = 10

show_table = display_table


def keys(obj):
    return list(obj.keys())


embed_source_styling()

T = True
F = False
