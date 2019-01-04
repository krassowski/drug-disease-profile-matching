from pandas import read_csv

from config import DATA_DIR


class CancerCensus:

    def __init__(self):
        self.census = read_csv(DATA_DIR + '/cosmic/census.csv')
        self.census.columns = [c.lower().replace(' ', '_') for c in self.census.columns]
