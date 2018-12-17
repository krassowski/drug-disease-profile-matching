from pandas import read_csv


class CancerCensus:

    def __init__(self):
        self.census = read_csv('data/cosmic/census.csv')
        self.census.columns = [c.lower().replace(' ', '_') for c in self.census.columns]
