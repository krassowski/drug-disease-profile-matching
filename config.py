from pathlib import Path

DATA_DIR = '/home/ubuntu/data/drug_repositioning'
pwd = Path(__file__).absolute().parent
third_party_dir = pwd / 'thirdparty'
assert third_party_dir.exists()
