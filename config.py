from pathlib import Path

pwd = Path(__file__).absolute().parent

DATA_DIR = (pwd / 'data').as_posix()

third_party_dir = pwd / 'thirdparty'
assert third_party_dir.exists()
