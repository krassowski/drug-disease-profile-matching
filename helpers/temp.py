from pathlib import Path
from os import system

TEMP_ROOT = '/tmp/drug_repositioning/'


def mount_root_in_ramdisk(root):

    if not Path(root).exists():
        system(f'mkdir -p {root}')
        system(f'sudo mount -t tmpfs -o size=2048M tmpfs {root}')


def create_tmp_dir(name, root=TEMP_ROOT):

    mount_root_in_ramdisk(root)
    path = Path(root) / name

    if not path.exists():
        system(f'mkdir -p {str(path)}')

    return str(path)