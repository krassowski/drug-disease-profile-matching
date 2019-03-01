from pathlib import Path
from os import system
from warnings import warn

TEMP_ROOT = '/tmp/drug_repositioning/'


def mount_root_in_ramdisk(root):

    path = Path(root)

    if not path.expanduser().exists() or not path.is_mount():
        creation_status = system(f'mkdir -p {root}')
        assert creation_status == 0
        mount_command = f'sudo mount -t tmpfs -o size=2G tmpfs {root}'
        mounting_status = system(mount_command)
        if mounting_status != 0:
            warn(
                f'Failed to mount {root} temporary file system.\n'
                f'Please run: {mount_command}\n'
                'and reload the offending module.'
            )


def create_tmp_dir(name, root=TEMP_ROOT):

    mount_root_in_ramdisk(root)
    path = Path(root) / name

    if not path.exists():
        system(f'mkdir -p {str(path)}')

    return str(path)
