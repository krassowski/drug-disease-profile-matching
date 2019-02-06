from helpers.temp import create_tmp_dir, mount_root_in_ramdisk

tmp_dir = create_tmp_dir('gsea/input')
gsea_home = mount_root_in_ramdisk('~/gsea_home')
