from snpy import fset
import os
path_to_filters = fset['B'].__dict__['file'].split('filters/filters/')[0]
if 'rvgps' in path_to_filters: raise Exception('Dont update for dev env')
print (path_to_filters)


tarfile = 'complete_filters.tgz'
os.system(f"tar -czvf {tarfile} filters")

print (f'Untarring {tarfile} at: {path_to_filters}')
os.system(f"tar -xvf {tarfile} --directory {path_to_filters}")

os.system(f"rm {tarfile}")
