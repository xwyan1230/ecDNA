import os
import shutil
from datetime import datetime

# input parameters
exp_name = '20211022_EVOS-M5000_ColoDM-Cas9_nucleofectionTest'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20211022_ColoDM-Cas9_nucleofectionTest/"
data_source = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20211022_ColoDM-Cas9_nucleofectionTest/EVOS_M5000/"

# log all the running info
if not os.path.exists("%slog.txt" % master_folder):
    f = open("%slog.txt" % master_folder, "w+")
else:
    f = open("%slog.txt" % master_folder, "a+")

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
f.write("exp_name: %s\n" % exp_name)
f.write("datetime: %s\n" % dt_string)
f.write("code_running: %s\n" % __file__)
f.write("master_folder: %s\n" % master_folder)
f.write("data_source: %s\n" % data_source)
f.write("Notes:\n")

# script start
f.write("creating moving folders...\n")
dirs = [x for x in os.listdir(data_source)]
move_path_GFP = ("%sGFP/" % data_source)
if not os.path.exists(move_path_GFP):
    os.makedirs(move_path_GFP)
move_path_TRANS = ("%sTRANS/" % data_source)
if not os.path.exists(move_path_TRANS):
    os.makedirs(move_path_TRANS)
move_path_merge = ("%smerge/" % data_source)
if not os.path.exists(move_path_merge):
    os.makedirs(move_path_merge)

f.write("moving files...\n")
for i in dirs:
    if 'GFP' in i:
        shutil.move(("%s%s" % (data_source, i)), move_path_GFP)
    elif 'TRANS' in i:
        shutil.move(("%s%s" % (data_source, i)), move_path_TRANS)
    else:
        shutil.move(("%s%s" % (data_source, i)), move_path_merge)

print("DONE!")
f.write("DONE!\n\n")
f.close()