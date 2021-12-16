import os
import shutil
from datetime import datetime

# input parameters
exp_name = '20211210_EVOS-M5000_CRISPR-KO_nucleofectionAndLipofectionTest'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20211210_CRISPR-KO_nucleofectionAndLipofectionTest/"
data_source = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20211210_CRISPR-KO_nucleofectionAndLipofectionTest/"
ch = ['DAPI', 'RFP', 'TRANS']

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
files = [x for x in os.listdir(data_source)]
for i in range(len(ch)):
    move_path = "%s%s/" % (data_source, ch[i])
    if not os.path.exists(move_path):
        os.makedirs(move_path)

f.write("moving files...\n")
for j in files:
    n = 0
    for i in ch:
        if i in j:
            shutil.move(("%s%s" % (data_source, j)), ("%s%s/" % (data_source, i)))
            n = 1
    if n == 0:
        shutil.move(("%s%s" % (data_source, j)), ("%smerge/" % data_source))

print("DONE!")
f.write("DONE!\n\n")
f.close()
