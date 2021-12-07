import pandas as pd
import shared.dataframe as dat
from datetime import datetime

# input parameters
exp_name = '20211207_EVOS-M5000_ColoDM-Cas9_nucleofectionTest'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20211207_nucleofection-GFPtest/"
cutoff = 13.72636103

# write log file
f = open("%slog.txt" % master_folder, "a+")
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
f.write("exp_name: %s\n" % exp_name)
f.write("datetime: %s\n" % dt_string)
f.write("code_running: %s\n" % __file__)
f.write("master_folder: %s\n" % master_folder)
f.write("Notes:\n")
f.write("Analyze intensity file to generate summary files grouped by either individual FOV or treatment.\n")

# load data
data_FOV = pd.read_csv('%ssummary_intensity.txt' % master_folder, na_values=['.'], sep='\t')
# reformat data
data_FOV['intensity'] = [dat.str_to_float(data_FOV['intensity'][i]) for i in range(len(data_FOV))]
data_FOV['number'] = [len(data_FOV['intensity'][i]) for i in range(len(data_FOV))]
# calculate positive
data_FOV = data_FOV.sort_values(by='sample', ignore_index=True)
data_FOV['cutoff'] = [cutoff]*int(len(data_FOV))
data_FOV['positive'] = [[x for x in data_FOV['intensity'][i] if x > data_FOV['cutoff'][i]] for i in range(len(data_FOV))]
data_FOV['number_positive'] = [len(data_FOV['positive'][i]) for i in range(len(data_FOV))]
data_FOV['percentage_positive'] = data_FOV['number_positive']*100.0/data_FOV['number']

save_path = master_folder
data_FOV.to_csv('%ssummary_FOV.txt' % save_path, index=False, sep='\t')

print("DONE!")
f.write("DONE!\n\n")
f.close()
