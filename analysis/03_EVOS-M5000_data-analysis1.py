import pandas as pd
import shared.dataframe as dat
from datetime import datetime

# input parameters
exp_name = '20211210_EVOS-M5000_CRISPR-KO_nucleofectionAndLipofectionTest'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20211210_CRISPR-KO_nucleofectionAndLipofectionTest/"
num_per_treatment = 6

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
data_FOV['treatment'] = [data_FOV['sample'][i][:-5] for i in range(len(data_FOV))]
data_FOV = data_FOV.sort_values(by='sample', ignore_index=True)

# extract information for generating summary data dataframe
data_treatment = list()
data_int = list()
data_number = list()
for i in range(int(len(data_FOV)/num_per_treatment)):
    data_treatment.append(data_FOV['treatment'][num_per_treatment*i])
    temp_int = []
    for j in range(num_per_treatment):
        temp_int = temp_int + data_FOV['intensity'][num_per_treatment*i+j]
    data_int.append(temp_int)
    data_number.append(len(temp_int))
data = pd.DataFrame({'treatment': data_treatment, 'intensity': data_int, 'number': data_number})

# extract information about controls
f.write("Maximum intensity in control group is used as a cutoff due to relative small size of this experiment.\n")
ctrl = max(list(data[data['treatment'] == 'ctrl']['intensity'])[0])
ctrl_count = int(data[data['treatment'] == 'ctrl']['number'])

# add essential data to data dataframe
data['cutoff'] = [ctrl]*int(len(data))
data['ctrl_count'] = [ctrl_count]*int(len(data))
data['positive'] = [[x for x in data['intensity'][i] if x > data['cutoff'][i]] for i in range(len(data))]
data['number_positive'] = [len(data['positive'][i]) for i in range(len(data))]
data['percentage_positive'] = data['number_positive']*100.0/data['number']

save_path = master_folder
data_FOV.to_csv('%ssummary_FOV.txt' % save_path, index=False, sep='\t')
data.to_csv('%ssummary_treatment.txt' % save_path, index=False, sep='\t')
data.to_excel('%ssummary_treatment.xlsx' % save_path, index=False)

print("DONE!")
f.write("DONE!\n\n")
f.close()
