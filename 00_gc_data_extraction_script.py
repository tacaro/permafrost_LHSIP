import numpy as np
import pandas as pd
import os

# A quick little python script for iteratively parsing .RES files from SRI GC and outputting the results as a .csv

# FUNCTION DEFINITION

def read_GC_RES_file(file_name, file_num):
    print("Reading file path: " + file_name)

    with open(file_name) as file:
        lines = file.readlines()

    sample_type = lines[34]
    column = lines[18]
    carrier = lines[20]

    file.close() # close the file

    file = open(file_name) # open the file again

    cleaned_array = [] # create an empty array

    for i, line in enumerate(file):

        if i < 42: # skip lines that have irrerlevant data
            continue

        else:
            line_array = line.strip().split(',')

            if len(line_array) == 1:

                # Skipping the rows that have an empty list (empty lines in the file)
                continue

            # List comprehension on each element, stripping it of quotes. Doesn't matter if it doesn't have quotes.
            new_array = [x.strip('\"') for x in line_array]

            # Append your new array to the storage array.
            cleaned_array.append(new_array)

    try:
        # Write the results
        data_array = np.array(cleaned_array)
        data_array = data_array[:,:5]
        # Write the relevant metadata
        df = pd.DataFrame(data_array, columns=['gas','retention_time','peak_area','peak_height','correction_factor'])
        df['file_number'] = file_num
        df['file_name'] = file_name.replace('\n', '') # add the file name, removing newlines
        df['sample_type'] = sample_type.replace('\n', '') # add the sample type (entered manually in program), removing newlines
        df['column'] = column.replace('\n', '') # add column info, removing newlines
        df['carrier'] = carrier.replace('\n', '') # add carrier info, removing newlines
        return df
    except IndexError:
        print("Data file: " + file_name + " encountered an error. Check the file!")
        pass
    except ValueError:
        print("Data file: " + file_name + " encountered an error. Check the file!")
        pass


# EXECUTE

GC_file_directory = "data/GC_Data/raw_data/" # If your GC data is not in a file like this, change this string!
GC_files = [f for f in os.listdir(GC_file_directory) if f.endswith(".RES") and "FID" in f]

appended_data = []

n = 1
for res in GC_files:
    data = read_GC_RES_file(GC_file_directory + res, file_num = n)
    n = n + 1
    # store DataFrame in a list
    appended_data.append(data)

# See pd.concat documentation for more info on this operation
appended_data = pd.concat(appended_data)
# write DataFrame to csv
appended_data.to_csv("data/GC_Data/GC_Data_Table3.csv")
