import os
import numpy as np
from scipy.io import savemat
from scipy.signal import resample

input_path = './outputs/30Hz/features/envelope2/'
output_path = './outputs/64Hz/features/envelope2/'
os.makedirs(output_path, exist_ok=True)

resample = False
oldfreq = 30
newfreq = 64

# iterate over numpy files in input_path
for filename in os.listdir(input_path):
    if filename.endswith( ".npy"):
        filename = filename.split('/')[-1].split('.')[0]
        print('processing file: ' + filename)

        # Load .npy file
        data = np.load(os.path.join(input_path, filename + '.npy'))
        
        # Resample
        if resample:    
            data = resample(data, int(data.shape[1] / oldfreq * newfreq), axis=1)

        # Save as .mat file
        savemat(os.path.join(output_path, filename + '.mat'), {'feature': data})