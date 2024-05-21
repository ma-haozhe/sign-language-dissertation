import cv2
import os
import numpy as np
from scipy.io import savemat

import matplotlib.pyplot as plt

'''
This version of IVC crops the video to a specific 'region of interest'
by using pixel_margin_height and pixel_margin_width.

Calculates the sum of squared differences in grayscale values between successive frames 
and stores these values in a list (envelope).

save in both npy (numpy) and mat format.

Q1 What does it mean by Giorgia "based on grayscale - 
must be extended to RGB"?

Q2 why does it crop the ROI region, what does the top right corner do
'''
# set this flag to enable reveral of envelope data.
reverse_envelopes = True
#input_path = '../../EEG_experiments/SignLanguage_PsychoPy/input/stimuli/'
#output_path = './outputs/features/envelope2/'
input_path = 'stimuli/'
#output_path = 'envelope'
output_path = 'outputs/features/envelope2/'
reverse_output_path = 'outputs/features/envelope1/'
# video 1 8096frames/269sec = 30fps
pixel_margin_heigh = 50
pixel_margin_width = 150

os.makedirs(output_path, exist_ok=True)
os.makedirs(reverse_output_path, exist_ok=True)

# Function to reverse and save envelope data
def reverse_and_save_envelope(filename, data, output_path):
    reversed_data = data[::-1]
    reversed_filename = filename.replace('V', 'R', 1)  # Replace the first occurrence of 'V' with 'R'
    np.save(os.path.join(output_path, reversed_filename + '.npy'), reversed_data)
    savemat(os.path.join(output_path, reversed_filename + '.mat'), {'feature': reversed_data})

for filename in os.listdir(input_path):
    if filename.endswith( ".mp4"):
        filename = filename.split('/')[-1].split('.')[0]
        cam = cv2.VideoCapture(os.path.join(input_path, filename + '.mp4'))
        fps = cam.get(cv2.CAP_PROP_FPS)

        # Init envelope list
        envelope = []

        # Read in just first frame
        ret,imgPrev = cam.read()   

        # Get the patch of interest   
        heigh = imgPrev.shape[0] 
        width = imgPrev.shape[1]
        imgPrev = imgPrev[pixel_margin_heigh:heigh,pixel_margin_width:width-pixel_margin_width,:] 
               
        # Iterate over frames
        f = 1
        while(True):
            ret,imgCur = cam.read()

            # if there remain frames in the video to read in
            if ret:  
                ''' 
                # Draw a rectangle around the ROI
                cv2.rectangle(imgPrev, 
                            (pixel_margin_width, pixel_margin_heigh), 
                            (width - pixel_margin_width, heigh - pixel_margin_heigh), 
                            (0, 0, 255),  # Red color in BGR
                            3)  # Thickness      
                # Display the image with ROI
                plt.imshow(cv2.cvtColor(imgPrev, cv2.COLOR_BGR2RGB))
                plt.title("Frame with ROI")
                plt.show()      
                '''
                
                # Init
                imgDif = np.zeros_like(imgPrev)
                Dif = 0

                # Get the patch
                imgCur = imgCur[pixel_margin_heigh:heigh,pixel_margin_width:width-pixel_margin_width,:]

                # Get the value                
                imgDif  = np.sum(np.sum((imgCur[:,:,0] - imgPrev[:,:,0]) ** 2))   #0 specifies grayscale value
                envelope.append(imgDif)
                np.save(os.path.join(output_path, filename + '.npy'), envelope) 
                savemat(os.path.join(output_path, filename + '.mat'), {'feature': envelope})               
                print (f, imgDif)

                imgPrev = imgCur
                f += 1
            else:
                break

        cam.release()
        cv2.destroyAllWindows()
        np.save(os.path.join(output_path, filename + '.npy'), envelope)
        savemat(os.path.join(output_path, filename + '.mat'), {'feature': envelope})

# Reverse the envelopes if the flag is set to True
if reverse_envelopes:
    for filename in os.listdir(output_path):
        if filename.endswith(".npy"):
            base_filename = filename.split('.')[0]
            envelope_data = np.load(os.path.join(output_path, filename))
            reverse_and_save_envelope(base_filename, envelope_data, reverse_output_path)
        elif filename.endswith(".mat"):
            base_filename = filename.split('.')[0]
            envelope_data = np.load(os.path.join(output_path, base_filename + '.npy'))  # Load the numpy file
            reverse_and_save_envelope(base_filename, envelope_data, reverse_output_path)
            