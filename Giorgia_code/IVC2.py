import cv2
import os
import numpy as np
from scipy.io import savemat

'''
This version of IVC crops the video to a specific 'region of interest'
by using pixel_margin_height and pixel_margin_width.

Calculates the sum of squared differences in grayscale values between successive frames 
and stores these values in a list (envelope).

save in both npy (numpy) and mat format.

What does it mean by Giorgia "based on grayscale - 
must be extended to RGB"?
'''

input_path = '../../EEG_experiments/SignLanguage_PsychoPy/input/stimuli/'
output_path = './outputs/features/envelope2/'
pixel_margin_heigh = 50
pixel_margin_width = 150

os.makedirs(output_path, exist_ok=True)
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
