import json
import numpy as np
from scipy.io import savemat
import os
from glob import glob
import re #regular expression

# NOTE: 
# this file is for process the keypoints detected by open pose
# and do the calculation for frame by frame difference
# similar to IVC, to see the improvements/changes in the envelope we get


# Define paths
input_path = 'result_without_hand_face'  # replace with the path to your JSON files
output_path = 'IVC_openpose_envelope'
os.makedirs(output_path, exist_ok=True)

# Function to extract frame number from filename
def extract_frame_number(filename):
    match = re.search(r'(\d+)_keypoints', filename)
    return int(match.group(1)) if match else -1

# Get the list of JSON files sorted by frame order
json_files = sorted(
    glob(os.path.join(input_path, '*.json')),
    key=lambda x: extract_frame_number(os.path.basename(x))
)

# Initialize the list to store the IVC results
ivc_values = []

# Initialize previous keypoints array
prev_keypoints = None

# Process each JSON file
for i, json_file in enumerate(json_files):
    with open(json_file, 'r') as f:
        data = json.load(f)
        people = data.get('people', [])
        if not people:
            continue  # Skip if no people are detected in the frame

        # Assuming one person per frame, adjust if necessary
        person = people[0]
        keypoints = np.array(person.get('pose_keypoints_2d', []))

        # Reshape keypoints array to (N,3) where N is number of keypoints
        keypoints = keypoints.reshape(-1, 3)

        # Calculate IVC if not the first frame
        if prev_keypoints is not None:
            # Calculate the sum of squared differences
            ivc_value = np.sum((keypoints[:, :2] - prev_keypoints[:, :2])**2)
            ivc_values.append(ivc_value)
            print(f"Processing frame {i}, IVC value: {ivc_value}")

        # Update previous keypoints
        prev_keypoints = keypoints

    # Save the IVC values after each iteration
    np.save(os.path.join(output_path, 'ivc_values.npy'), ivc_values)
    savemat(os.path.join(output_path, 'ivc_values.mat'), {'ivc_values': ivc_values})

# After processing all frames, save the final IVC values
np.save(os.path.join(output_path, 'ivc_values_final.npy'), ivc_values)
savemat(os.path.join(output_path, 'ivc_values_final.mat'), {'ivc_values': ivc_values})
