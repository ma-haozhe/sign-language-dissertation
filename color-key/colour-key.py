import cv2
import numpy as np

video_path = "/Users/haozhema/sign-language-dissertation-tcd-2023/stimuli/V01.mp4"
background_image_path = "/Users/haozhema/sign-language-dissertation-tcd-2023/color-key/black.png"
output_video_path = "/Users/haozhema/sign-language-dissertation-tcd-2023/stimuli/keyed_V01.mp4"  # MP4 for H.264

video = cv2.VideoCapture(video_path)

# Check if video opened successfully
if not video.isOpened():
    print("Error opening video file")

# Load the background image, in this case, it should be black static image.
background_image = cv2.imread(background_image_path)

# Get video properties
frame_width = int(video.get(cv2.CAP_PROP_FRAME_WIDTH))
frame_height = int(video.get(cv2.CAP_PROP_FRAME_HEIGHT))
frame_rate = int(video.get(cv2.CAP_PROP_FPS))

# Stretch the background image to match the video resolution
background_image = cv2.resize(background_image, (frame_width, frame_height))

# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # You can also use 'avc1' for H.264

out = cv2.VideoWriter(output_video_path, fourcc, frame_rate, (frame_width, frame_height))

while True:
    # Read a frame from the video
    ret, frame = video.read()
    if not ret:
        break  # Exit the loop if no more frames are read

    # Define the upper and lower bounds for the green color
    upper_green = np.array([163, 143, 131])
    lower_green = np.array([98, 83, 67])

    # Create a mask to extract the green color from the frame
    mask = cv2.inRange(frame, lower_green, upper_green)

    # Apply the mask to the frame to isolate the green color
    green_parts = cv2.bitwise_and(frame, frame, mask=mask)

    # Calculate the parts of the frame that are not green
    non_green_parts = frame - green_parts

    # Replace the non-green parts with the background image
    final_frame = np.where(non_green_parts == 0, background_image, non_green_parts)

    # Write the frame into the file 'output_video.mp4'
    out.write(final_frame)

    # Display the original video frame and the modified frame with the replaced background
    cv2.imshow("Original Video", frame)
    cv2.imshow("Modified Video", final_frame)

    # Exit the loop when the 'Esc' key is pressed (ASCII code 27)
    if cv2.waitKey(25) == 27:
        break

# Release the video capture object, the video writer object and close all windows
video.release()
out.release()
cv2.destroyAllWindows()
