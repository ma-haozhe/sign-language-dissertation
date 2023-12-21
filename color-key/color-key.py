import cv2
import numpy as np

video = cv2.VideoCapture("/Users/haozhema/sign-language-dissertation-tcd-2023/stimuli/V01.mp4")

# Load the background image, in this case, it should be black static image. 
background_image = cv2.imread("/Users/haozhema/sign-language-dissertation-tcd-2023/color-key/black.png")

while True:
    # Read a frame from the video
    ret, frame = video.read()

    # Resize the frame and background image to the same dimensions
    frame = cv2.resize(frame, (640, 480))
    background_image = cv2.resize(background_image, (640, 480))

    # Define the upper and lower bounds for the green color
    upper_green = np.array([163,143,131])
    lower_green = np.array([98, 83, 67])

    # Create a mask to extract the green color from the frame
    mask = cv2.inRange(frame, lower_green, upper_green)

    # Apply the mask to the frame to isolate the green color
    green_parts = cv2.bitwise_and(frame, frame, mask=mask)

    # Calculate the parts of the frame that are not green
    non_green_parts = frame - green_parts

    # Replace the non-green parts with the background image
    final_frame = np.where(non_green_parts == 0, background_image, non_green_parts)

    # Display the original video frame and the modified frame with the replaced background
    cv2.imshow("Original Video", frame)
    cv2.imshow("Modified Video", final_frame)

    # Exit the loop when the 'Esc' key is pressed (ASCII code 27)
    if cv2.waitKey(25) == 27:
        break

# Release the video capture object and close all windows
video.release()
cv2.destroyAllWindows()
