import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Define the file paths for your two images
image1_path = 'envelope/plot-original-ivc.jpg'
image2_path = 'IVC_openpose_envelope/plot_with_openpose.jpg'

# Load the images
image1 = mpimg.imread(image1_path)
image2 = mpimg.imread(image2_path)

# Create a subplot with 2 rows and 1 column
plt.subplot(2, 1, 1)
plt.imshow(image1)
plt.title(image1_path)

plt.axis('off')

plt.subplot(2, 1, 2)
plt.imshow(image2)
plt.title(image2_path)

plt.axis('off')

# Adjust spacing between plots for readability
plt.tight_layout()

# Display the combined plot
plt.show()