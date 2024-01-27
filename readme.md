# Sign Language Dissertation 2023-24

## Todo

- output the openpose to json format
- preprocess the json format of body/face/hand tracking to calculate point changes between each frames (this will be the new metrics)
- write up for introduction and literture review part. 


## Files
Folder - `giorgia_code`: files Giorgia sent on Oct 20, 2023's email,
contains: 
- `segment_and_preprocess_EEG.py` EEG segmentation and preprocessing code
- `IVC2.py`code for extracting the IVC or 'envelope' from the video (based on grayscale - must be extended to RGB) *(waht does this mean)*
- code for generating the CND data structure
- code for running the cross-validation
- code for running some basic analysis (this is a bit messy and not well documented)
- `result_without_hand_face` folder contains json files of body pose detection from open pose. Every file is a a frame in the video. 8000+ files.

### Files in the EEG dataset
.BDF files: BioSemi Data Format

## Code and implemnentation part

1. **IVC Implementation - Compare and Contrast**

2. **Human Pose Estimation**

   - [Read more](https://www.analyticsvidhya.com/?url=https%3A%2F%2Fwww.analyticsvidhya.com%2Fblog%2F2022%2F01%2Fa-comprehensive-guide-on-human-pose-estimation%2F)
   - A Comprehensive Guide on Human Pose Estimation ([analyticsvidhya.com](https://www.analyticsvidhya.com/blog/2022/01/a-comprehensive-guide-on-human-pose-estimation/))

    Pose detection used Open pose as the solution 


## Notes for Presentation

1. **Apple HDR Animation**

   - [Watch the video](https://www.youtube.com/live/KR0g-1hnQPA?si=MJJxT0utcJuKFp1p&t=3254)

2. **3D Rotate Test PowerPoint**
   
3. **PNG to GLB Conversion**
   - Used an online tool to convert PNG to GLB format.