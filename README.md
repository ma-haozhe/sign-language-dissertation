# Neural Correlates of Sign Language Movements: A Multimodal EEG and Computer Vision Study

**MSc Research Dissertation | Trinity College Dublin (2023-24)**
**Author: Haozhe Ma**
**Supervisor: Prof. Giovanni Di Liberto (Trinity College Dublin)**

A dissertation project investigating the relationship between sign language movements and brain activity using multimodal neuroimaging and computer vision techniques.

## Overview

This research explores how the brain processes sign language movements by combining two powerful approaches:

1. **Electroencephalography (EEG)** - Capturing real-time brain electrical activity with high temporal resolution
2. **Computer Vision** - Quantifying movement dynamics through pose estimation and motion analysis

The goal is to understand how neural responses correlate with the kinematic features of sign language gestures, potentially informing both cognitive neuroscience research and assistive communication technologies.

## Research Problem

Sign language represents a unique window into human communication and cognition, serving as a fully developed natural language that uses visual-spatial gestures rather than spoken words. Understanding how the brain processes these movements can provide insights into:

- The neural basis of language processing in the visual modality
- How movement complexity influences cortical responses
- Potential applications in brain-computer interfaces for communication

## Research Questions

1. **Can movement kinematics of sign language be quantified using computer vision?** - Using OpenPose for body, face, and hand keypoint detection to extract meaningful movement features.

2. **Do EEG responses correlate with quantified movement features?** - Investigating whether Temporal Response Functions (TRF) can predict neural activity from movement-based stimulus envelopes.

3. **How do semantic content vs. motor complexity affect neural processing?** - Comparing neural responses between LSF (Linguistic Semantic Features) and FR (Fingerspelling Recognition) conditions.

## Methodology

### 1. EEG Data Acquisition and Processing

**Data Collection:**
- 64-channel EEG recordings using BioSemi ActiveTwo system
- 22 participants: 6 fluent LSF signers (Deaf) + 16 French non-signers (hearing)
- Participants viewed sign language video stimuli during recording
- Two experimental conditions: Linguistic Semantic Features (LSF) and Fingerspelling Recognition (FR)

**Preprocessing Pipeline:**
```
Raw BDF Files
    -> Notch Filter (50Hz)
    -> Bandpass Filter (1-4Hz, 1-8Hz, 1-30Hz)
    -> Rereferencing (Mastoids)
    -> Downsampling (64Hz)
    -> Segmentation by stimulus onset
```

**Analysis Approach:**
- **Temporal Response Functions (TRF)** modeling using mTRF Toolbox
- Global Field Power (GFP) analysis for response magnitude quantification
- Forward and reverse encoding models to capture neural dynamics

### 2. Computer Vision for Motion Quantification

**Pose Estimation:**
- OpenPose library for body, face, and hand keypoint detection
- 25-body, 70-face, and 21-hand keypoints per frame
- Confidence-based filtering for robust tracking

**Movement Quantification - IVC (Intensity Vector Changes):**
```
IVC = Sum of squared differences between consecutive frame keypoints
    = Σ[(x_i(t+1) - x_i(t))² + (y_i(t+1) - y_i(t))²]
```

This metric captures the total movement energy across all detected body points, providing a continuous measure of gesture dynamics.

### 3. Data Integration

The IVC time series served as the stimulus envelope for TRF modeling, enabling direct comparison between movement kinematics and neural responses.

## Results

### Key Findings

1. **Successful Neural-Movement Correlation:** The study demonstrated that EEG responses can be modeled using the IVC stimulus envelope derived from computer vision analysis.

2. **Temporal Dynamics:** Global Field Power (GFP) peaks were observed at approximately **150-280ms latency** following stimulus onset, consistent with classical ERP findings in language processing.

3. **Frequency-Band Specificity:** Different frequency bands (1-4Hz, 1-8Hz, 1-30Hz) revealed varying strengths of movement-neural correlations, suggesting frequency-dependent processing of movement features.

4. **Condition-Dependent Responses:** LSF and FR conditions showed differential neural patterns, indicating that semantic content and motor complexity influence cortical processing.

5. **Original vs. Time-Reversed:** Neural responses to original playback showed distinct peaks that were absent or altered in time-reversed conditions, confirming sensitivity to temporal order of movements.

### Key Results with Numbers

| Finding | Metric | Value | Interpretation |
|---------|--------|-------|----------------|
| Sample Size | Participants | N = 22 | 6 LSF signers + 16 non-signers |
| GFP Peak Latency | Time (ms) | ~150-280ms | Consistent with language ERP literature |
| EEG Channels | System | 64-channel | BioSemi ActiveTwo |
| Sampling Rate | Original | Up to 16kHz | Downsampled to 64Hz |
| Trials per Participant | Conditions | 14 trials | 11 correct playback, 3 time-reversed |
| Frequency Bands | Analysis | 1-4Hz, 1-8Hz, 1-30Hz | Delta to Theta range |

### Experimental Design Summary

- **Stimuli:** 14 French Sign Language (LSF) videos, 2-8 minutes each
- **Video Resolution:** 1280x720 pixels at 30 fps
- **Control Condition:** Time-reversed videos to test temporal order sensitivity
- **Statistical Tests:** Wilcoxon Signed-Rank (within-subject), Wilcoxon Rank-Sum (between conditions)

## Visualizations

### 1. Research Framework

![Forward TRF Sign Language](plots/forwardTRF-sign-language.png)

*Overview of the Temporal Response Function (TRF) analysis framework for sign language processing. This pipeline illustrates the multimodal approach combining behavioral features with neural data.*

---

### 2. EEG Analysis: Butterfly Plots & GFP

#### LSF Condition (1-8Hz)
![Butterfly GFP LSF](plots/IVC/1-8Hz/LSF/butterfly-and-GFP.png)

*EEG butterfly plot with Global Field Power (GFP) for the LSF (Linguistic Semantic Features) condition. Left panels show overlaid responses across all 64 channels; the black trace represents GFP, highlighting temporal dynamics of neural responses.*

#### FR Condition (1-8Hz)
![Butterfly GFP FR](plots/IVC/1-8Hz/FR/butterfly-and-GFP.png)

*EEG butterfly plot with GFP for the FR (Fingerspelling Recognition) condition. Comparing LSF and FR reveals differences in neural response magnitude and timing.*

#### GFP Comparison: Original vs. Time-Reversed
![GFP Comparison](plots/IVC/1-8Hz/GFP-comparison.png)

*Global Field Power comparison across conditions. Top row: LSF group (Original vs. Time-Reversed playback). Bottom row: FR group. Clear peaks at ~150ms in original conditions demonstrate neural encoding of movement direction.*

---

### 3. TRF Weight Analysis

#### TRF Weights - LSF (1-30Hz)
![TRF Weights LSF](plots/IVC/1-30Hz/LSF/weights.png)

*Temporal Response Function weights across EEG channels and time lags for LSF condition. Heat map showing neural encoding strength at different post-stimulus latencies.*

#### Average TRF Response
![Average TRF](plots/IVC/1-30Hz/LSF/avgTRF.png)

*Grand-averaged TRF response across all channels and time lags, demonstrating consistent neural encoding patterns.*

#### Fz-Cz-Pz Electrodes
![Fz-Cz-Pz TRF](plots/IVC/1-30Hz/LSF/Fz-Cz-Pz.png)

*TRF time course extracted from centro-parietal electrode cluster (Fz, Cz, Pz), showing response peaks at specific latencies.*

#### Topographic Maps Over Time
![Topography Series](plots/IVC/1-30Hz/LSF/3.png)
![Topography Series](plots/IVC/1-30Hz/LSF/4.png)
![Topography Series](plots/IVC/1-30Hz/LSF/5.png)
![Topography Series](plots/IVC/1-30Hz/LSF/6.png)
![Topography Series](plots/IVC/1-30Hz/LSF/7.png)
![Topography Series](plots/IVC/1-30Hz/LSF/8.png)

*Spatiotemporal evolution of neural activation across the scalp. Each topographic map shows the spatial distribution of TRF weights at different time lags post-stimulus onset.*

---

### 4. Topographic Analysis

#### Topoplot Average - LSF (1-8Hz)
![Topoplot Average LSF](plots/IVC/1-8Hz/LSF/topoplot-avg.png)

*Average topographic map across the response window for LSF condition, showing spatial distribution of neural activity.*

#### Topoplot Weights - FR (1-8Hz)
![Topoplot Weights FR](plots/IVC/1-8Hz/FR/topoplot-weights.png)

*Topographic weight map for FR condition, comparing spatial activation patterns.*

---

### 5. IVC (Intensity Vector Changes) Analysis

![IVC Analysis Overview](plots/IVC.png)

*IVC methodology overview showing the correlation between OpenPose-derived movement features and neural responses. The video panel shows frame-by-frame movement quantification.*

#### IVC Analysis - FR (1-4Hz)
![IVC FR 1-4Hz](plots/IVC/1-4Hz/FR/IVC%20analysis.png)

*Intensity Vector Changes analysis for FR condition in the 1-4Hz frequency band, showing feature-by-feature correlations with neural data.*

#### TRF Weights Topography Comparison (1-8Hz)
![TRF Weights Topography IVC](plots/IVC/1-8Hz/TRF-weights-topoplot-comparison-IVC.png)

*Comparison of TRF weight topographies across IVC-derived features, showing which movement components drive neural responses.*

---

### 6. Statistical Comparisons

#### Prediction Correlation Across Conditions
![Prediction Correlation](plots/IVC/1-8Hz/prediction-correlation-all-conditions.jpg)

*Box plots comparing prediction correlations (r-values) across four conditions: FR-V, FR-R, LSF-V, LSF-R. The LSF-V condition shows the highest predictive performance, indicating strongest neural-movement correlation.*

#### Wilcoxon Signed-Rank Test: GFP
![Wilcoxon Signed Rank Test](plots/IVC/1-8Hz/wilcoxon-signed-rank-test.png)

*Statistical comparison of GFP between Original and Time-Reversed playback conditions. Top panel: LSF group. Bottom panel: FR group. Demonstrates significant differences in neural processing based on temporal order.*

#### Wilcoxon Rank-Sum Test
![Wilcoxon Rank Sum Test](plots/IVC/1-8Hz/wilcoxon-rank-sum-test.png)

*Non-parametric statistical test comparing distributions across experimental conditions.*

---

### 7. OpenPose Computer Vision Visualizations

#### Motion Tracking with ROI
![Motion Tracking ROI](IVC%20frame%20with%20ROI.png)

*Video frame with Region of Interest (ROI) bounding box highlighting the detected signer. OpenPose keypoints overlay the body, face, and hands, enabling precise motion tracking for IVC calculation.*

#### Frame with ROI
![Frame with ROI](frame_with_ROI.png)

*Annotated video frame showing OpenPose skeleton detection and ROI definition for motion analysis.*

#### Gesture Progression Analysis
![Gesture Progression](plots/Slide1.png)

*Three-frame sequence showing sign language gesture progression with motion tracking overlay. This visualization demonstrates how continuous movement is captured and quantified for correlation with neural data.*

---

### 8. OpenPose IVC Results

#### TRF and GFP Comparison - OpenPose IVC
![TRF GFP OpenPose](plots/openpose_IVC/TRF-and-GFP-comparison-openposeIVC.png)

*Temporal Response Functions and GFP computed using OpenPose-derived stimulus envelope. Comparison between original and time-reversed conditions.*

#### TRF Weights Topography - OpenPose IVC
![TRF Weights Topography OpenPose](plots/openpose_IVC/TRF-weights-topography-comparison-openposeIVC.png)

*Topographic comparison of TRF weights when using OpenPose features as stimulus envelope.*

#### FR Condition - OpenPose IVC
![OpenPose IVC FR](plots/openpose_IVC/FR/weight_and_avg.jpg)

*OpenPose IVC analysis for Fingerspelling Recognition condition showing weight distributions and average responses.*

---

## Project Structure

```
sign-language-dissertation-tcd-2023/
├── EEG_raw/                 # Raw EEG data in BioSemi BDF format
├── giorgia_code/           # EEG preprocessing and analysis scripts
│   ├── segment_and_preprocess_EEG.py
│   ├── IVC2.py
│   └── npy2resampledmat.py
├── plots/                  # Generated figures and visualizations
│   ├── IVC/                # IVC analysis results
│   │   ├── 1-4Hz/
│   │   ├── 1-8Hz/
│   │   └── 1-30Hz/
│   ├── openpose_IVC/       # OpenPose-based IVC results
│   └── forwardTRF-sign-language.png
├── openpose-output/        # OpenPose detection results
├── result_without_hand_face/  # JSON keypoint files from OpenPose
├── IVC_for_keypoints.py    # IVC calculation from OpenPose data
├── TCD-dissertation/       # Dissertation documents
└── stimuli/                # Sign language video stimuli
```

## Getting Started

### Prerequisites

- **Python 3.8+**
- **MATLAB** (for mTRF Toolbox)
- **MNE-Python** for EEG preprocessing
- **OpenPose** for pose estimation
- **NumPy, SciPy, pandas** for data processing

### Installation

1. Clone the repository:
```bash
git clone https://github.com/ma-haozhe/sign-language-dissertation.git
```

2. Install Python dependencies:
```bash
pip install mne numpy scipy pandas matplotlib opencv-python
```

3. Set up the mTRF Toolbox (included in CNSP-resources):
```bash
cd CNSP-resources/CNSP/libs/mTRF-Toolbox_v2
```

### Running the Analysis

**Preprocess EEG data:**
```bash
cd giorgia_code
python segment_and_preprocess_EEG.py
```

**Compute IVC from video:**
```bash
python IVC_for_keypoints.py
```

## Technologies Used

| Technology | Purpose |
|------------|---------|
| **MNE-Python** | EEG data preprocessing and analysis |
| **mTRF Toolbox** | Temporal Response Function modeling |
| **OpenPose** | Body, face, and hand pose estimation |
| **NumPy/SciPy** | Numerical computations and signal processing |
| **Matplotlib** | Data visualization |
| **BioSemi BDF** | EEG data format |
| **MATLAB** | Statistical analysis and mTRF toolbox |

## Key Achievements

- Designed and implemented a multimodal framework combining EEG and computer vision
- Developed novel methodology for quantifying sign language movements using IVC metrics
- Achieved correlation between movement kinematics and neural responses
- Validated findings across multiple frequency bands and experimental conditions
- Created reproducible preprocessing and analysis pipelines
- Demonstrated neural sensitivity to temporal order of sign language movements

## Limitations & Future Work

### Limitations

- **Sample Size:** N = 22 participants (6 LSF signers + 16 non-signers), short of target of 50
- **Group Imbalance:** Uneven distribution between signer and non-signer groups may affect statistical power
- **Single Language:** Study focused on French Sign Language (LSF); findings may not generalize to other sign languages
- **Pilot Study:** Initial data collection phase; results are preliminary
- **Linguistic Depth:** Analysis focused on movement features rather than deeper linguistic properties

### Future Directions

- **Complete Data Collection:** Finish recruiting to reach target sample size with balanced groups
- **Linguistic Transcription:** Integrate word-by-word transcription to analyze linguistic features
- **Separate Body Part Analysis:** Independent analyses for hand movements vs. facial expressions
- **Extended Frequency Analysis:** Further investigate theta band (4-8Hz) correlations
- **Cross-Linguistic Studies:** Replicate with ASL, BSL, or other sign languages
- **BCI Applications:** Apply findings to develop communication technologies for deaf individuals

## License

This project is part of an academic dissertation. Please contact the author for usage permissions.

## Acknowledgments

- Trinity College Dublin - School of Computer Science and Statistics
- Computational Neuroscience and Psychophysiology (CNSP) resources
- mTRF Toolbox developers for providing the temporal response function analysis framework
