# Advanced Filtering of EEG Signals for Neurobehavioral Analysis

This project focuses on cleaning and analyzing EEG signals using Digital Signal Processing (DSP) techniques in MATLAB. The goal is to isolate critical brainwave frequencies (delta, theta, alpha, sigma, and beta) by implementing and evaluating various filters, enabling clearer interpretation of neurobehavioral traits.


## Repository Contents
- `Complete_Project_Cleaning EEG Data Using Filters`: MATLAB code for all filtering and analysis tasks.
- Dataset file: `825_2_PD_REST.mat`
- Project report:`Project_Report.pdf`


---

## Table of Contents
- [Introduction](#introduction)
- [Dataset Description](#dataset-description)
- [Implemented Techniques](#implemented-techniques)
- [Project Workflow](#project-workflow)
- [Results and Analysis](#results-and-analysis)
- [Conclusion and Future Work](#conclusion-and-future-work)
- [How to Run](#how-to-run)
- [References](#references)
- [Appendix](#appendix)

---

## Introduction
Electroencephalography (EEG) captures brain activity through electrodes placed on the scalp, offering insights into cognitive and emotional states. This project implements noise reduction techniques to clean EEG signals and filter specific brainwave bands for behavioral interpretation.

---

## Dataset Description
- **Source**: Multichannel EEG recordings
- **Structure**: Matrix of dimensions `[channels x samples]`
- **Sampling Rate**: `Fs` Hz, indicating temporal resolution
- **Key Components**:
  - Number of channels: `nbchannels`
  - Data: `EEG.data`
  - Sampling rate: `EEG.srate`

---

## Implemented Techniques
### Filters Applied
1. **Notch Filter**: Removes line noise at 60 Hz and harmonics.
2. **Frequency Band Filters**:
   - Delta (0.5–4 Hz)
   - Theta (4–7 Hz)
   - Alpha (8–12 Hz)
   - Sigma (12–16 Hz)
   - Beta (13–30 Hz)

### Signal Analysis
- Power Spectral Density (PSD) computation for noise identification and validation.
- Behavioral trait mapping for brain activity interpretation.

---

## Project Workflow
1. **Preliminary Analysis**:
   - Analyze raw EEG data for noise and artifacts.
   - Identify line noise peaks in PSD.
2. **Filtering**:
   - Apply notch filters for line noise removal.
   - Implement band-specific filters for brainwave extraction.
3. **Analysis and Visualization**:
   - Plot filtered EEG signals across time and frequency domains.
   - Interpret cognitive and emotional states based on frequency bands.

---

## Results and Analysis
1. **Noise Reduction**:
   - Significant attenuation of line noise observed in PSD.
2. **Filtered Data**:
   - Cleaned signals exhibit reduced artifacts, enabling better interpretation.
3. **Behavioral Insights**:
   - Delta: Associated with deep sleep.
   - Theta: Linked to relaxation and meditative states.
   - Alpha: Indicative of stress reduction and learning.
   - Sigma: Related to memory consolidation during sleep.
   - Beta: Reflects active thinking and concentration.

---

## Conclusion and Future Work
### Summary
- Effective noise reduction and frequency band isolation were achieved, enabling clearer understanding of brain activity.

### Future Directions
- **Adaptive Filtering**: Tailor filters to individual signal characteristics.
- **Machine Learning**: Use deep learning models for noise classification and artifact removal.
- **Real-Time Processing**: Implement dynamic filtering for live EEG analysis.

---

## How to Run
### Prerequisites
- MATLAB installed on your system.
- EEGLAB added to MATLAB's path.
    Details:
     
1.	Download: Visit the official EEG Lab website and navigate to the download section. Select the appropriate version of the toolbox compatible with our MATLAB version and operating system. https://sccn.ucsd.edu/eeglab/downloadtoolbox.php/
2.	Extract Files: After downloading the EEG Lab toolbox archive, extract its contents to a designated folder on our computer.
3.	Add Path: Open MATLAB and navigate to the "Set Path" option in the "File" menu. Add the path to the extracted EEG Lab toolbox folder to MATLAB's search path to enable access to its functions and scripts.
4.	Verify Installation: To verify the successful installation of the EEG Lab toolbox, type "eeglab" in the MATLAB command window. If installed correctly, EEG Lab's graphical user interface (GUI) should open without any errors.


### Steps
1. Clone this repository:
   ```bash
              git clone https://github.com/your-repo/EEG-Filtering-Analysis.git
  (https://github.com/ahtisham73/Advanced-Filtering-of-EEG-Signals-for-Enhanced-Neurobehavioral-Analysis.git)    
 

  
3. Add EEGLAB path in MATLAB:
   ```bash
   addpath('path_to_eeglab');
   `

4. Load the dataset:
     ```bash
    load('825_2_PD_REST.mat');



5. Run the MATLAB script:
     ```bash
      run('Complete_Project_Cleaning EEG Data Using Filters.m');




