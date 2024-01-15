# Active-Sensing

This repository contains the toolkit of data processing and MATLAB scripts to reproduce the main and extended figures of the following article: Capturing Human Active Sensing in Real-world Driving Tasks with Behavioural Causality

# Data Processing Guide

This guide covers the steps for processing the DR(eye)VE open-source dataset, including code and workflow. Please note that the BIT dataset is a reference dataset shared by BIT and is temporarily not available. Therefore, this part focuses on the processing methods for the DR(eye)VE dataset. The data processing steps for BIT are not within the scope of this guide. This part is dedicated to simple data integration and processing.

Our goal is to convert the DR(eye)VE dataset into MAT format for further analysis and research. Below is an overview of the entire processing workflow, including the required environment setup and detailed instructions for each processing step.

## 0. Environment Configuration

We use libraries such as opencv-python and pytorch for data processing. Use the following command to set up the required environment:

```shell
pip install -r requirements.txt
```

## 1. Processing the DR(eye)VE Dataset

Each small folder in this dataset contains multiple video files and text files with gaze data. Our first step is to extract relevant segments where the driver is looking at the rearview mirror.

### 1.1 Obtaining Gaze Data

Use the following command to obtain gaze data:

```shell
python getoutfix.py
```

This will generate a file named `outfix.txt`, where the first column represents the frame number, and the second column indicates whether the driver is looking at the rearview mirror.

### 1.2 Integrating Gaze Data

Use the following command to integrate gaze data, combining consecutive gaze frame numbers into continuous gaze segments:

```shell
python dealoutfix.py
```

This will generate a file named `fixing.txt`, with the first column indicating the starting frame and the second column indicating the ending frame of each gaze segment.

### 1.3 Further Segment Integration

Use the following command to merge gaze segments that are too close together:

```shell
python fitdata.py
```

### 1.4 Extracting Video Segments

Use the following command to extract video segments from the `video_etg.avi` file based on the gaze segments and copy them to the `Splits` folder:

```shell
python getCuts.py
```

## 2. Processing the Splits Folder

In the `Splits` folder, you can use the following code for processing:

```shell
python player.py
```

This is a simple player code that can be used to work with the previously extracted video and text files. The video display shows the current gaze location.

- You can fine-tune the extracted segments using the slider, frame numbers, and video display.
- The "Save" button is used to save the `itemvalue` values for "road" and "mirror." You can refer to the code for value assignment rules.
- The "Yes" and "No" buttons are for saving the final choices.
- Other buttons allow for video playback, switching, and deleting video segments.

## 3. Processing Gaze Data

Use the following command to process video segments and obtain gaze segments in the required format:

```shell
python gazedata.py
```

## 4. Generating MAT Files

Use the following command to consolidate all the previously obtained data into a `.mat` file, which will be used by subsequent model code for further processing:

```shell
python getMat.py
```

This guide provides a comprehensive overview of the steps involved in processing the DR(eye)VE dataset and converting it into MAT format. If you encounter any issues or require more detailed information, please feel free to reach out.

# Figure Generating Guide

Description of data structure fields in the generated .mat file:
- trialnum: trial number 
- fixdur: fixation duration [s] per each item fixation
- fixitem: fixated item (item1-RV or item2-RV)
- itemval: perceptual states value associated with each item 
- choice: decision at end of trial (1-Lane Changing or 2-Lane Keeping)
- rt: response time [s]
- tItem: total fixation time [s] spent on either item

## Installation

To run the Matlab scripts, download and extract them into a folder of your choice, and navigate to this folder within MATLAB.

The [`functions`] folder contains all the custom MATLAB functions that are used in the main scripts. This directory will be added in the beginning of each script.

All scripts have been tested and run under MATLAB R2023a in Windows10 system. There are four scripts that plot Figure2, Figure3, Figure4 in the result section of the article and FigureA2 in the extended section. These four scripts should utilize the functions included in this package to plot all the figures. Simply set the `rootdir` variable, that is declared at the beginning of each script, to the directory that includes the scripts, and they should run correctly.

1. [`Lanechange2023_fig2_rdv.m`]
2. [`Lanechange2023_fig3_boxplot.m`]
3. [`Lanechange2023_fig4_times.m`]
4. [`Lanechange2023_figextend_behavior.m`]

