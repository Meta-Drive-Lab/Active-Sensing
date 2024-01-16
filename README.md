# Active-Sensing

This repository contains the toolkit of data processing and MATLAB scripts to reproduce the main and extended figures of the following article: Capturing Human Active Sensing in Real-world Driving Tasks with Behavioural Causality

# Data Processing Guide

This guide covers the steps for processing the DR(eye)VE open-source dataset, including code and workflow. Please note that the BIT dataset is a reference dataset shared by BIT and is temporarily not available. Therefore, this part focuses on the processing methods for the DR(eye)VE dataset. The data processing steps for BIT are not within the scope of this guide. This part is dedicated to simple data integration and processing.

Our goal is to convert the DR(eye)VE dataset into MAT format for further analysis and research. Below is an overview of the entire processing workflow, including the required environment setup and detailed instructions for each processing step. The relevant files are in the [`OpenCode_Eye&Behavior_HKUST-GZ`] folder.

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

To execute the Matlab scripts, you need to download and extract them into a folder of your preference, and then navigate to this folder using MATLAB. The relevant files are in the [`FigurePlot_HKUST-GZ`] folder.

The [`data`] folder contains the human driving behavioural data.
The [`functions`] folder contains all the custom MATLAB functions utilized in the main scripts. This directory will be included at the start of each script.

All scripts have been tested and run under MATLAB R2023a on a Windows10 system. There are four scripts responsible for generating Figure2-5 in the results section of the article, and FigureA2 in the extended section. These four scripts are designed to utilize the functions provided within this package to generate all the figures. By simply setting the rootdir variable, which is declared at the beginning of each script, to the directory containing the scripts, they should execute correctly.
1. [`Lanechange2023_rdv.m`]
2. [`Lanechange2023_times.m`]
3. [`Lanechange2023_behaviour.m`]
4. [`Lanechange2023_extended.m`]

# License

Copyright (c) 2024, Hongliang LU & Yunmeng LIU

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


