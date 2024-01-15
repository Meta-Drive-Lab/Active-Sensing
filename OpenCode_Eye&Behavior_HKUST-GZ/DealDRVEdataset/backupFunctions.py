try:
    import torch
except ImportError:
    pass  # If the import fails, continue without raising an error
try:
    from ultralytics import YOLO
except ImportError:
    pass
import numpy as np
import cv2
from time import time
import sys
import math
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import shutil


def deallabels(labelfilename):
    # Function to process label data from a file and return it as a dictionary
    data = {}
    with open(labelfilename, "r") as f:
        for line in f:
            parts = line.strip().split()
            etgdata = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
            if etgdata not in data:
                data[etgdata] = []
            if not math.isnan(x) and not math.isnan(y):
                data[etgdata].append((int(x), int(y)))
    return data
'''
Usage of deallabels:
if CurrentframeCount in currentdata:
    for x, y in currentdata[CurrentframeCount]:
        cv2.circle(srcframe, (x, y), 5, color, -1)
'''


'''
We initially used yolov5 and later used the more advanced yolov8. 
In this project, we did not make any changes to the official code of yolo.
'''
def modelinit():
    # Initialize the YOLO model (YOLOv5 in this case)
    model = torch.hub.load('ultralytics/yolov5', 'custom', 'best.pt')
    device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
    model.to(device)
    return model

def modelinitV8():
    # Initialize the YOLOv8 model
    model = YOLO('bestv8.pt')
    return model

def detectframeV8(image,model):
    # Detect objects in an image using YOLOv8 model
    results = model(image)
    # Process results list
    datalist = []
    for result in results:
        data = result.boxes.cpu().numpy()
        for box in data:
            x1,y1,x2,y2 = box.xyxy[0]
            cls = box.cls
            cls = int(cls[0])            
            thisdata = [int(x1),int(y1),int(x2),int(y2),cls]
            datalist.append(thisdata)
    return datalist

'''
Usage:
for boxes in datalist:
    x1,y1,x2,y2 = boxes[0],boxes[1],boxes[2],boxes[3]
    item = boxes[4]
'''

'''
Usage:
        m_labels, m_cord = detectframe(frame, v5model)
        n = len(m_labels)
        x_shape, y_shape = frame.shape[1], frame.shape[0]

        for i in range(n):
            row = m_cord[i]
            item = classes[int(m_labels[i])]
            print(item)
            x1, y1, x2, y2 = int(row[0] * x_shape), int(row[1] * y_shape), int(row[2] * x_shape), int(row[3] * y_shape)
            bgr = (0, 0, 255)
'''

def detectframe(image, model):
    # Detect objects in an image using YOLOv5 model
    results = model([image])
    labels, cord = results.xyxyn[0][:, -1].to('cpu').numpy(), results.xyxyn[0][:, :-1].to('cpu').numpy()
    return labels, cord

def getDis(x1,y1,x2,y2):
    # Calculate the distance between two points
    xdis = abs(x1-x2)
    ydis = abs(y1-y2)
    dis = int(math.sqrt(xdis*xdis + ydis*ydis))
    return dis


def initvideo(videofilename):
    # Initialize video capture from a file and retrieve video properties
    cap = cv2.VideoCapture(videofilename)
    fps = int(cap.get(cv2.CAP_PROP_FPS))
    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    return cap,fps,width,height,fourcc

def getSpeedData(SpeedData_filename):
    # Read and process speed data from a file
    data = []
    with open(SpeedData_filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            gardata = int(parts[0])
            speed = float(parts[1]) * 1000 / 3600
            course = (90 - float(parts[2])) * np.pi / 180
            if not math.isnan(speed) and not math.isnan(course):
                data.append((gardata,speed, course))
    return data

def dealturnsplits(turn_data_filename):
    # Read and process turn/split data from a file
    data = []
    with open(turn_data_filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            start = int(parts[0])
            end = int(parts[1])
            if not math.isnan(start) and not math.isnan(end):
                data.append((start,end))
    return data

def processSpeedData(data):
    # Process speed data to calculate points and positions
    points = []
    prev_point = np.array([0, 0])
    for gardata,speed, course in data:
        # calculate points and positions
        point = prev_point + 0.04 * np.array([speed * np.cos(course), speed * np.sin(course)])
        points.append((gardata,point))
        prev_point = point
    return points

def getFigdata(points):
    # Create a figure for plotting the route
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot([p[1][0] for p in points], [p[1][1] for p in points], 'r-', linewidth=1)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True)
    return fig,ax

def drawRoadmap(framecount,points,fig,ax,Smallframe,circle,change):
    # Draw the roadmap on the figure
    if change:
        patches = ax.patches
        circles = [patch for patch in patches if isinstance(patch, Circle)]
        for patch in circles:
            patch.remove()
        circle.center = (points[framecount][1][0], points[framecount][1][1])
        ax.add_patch(circle)
        if Smallframe:
            ax.set_xlim([points[framecount][1][0]-100, points[framecount][1][0]+100])
            ax.set_ylim([points[framecount][1][1]-100, points[framecount][1][1]+100])
    fig.canvas.draw()
    frame = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8).reshape(fig.canvas.get_width_height()[::-1] + (3,))
    frame = cv2.cvtColor(frame, cv2.COLOR_RGB2BGR)
    return frame

def waitkeyX():
    # Wait for a key press event
    press = cv2.waitKey(1) & 0xFF
    return press


