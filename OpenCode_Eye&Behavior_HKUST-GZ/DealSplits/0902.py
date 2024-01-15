'''
derived from gaze
1 2 1
1.03 0.37 1.17
Data in this format
'''

from backupFunctions import *

# Get a list of video files in the 'etgvideos' directory with supported extensions
video_files = sorted([f for f in os.listdir('./etgvideos') if f.endswith(('.avi', '.mp4', '.flv', '.mkv'))])

try:
    v8model = modelinitV8() # Initialize the YOLOv8 model
except NameError:
    pass

# Define a function to compress a sequence of values
def compress_sequence(sequence):
    if not sequence:
        return [], []

    values = [sequence[0]]
    counts = [1]

    for num in sequence[1:]:
        if num == values[-1]:
            counts[-1] += 1
        else:
            values.append(num)
            counts.append(1)

    return values, counts

# Iterate through each video file
for video_file in video_files:
    video_path = os.path.join('./etgvideos', video_file)
    out_filename = os.path.splitext(video_file)[0] + ".txt"
    outfile_path = os.path.join('./fix', out_filename)
    etgfile_path = os.path.join('./etgdatas', out_filename)
    startendfile_path = os.path.join('./startend', out_filename)
    currentcap,currentfps,currentwidth,currentheight,currentfourcc = initvideo(video_path)

    # Read the start and end frames for annotation from a file
    with open(startendfile_path,'r') as f:
        line = f.readline().strip().split()
        start = int(line[0])
        end = int(line[1])

    etgdatas = deallabels(etgfile_path)
    datalist = []

    while currentcap.isOpened():
        CurrentframeCount = int(currentcap.get(cv2.CAP_PROP_POS_FRAMES))#从第0帧开始算的
        ret, srcframe = currentcap.read()

        if not ret:
            break
        boxes = detectframeV8(srcframe,v8model)
        value = 0
        for box in boxes:
            x1,y1,x2,y2 = box[0],box[1],box[2],box[3]# Extract coordinates of the bounding box
            item = str(box[4])# Get the class label of the detected object
            
            cv2.rectangle(srcframe,(x1,y1),(x2,y2),(0, 0, 255),1)
            cv2.putText(srcframe, item, (x1,y1), cv2.FONT_HERSHEY_SIMPLEX, 0.75, (0, 0, 255), 2)

            ydis = abs(y1-y2)
            x_center = int((x1 + x2)/2)
            y_center = int((y1 + y2)/2)
            item = str(box[4])
            if CurrentframeCount in etgdatas:# Iterate through the eye-tracking data for the current frame
                for x, y in etgdatas[CurrentframeCount]:
                    cv2.circle(srcframe, (x, y), 5, (255, 0, 0), -1)
                    dis = getDis(x,y,x_center,y_center)
                    if dis <= 1.1 * ydis:
                        value = 1

        # Check if the frame is within the specified start and end frames
        if start <= CurrentframeCount <= end:
            datalist.append(value)
        print(value)
        cv2.imshow("Video", srcframe)
        cv2.waitKey(2)


    # Compress the sequence of values and counts
    values, counts = compress_sequence(datalist)
    values = [(value + 1) for value in values]
    counts = [round(count / 30.0, 2) for count in counts]
    with open(outfile_path, "w") as f:
        # Write the 'values' list to the first line of the output file
        f.write(" ".join(map(str, values)) + "\n")
        # Write the 'counts' list to the second line of the output file
        f.write(" ".join(map(str, counts)) + "\n")
    currentcap.release()

    