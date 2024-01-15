'''
Based on the gaze situation and the rearview mirror position predicted by the model, 
the situation of all frames looking at the rearview mirror is obtained.
'''


from backupFunctions import *

# Get a list of directory paths with numeric names
dir_list = sorted([f for f in os.listdir('.') if os.path.isdir(f) and f.isdigit()])

try:
    v8model = modelinitV8() # Initialize the YOLOv8 model
except NameError:
    pass

for dir_name in dir_list:
    # Get the paths for video and text data files
    video_path = os.path.join(dir_name, 'video_etg.avi')
    text_path = os.path.join(dir_name, 'etg_samples.txt')
    outfixpath = os.path.join(dir_name, 'outfix.txt')

    # Initialize video capture and other properties
    currentcap,currentfps,currentwidth,currentheight,currentfourcc = initvideo(video_path)

    # Load eye-tracking data
    etgdatas = deallabels(text_path)
    is_paused = False
    datalist = []

    while currentcap.isOpened():
        CurrentframeCount = int(currentcap.get(cv2.CAP_PROP_POS_FRAMES))#从第0帧开始算的
        if not is_paused:
            ret, srcframe = currentcap.read()
            if not ret:
                break

            # Detect objects in the frame using YOLOv8 model
            boxes = detectframeV8(srcframe,v8model)
            if CurrentframeCount in etgdatas:
                for x, y in etgdatas[CurrentframeCount]:
                    cv2.circle(srcframe, (x, y), 5, (255, 0, 0), -1)
            value = 0

            # Iterate through each detected object in the 'boxes' list
            for box in boxes:
                x1,y1,x2,y2 = box[0],box[1],box[2],box[3]
                item = str(box[4])
                cv2.rectangle(srcframe,(x1,y1),(x2,y2),(0, 0, 255),1)
                cv2.putText(srcframe, item, (x1,y1), cv2.FONT_HERSHEY_SIMPLEX, 0.75, (0, 0, 255), 2)
                ydis = abs(y1-y2)
                x_center = int((x1 + x2)/2)
                y_center = int((y1 + y2)/2)
                item = str(box[4])

                # Check if the current frame count is in the 'etgdatas' dictionary
                if CurrentframeCount in etgdatas:
                    for x, y in etgdatas[CurrentframeCount]:
                        cv2.circle(srcframe, (x, y), 5, (255, 0, 0), -1)
                        dis = getDis(x,y,x_center,y_center)
                        # If the distance is within 1.2 times the height of the bounding box, set 'value' to 1
                        if dis <= 1.2 * ydis:
                            value = 1
                            
            print(CurrentframeCount,value)
            currentdata = [CurrentframeCount,value]
            datalist.append(currentdata)
            srcframe = cv2.resize(srcframe, None, fx=0.5, fy=0.5, interpolation=cv2.INTER_NEAREST)
            cv2.imshow("srcFrame",srcframe)

    with open(outfixpath, "w") as f:
        for item in datalist:
            line = ' '.join(str(i) for i in item) + '\n'  
            f.write(line)  

    currentcap.release()
    cv2.destroyAllWindows()