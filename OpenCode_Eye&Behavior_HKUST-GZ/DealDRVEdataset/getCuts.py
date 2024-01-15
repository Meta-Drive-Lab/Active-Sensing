'''
According to the interception interval of newfixing.txt, intercept the video and etgdata
'''

from backupFunctions import *

dir_list = sorted([f for f in os.listdir('.') if os.path.isdir(f) and f.isdigit()])
classes = ['R', 'fan', 'armor', 'mid', 'left', 'right']


for dir_name in dir_list:
    video_path = os.path.join(dir_name, 'video_etg.avi')
    text_path = os.path.join(dir_name, 'etg_samples.txt')
    turnlabel_path = os.path.join(dir_name, 'newfixing.txt')

    currentcap,currentfps,currentwidth,currentheight,currentfourcc = initvideo(video_path)
    total_frames = int(currentcap.get(cv2.CAP_PROP_FRAME_COUNT))
    etgdatas = deallabels(text_path)
    turndatas = dealturnsplits(turnlabel_path)
    print(turndatas)

    while currentcap.isOpened():
        for turnrow in turndatas:
            start = turnrow[0] - 30
            end = turnrow[1] + 200
            print(turnrow)
            if start < 0:
                start = 0
            if end > total_frames:
                end = total_frames
            
            labelfilename = f"output_{int(start)}.txt"
            out_label_path = os.path.join(dir_name, labelfilename)
            filename = f"output_{int(start)}.avi"
            out_video_path = os.path.join(dir_name, filename)
            print("now is writing:",out_video_path)
            out = cv2.VideoWriter(out_video_path, currentfourcc, currentfps, (currentwidth, currentheight))
            frames_to_write = end - start + 1
            frames_written = 0
            currentcap.set(cv2.CAP_PROP_POS_FRAMES, start)
            with open(out_label_path, "w") as f:
                while frames_written < frames_to_write:
                    ret, frame = currentcap.read()
                    if not ret:
                        break
                    frame_num = int(currentcap.get(cv2.CAP_PROP_POS_FRAMES))
                    # frame = cv2.resize(frame, None, fx=0.5, fy=0.5, interpolation=cv2.INTER_NEAREST)
                    # cv2.imshow("srcFrame",frame)
                    # cv2.waitKey(1)
                    if frame_num in etgdatas:
                        for x, y in etgdatas[frame_num]:
                            f.write("{} {} {}\n".format(frame_num - start, x, y))
                    out.write(frame)
                    frames_written += 1
        currentcap.release()

    cv2.destroyAllWindows()

