from backupFunctions import *
'''
Consolidate a series of image files into a video file
'''


dir_list = sorted([f for f in os.listdir('.') if os.path.isdir(f) and f.isdigit()])

output_garvideo = 'etg_video.mp4'
fps = 30

for dir_name in dir_list:
    camframesfolder_path = os.path.join(dir_name,'eyeframe')
    output_video = os.path.join(dir_name,output_garvideo)
    camframes = sorted([f for f in os.listdir(camframesfolder_path) if f.endswith('.jpg')])

    tempimg = cv2.imread(os.path.join(camframesfolder_path, camframes[0]))
    height, width, layers = tempimg.shape
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(output_video, fourcc, fps, (width, height))

    for image_file in camframes:
        image_path = os.path.join(camframesfolder_path, image_file)
        img = cv2.imread(image_path)
        video.write(img)
        print(f'writing {image_file}')

    video.release()
    cv2.destroyAllWindows()
    print(f"Video generated: {output_video}")
