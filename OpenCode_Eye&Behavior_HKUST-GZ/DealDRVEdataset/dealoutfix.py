'''
Merge close frames into an interval
'''

from backupFunctions import *
dir_list = sorted([f for f in os.listdir('.') if os.path.isdir(f) and f.isdigit()])

for dir_name in dir_list:
    outfixpath = os.path.join(dir_name, 'outfix.txt')
    outmirrorpath = os.path.join(dir_name, 'fixing.txt')
    segments = []

    with open(outfixpath, 'r') as file:
        lines = file.readlines()

    start_frame = None
    end_frame = None

    for line in lines:
        thisstart, thisend = map(int, line.strip().split())
        if start_frame is None:
            start_frame = thisstart
            end_frame = thisend
        else: #The previous segment is not over
            if thisstart - end_frame < 90:
                end_frame = thisend
            elif thisstart - end_frame >= 90:
                segments.append([start_frame, end_frame])
                start_frame = thisstart
                end_frame = thisend


    with open(outmirrorpath, "w") as f:
        for item in segments:
            line = ' '.join(str(i) for i in item) + '\n' 
            f.write(line) 