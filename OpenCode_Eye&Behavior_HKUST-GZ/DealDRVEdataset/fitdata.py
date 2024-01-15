'''
According to used.txt, delete the invalid interval in fixing.txt
'''

import os

with open('used.txt', 'r') as f:
    used_frames = [line.strip().split()[1] for line in f.readlines()]

dir_list = sorted([f for f in os.listdir('.') if os.path.isdir(f) and f.isdigit()])

for dir_name in dir_list:
    fixing_path = os.path.join(dir_name, 'fixing.txt')
    newfixing_path = os.path.join(dir_name, 'newfixing.txt')
    
    if os.path.exists(fixing_path):
        with open(fixing_path, 'r') as f:
            lines = f.readlines()
        
        # fit
        with open(newfixing_path, 'w') as f:
            for line in lines:
                start_frame = line.strip().split()[0]
                if start_frame in used_frames:
                    print(dir_name,start_frame)
                    f.write(line)
