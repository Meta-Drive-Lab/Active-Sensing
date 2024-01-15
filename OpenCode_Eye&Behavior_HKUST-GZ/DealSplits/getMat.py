'''
Integrate various previously generated txt files into mat format files
'''

import scipy.io as sio
import numpy as np
from backupFunctions import *

"""
fixitem 1:road 2:mirror
"""

"""
the item value cases
road
     3: The car in front is very close and cannot change lanes
     2: There is no car ahead and you can change lanes
     1: The car in front is far away and you can change lanes

mirror
     3: There are no cars on the side and you can change lanes
     2: There is a car on the side, it is far away, and you can change lanes
     1: There is a car on the side. It is very close and cannot change the lane.
"""
 
data_struct = np.empty((4,), dtype=[
    ('trialnum', float, (36,1)),
    ('fixdur', object, (36,)),
    ('fixitem', object, (36,)),
    ('itemval', float, (36,2)),
    ('choice', float, (36,1)),
    ('rt', float, (36,1)),
    ('tItem', float, (36,2))
])

outfiles = sorted([f for f in os.listdir('./out') if f.endswith(('.txt'))])
trialnum = 1
roadcnt = 0


for x in range(4):  # Loop through the 4 experiments
    for i in range(36):  # Each experiment has 36 trials
        outfile = outfiles.pop(0)  # Get and remove the first element from the sorted list
        outfilepath = os.path.join('./out', outfile)
        fixfilepath = os.path.join('./fix', outfile)

        with open(fixfilepath, "r") as f:
            lines = f.readlines()
        fixitem = list(map(int, lines[0].strip().split()))
        fixdur = list(map(float, lines[1].strip().split()))

        # if len(fixdur) > 8:
        #     trialnum, fixdur, fixitem, itemval, choice, rt, tItem = lasttrialnum, lastfixdur, lastfixitem, lastitemval, lastchoice, lastrt, lasttItem

        with open(outfilepath, "r") as f:
            lines = f.readlines()
        choicedata = lines[1].strip()
        choice = 1 if choicedata == 'no' else 2
        
        road_value = int(lines[2].strip().split(": ")[1])
        mirror_value = int(lines[3].strip().split(": ")[1])
        itemval = [road_value, mirror_value]
        
        rt = round(sum(fixdur) + 1.5, 2)
        tItem = [round(sum(fixdur[i] for i, item in enumerate(fixitem) if item == 1), 2),
                 round(sum(fixdur[i] for i, item in enumerate(fixitem) if item == 2), 2)]

        lasttrialnum, lastfixdur, lastfixitem, lastitemval, lastchoice, lastrt, lasttItem = trialnum, fixdur, fixitem, itemval, choice, rt, tItem
        
        # Store the values in the numpy structured array
        data_struct['trialnum'][x][i] = trialnum
        data_struct['fixdur'][x][i] = fixdur
        data_struct['fixitem'][x][i] = fixitem
        data_struct['itemval'][x][i] = itemval
        data_struct['choice'][x][i] = choice
        data_struct['rt'][x][i] = rt
        data_struct['tItem'][x][i] = tItem
        
        trialnum += 1


# 保存为MAT文件
sio.savemat('datastruct_hkust06281.mat', {'dstruct_real': data_struct}, appendmat=False)