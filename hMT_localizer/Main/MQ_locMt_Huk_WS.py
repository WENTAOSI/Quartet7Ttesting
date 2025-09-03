"""
Motion Localiser according to Huk

Sequence used: TR=1s
In the MotionQuartet experiment this localizer runs at the beginning of
the scanning session for 4 min 50s [290 volumes]. 
Block design: moving dots (10s), static dots (10s). [14 repetitions]
7 target are showed.

If TR=2s 
period=10 yields 5 slices per period; 4 min 50s with 145 volumes
period=20 yield 10 slices per period but double the scan time (9 min 40s) with 290 volumes


Psychopy3 (v2020.2.4)
Based on https://github.com/MSchnei/motion_quartet_scripts (@author: Marian.Schneider)

@author: siwentao based on @author: Alessandra Pizzuti adopted from @author Marian.Schneider
"""
from psychopy import visual, event, core,  monitors, logging, gui, data, misc
from psychopy.tools.coordinatetools import pol2cart, cart2pol
import numpy as np
import os
import csv
import sys
import time

###############################################################
TR=2 # 2 sec TR
TRIGGER_KEY='quoteleft' # psychopy accepts "quoteleft" as "`" 
period = 16 #sec of each period 

# set global offset
ho_dva = 0
vo_dva = 0
global_offset = (ho_dva, -vo_dva)
 
def apply_global_offset(base_pos=(0,0), global_offset=global_offset):
    return (base_pos[0] + global_offset[0], base_pos[1] + global_offset[1])
'''
          ↑ +Y
          |
 (-,+ )   |   (+,+)
          |
←--------(0,0)--------→ +X
          |
 (-,-)    |   (+,-)
          |
          ↓ -Y
'''
###############################################################

#%% GENERAL PARAMETERS """
targetDur = 0.3  # target lenght in s
# set properties for moving dots
# The number of dots
nDots = 200
# specify speed in units per frame
dotSpeed = 8  # deg per s [8]
# dot Life, how long should a dot life
dotLife = 10  # number of frames [10]
# The size of the dots [diameter]
dotSize = 0.2  # in deg
# misc.deg2pix(0.2, myWin.monitor)
# The size of the field.
FieldSizeRadius = 4  # radius in deg
# How close can dots come to center of screen (i.e. fixation cross)
innerBorder = 0.5  # distance from screen center in deg

# set background color
backgrColor = [-0.5, -0.5, -0.5]
# set dot color
dotColor = np.multiply(backgrColor, -1)  # from -1 (black) to 1 (white)

#%%""" SAVING and LOGGING """
# Store info about experiment and experimental run
expName = 'locMt_MQ'  # set experiment name here
expInfo = {
    'run': '1',
    'participant': 'test',
    'display': ['Vanderbilt7T', 'dbic']
    }

# Create GUI at the beginning of exp to get more expInfo
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False: core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName

# get current path and save to variable _thisDir
_thisDir = os.path.dirname(os.path.abspath(__file__))
parentDir = os.path.dirname(_thisDir)
os.chdir(parentDir)  # change directory to this pat

# get parent path and move up one directory
str_path_parent_up = os.path.abspath(
    os.path.join(os.path.dirname( __file__ ), '..','..'))

# Name and create specific subject folder
subjFolderName = str_path_parent_up + os.path.sep + '%s_SubjData' \
    % expInfo['participant'] + os.path.sep + 'hMT+localizer'

if not os.path.isdir(subjFolderName):
    os.makedirs(subjFolderName)

# Name and create specific folder for logging results
logFolderName = subjFolderName + os.path.sep + 'Logging'
if not os.path.isdir(logFolderName):
    os.makedirs(logFolderName)
logFileName = logFolderName + os.path.sep + '%s_%s_Run%s_%s' % (
    expInfo['participant'], expInfo['expName'], expInfo['run'],
    expInfo['date'])

# Name and create specific folder for pickle output
outFolderName = subjFolderName + os.path.sep + 'Output'
if not os.path.isdir(outFolderName):
    os.makedirs(outFolderName)
outFileName = outFolderName + os.path.sep + '%s_%s_Run%s_%s' % (
    expInfo['participant'], expInfo['expName'], expInfo['run'],
    expInfo['date'])

# Name and create specific folder for BV protocol files
prtFolderName = subjFolderName + os.path.sep + 'Protocols'
if not os.path.isdir(prtFolderName):
    os.makedirs(prtFolderName)

# save a log file and set level for msg to be received
logFile = logging.LogFile(logFileName+'.log', level=logging.INFO)
logging.console.setLevel(logging.WARNING)  # console receives warnings/errors

# create array to log key pressed events
TimeKeyPressedArray = np.array([], dtype=float)
print(logFileName)
print(expInfo['run'])

#%%"""MONITOR AND WINDOW"""
# source:https://www.dartmouth.edu/dbic/research_infrastructure/peripherals.html
# set monitor information:
if expInfo['display'] == 'dbic':
    distanceMon = 128.7  # cm
    widthMon = 42.8  # cm
    PixW = 1920  # cm
    PixH = 1080 # cm

elif expInfo['display'] == 'Vanderbilt7T':
    distanceMon = 48  # cm
    widthMon = 17  # cm
    PixW = 1024  # cm
    PixH = 768 # cm
    
moni = monitors.Monitor('testMonitor', width=widthMon, distance=distanceMon)
moni.setSizePix([PixW, PixH])

# log monitor info
logFile.write('MonitorDistance=' + str(distanceMon) + 'cm' + '\n')
logFile.write('MonitorWidth=' + str(widthMon) + 'cm' + '\n')
logFile.write('PixelWidth=' + str(PixW) + '\n')
logFile.write('PixelHeight=' + str(PixH) + '\n')

# set screen:
if expInfo['display'] == 'Vanderbilt7T':
    screen=1
elif expInfo['display'] == 'dbic':
    screen=0

myWin = visual.Window(
    size=(PixW, PixH),
    screen=screen,
    winType='pyglet',  # winType : None, ‘pyglet’, ‘pygame’
    allowGUI=False,
    allowStencil=True,
    fullscr=True,  # for psychoph lab: fullscr = True
    monitor=moni,
    color=backgrColor,
    colorSpace='rgb',
    units='deg',
    blendMode='avg'
    )

#%%"""CONDITIONS AND DURATIONS"""
# Set parameters
NumOfCond = 1  # Number of different conditions [ originally there were 3: Right, Center and Left; currently Center is used]
NumRepCond = 14 # Number of repetitions per condition per run [27]
NrOfTargets = 7  # Number of targets that participants need to detect [20]

StimDur = period/TR  # Dur moving dots, in TRs
IsiDur = np.array([period/TR])
FixDur = period/TR  # Dur fixation in beginning and end

# Define conditions
conditions = []
for ind in range(0, NumRepCond):
    CondOrder = np.arange(1, NumOfCond+1)
    np.random.shuffle(CondOrder)
    BaseOrder = np.zeros(NumOfCond, dtype=int)
    block_elem = np.insert(BaseOrder, np.arange(len(CondOrder)), CondOrder)
    conditions = np.hstack((conditions, block_elem))

print(f"CONDITIONS: {conditions}")
# add -1 in beginning and end, for fixation dot
conditions = np.hstack(([-1], conditions))
# make sure array is int
conditions = conditions.astype(int)

# Define durations of stimulus and rest
durations = np.ones(len(conditions), dtype=int)*StimDur
# tile IsiDur so it matches the NumRepCond
if NumRepCond % len(IsiDur) != 0:
    print('WARNING: NumRepCond not exact multiples of IsiDur')
    IsiDur = np.tile(IsiDur, int(NumRepCond/len(IsiDur)))

for ind in CondOrder:
    Pos = np.where(conditions == ind)[0]+1
    np.random.shuffle(IsiDur)
    durations[Pos] = IsiDur
Pos = np.where(conditions == -1)
durations[Pos] = FixDur  # Dur fixation

conditions[Pos] = 0 # for prt saving and stim presentation

# Define random target
lgcRep = True
while lgcRep:

    targets = np.random.choice(np.arange(FixDur, np.sum(durations)-FixDur),
                               NrOfTargets, replace=False)
    # check that two Targets do not follow each other immediately
    lgcRep = np.greater(np.sum(np.diff(np.sort(targets)) <= 1), 0)

targets = np.sort(targets)

# Calculate timestamps (cumulative time at the end of each condition)
timestamps = np.cumsum(durations)

# Calculate Onset as timestamps - durations
onsets = np.zeros(len(timestamps), dtype=int)
for ind in range(0, len(timestamps)-1):
    onsets[ind] = timestamps[ind] - durations[ind]


# group condition info 
condition_info = {
'Conditions': conditions,
'Durations': durations,
'TimeStamps': timestamps,
'Targets': targets
}

#filename_npz = os.path.join(prtFolderName, f'MtLoc_MQ_4min_4_run{expInfo["run"]}')
# save the conditions as .npy into protocal
#np.savez(filename_npz, **condition_info) #save as filename_npz.npz

# Save the conditions as .csv
filename_csv = os.path.join(prtFolderName, f'MtLoc_MQ_4min_4_run{expInfo["run"]}.csv')

# Save .csv for inspection and record
with open(filename_csv, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write headers
    writer.writerow(["Index", "Condition", "Duration", "Timestamp", "Onset", "Target"])
    
    # Iterate through rows and check target mapping
    for idx, (cond, dur, cum_duration, onset) in enumerate(zip(conditions, durations, timestamps, onsets)):
        # Determine the target for the current row
        if idx == 0:  # First row has no previous cumulative duration
            prev_cum_duration = 0
        else:
            prev_cum_duration = timestamps[idx - 1]

        # Check if a target falls in this range
        target_in_range = [target for target in targets if prev_cum_duration < target <= cum_duration]
        target_value = target_in_range[0] if target_in_range else "None"
        
        # Write row
        writer.writerow([idx, cond, dur, cum_duration, onset, target_value])
print(f"CSV file saved at {filename_csv}")

logFile.write('Conditions=' + str(conditions) + '\n')
# load durations of stimulus and rest
durations = condition_info['Durations']
logFile.write('Durations=' + str(durations) + '\n')
# load the target onsets
targets = condition_info['Targets']
logFile.write('Targets=' + str(targets) + '\n')

#%%"""TIME, TIMING AND CLOCKS"""
# parameters
totalTime = np.sum(durations)
logFile.write('TotalTime=' + str(totalTime) + '\n')

# give system time to settle before it checks screen refresh rate
core.wait(1)

# get screen refresh rate
refr_rate = myWin.getActualFrameRate()

if refr_rate is not None:
    print("refresh rate: %i" % refr_rate)
    frameDur = 1.0 / round(refr_rate)
    print("actual frame dur: %f" % frameDur)
else:
    # couldnt get reliable measure, guess
    refr_rate = 60.0
    frameDur = 1.0 / 60.0
    print("fake frame dur: %f" % frameDur)

logFile.write('RefreshRate=' + str(refr_rate) + '\n')
logFile.write('FrameDuration=' + str(frameDur) + '\n')

# define clock
clock = core.Clock()
logging.setDefaultClock(clock)

#%%"""STIMULI"""
# divide the dotSpeed by the refresh rate
dotSpeed = dotSpeed / refr_rate

logFile.write('nDots=' + str(nDots) + '\n')
logFile.write('speed=' + str(dotSpeed) + '\n')
logFile.write('dotSize=' + str(dotSize) + '\n')
logFile.write('fieldSizeRadius=' + str(FieldSizeRadius) + '\n')

# initialise moving dot stimuli
dotPatch = visual.ElementArrayStim(
    myWin,
    fieldPos=apply_global_offset((0.0, 0.0)),
    autoLog=False,
    elementTex=None,
    name='dotPatch',
    fieldShape='circle',
    elementMask='circle',
    nElements=int(nDots),
    sizes=dotSize,
    units='deg',
    fieldSize=FieldSizeRadius*2,
    colors=dotColor
    )

# fixation dot
dotFix = visual.Circle(
    myWin,
    autoLog=False,
    name='dotFix',
    units='deg',
    radius=0.05,
    pos=apply_global_offset((0, 0)),
    fillColor=[1.0, 0.0, 0.0],
    lineColor=[1.0, 0.0, 0.0],
    )

dotFixSurround = visual.Circle(
    myWin,
    autoLog=False,
    name='dotFix',
    units='deg',
    radius=0.1,
    pos=apply_global_offset((0, 0)),
    fillColor=[0.5, 0.5, 0.0],
    lineColor=[0.0, 0.0, 0.0],
    )

# control text
controlText = visual.TextStim(
    win=myWin,
    colorSpace='rgb',
    color=[1.0, 1.0, 1.0],
    height=0.5,
    pos=apply_global_offset((0.0, -4.0)),
    autoLog=False,
    )

# text at the beginning of the experiment
triggerText = visual.TextStim(
    win=myWin,
    colorSpace='rgb',
    color=[1.0, 1.0, 1.0],
    height=0.5,
    pos=apply_global_offset((0, 0)),
    text='Experiment will start soon. Waiting for scanner'
    )

instruct_text = visual.TextStim(
    win=myWin,
    colorSpace='rgb',
    color=[1.0, 1.0, 1.0],
    height=0.5,
    pos=apply_global_offset((0, 0)),
    text='Please Fixate at the CENTER \n\
        Press 1 Immediately If Fixation Changes Color \n\
             Now Press 1 to CONTINUE'
    )
#%%"""FUNCTIONS"""
# function to determine initial dot positions
def dots_init(nDots):
    # specify the angle for each dot
    dotsTheta = np.random.rand(nDots)*360
    # specify the distance to the centre
    dotsRadius = (np.random.rand(nDots)**0.5)*FieldSizeRadius
    # convert
    dotsX, dotsY = pol2cart(dotsTheta, dotsRadius)
    # create array frameCount
    frameCount = np.random.uniform(0, dotLife, size=len(dotsX)).astype(int)
    return dotsX, dotsY, frameCount

def dots_update(dotsX, dotsY, frameCount, dotSpeed=dotSpeed,
                frameDeathAfter=np.inf):
    # convert to polar coordinates
    dotsTheta, dotsRadius = cart2pol(dotsX, dotsY)
    # update radius
    dotsRadius = (dotsRadius+dotSpeed)
    # decide which dots die
    lgcOutFieldDots = np.zeros(len(dotsTheta), dtype='bool')
    if dotSpeed > 0:
        # create lgc for elems where radius too large
        lgcOutFieldDots = (dotsRadius >= FieldSizeRadius)
    elif dotSpeed < 0:
        # create lgc for elems where radius too small
        lgcOutFieldDots = (dotsRadius <= innerBorder)
    # create logical for where frameCount too high
    lgcFrameDeath = (frameCount >= frameDeathAfter)
    # combine logicals
    lgcDeath = np.logical_or(lgcOutFieldDots, lgcFrameDeath)
    # replace dead dots
    dotsRadius[lgcDeath] = np.random.uniform(innerBorder, FieldSizeRadius,
                                             size=sum(lgcDeath))
    dotsTheta[lgcDeath] = np.random.rand(sum(lgcDeath))*360
    # convert
    dotsX, dotsY = pol2cart(dotsTheta, dotsRadius)
    # increase frameCount for every elements
    frameCount += 1
    # set the counter for newborn dots to zero
    frameCount[lgcDeath] = 0
    return dotsX, dotsY, frameCount

# target function
nrOfTargetFrames = int(targetDur/frameDur)
print("number of target frames")
print(nrOfTargetFrames)


#%%"""RENDER_LOOP"""
# Create Counters
i = 0  # counter for blocks
nTR= 1; # total TR counter
# draw dots for the first time [inward dots]
dotsX, dotsY, frameCntsIn = dots_init(nDots)
# set x and y positions to initialized values
dotPatch.setXYs(np.array([dotsX, dotsY]).transpose())

# define DirTime and calculate FlickerCycle
DirTime = 1  # move in one dir before moving in opposite
AxisTime = DirTime*2  # because we have 2 dir states (back and forth)

# give system time to settle before stimulus presentation
core.wait(1.0)

# Present Instruction
instruct_text.draw()
myWin.flip()
event.waitKeys(keyList=['1'], timeStamped=False)

# wait for scanner trigger
triggerText.draw()
myWin.flip()
event.waitKeys(keyList=[TRIGGER_KEY], timeStamped=False)
# reset clock
clock.reset()
logging.data('StartOfRun' + str(expInfo['run']))

runExp = True
once = False
mtargetCounter = 0
counter = 0
totalTRs = sum(durations)
print('TRs ' + str(totalTRs))
print('Targets ' + str(targets))

logging.data('TR ' + str(nTR))



while runExp:

    if nTR < totalTRs :

        # low-level rest (only central fixation dot)
        if conditions[i] == -1:
            # set loopDotSpeed to zero
            loopDotSpeed = 0
            # set loopDotLife to inf
            loopDotLife = np.inf
            # set opacities
            dotPatch.opacities = 0
            logging.data('CondFix' + '\n')

        # static dots rest
        elif conditions[i] == 0:
            # set opacity to 1 for all static
            dotPatch.opacities = 1
            # set loopDotSpeed to zero
            loopDotSpeed = 0
            # set loopDotLife to inf
            loopDotLife = np.inf
            # find out what sort of rest
            if conditions[i-1] == 1:  # central static
                conditions[i] = 4
            elif conditions[i-1] == 2:  # left static
                conditions[i] = 5
            elif conditions[i-1] == 3:  # right static
                conditions[i] = 6
            logging.data('CondStat' + '\n')
        # central motion
        elif conditions[i] == 1:
            # set loopDotSpeed to dotSpeed
            loopDotSpeed = dotSpeed
            # set loopDotLife to dotLife
            loopDotLife = dotLife
            # set opacaities
            dotPatch.opacities = 1
            dotPatch.fieldPos = apply_global_offset((0.0, 0.0))
            logging.data('CondCenter' + '\n')
        # left motion
        elif conditions[i] == 2:
            # set loopDotSpeed to dotSpeed
            loopDotSpeed = dotSpeed
            # set loopDotLife to dotLife
            loopDotLife = dotLife
            # set opacaities
            dotPatch.opacities = 1
            dotPatch.fieldPos = [(-5.0-FieldSizeRadius), 0.0]
            logging.data('CondLeft' + '\n')
        # right motion
        elif conditions[i] == 3:
            # set loopDotSpeed to dotSpeed
            loopDotSpeed = dotSpeed
            # set loopDotLife to dotLife
            loopDotLife = dotLife
            # set opacaities
            dotPatch.opacities = 1
            dotPatch.fieldPos = [(5.0+FieldSizeRadius), 0.0]
            logging.data('CondRight' + '\n')

        while nTR < np.sum(durations[0:i+1]): #clock.getTime() < np.sum(durations[0:i+1]):
            # update dots
            t = clock.getTime()

            if t % AxisTime < DirTime:
                dotsX, dotsY, frameCntsIn = dots_update(
                    dotsX, dotsY, frameCntsIn, dotSpeed=loopDotSpeed,
                    frameDeathAfter=loopDotLife)

            elif t % AxisTime >= DirTime and t % AxisTime < 2*DirTime:
                dotsX, dotsY, frameCntsIn = dots_update(
                    dotsX, dotsY, frameCntsIn, dotSpeed=-loopDotSpeed,
                    frameDeathAfter=loopDotLife)

            dotPatch.setXYs(np.array([dotsX, dotsY]).transpose())
            # draw dots
            dotPatch.draw()
            # update target
            if counter < len(targets):
                if  targets[counter] == nTR:
                    if once:
                        mtargetCounter = 0
                        print(counter)
                        once = False
                        once2 = True
                        # below number of target frames? display target!
                    if mtargetCounter < nrOfTargetFrames:
                        # change color fix dot surround to red
                        dotFixSurround.fillColor = [0.5, 0.0, 0.0]
                        dotFixSurround.lineColor = [0.5, 0.0, 0.0]

                    # above number of target frames? dont display target!
                    else:
                        # keep color fix dot surround yellow
                        dotFixSurround.fillColor = [0.5, 0.5, 0.0]
                        dotFixSurround.lineColor = [0.5, 0.5, 0.0]
                        if once2:

                            counter = counter + 1
                            once2 = False

            # update mtargetCounter
            mtargetCounter = mtargetCounter + 1
            #print(mtargetCounter)

            # draw fixation point surround
            dotFixSurround.draw()

            # draw fixation point
            dotFix.draw()

            # draw control text
            # controlText.setText(clock.getTime())
            # controlText.draw()

            myWin.flip()

            # handle key presses each frame
            for keys in event.getKeys():
                if keys in ['escape', 'q']:
                    #myWin.close()
                    #core.quit()
                    pass
                elif keys in ['1']:
                    TimeKeyPressedArray = np.append([TimeKeyPressedArray],[nTR])
                    logging.data(msg='Key1 pressed')
                elif keys == TRIGGER_KEY:
                    nTR = nTR +1
                    once = True
                    logging.data('TR ' + str(nTR))

        # update counter
        i = i + 1

    else:
        runExp = False
# log end of run
print('endofrun')
logging.data('EndOfRun' + str(expInfo['run']))

#%%"""TARGET DETECTION RESULTS"""
# calculate target detection results
# create an array 'targetDetected' for showing which targets were detected
targetDetected = np.zeros(len(targets))
if len(TimeKeyPressedArray) == 0:
    # if no buttons were pressed
    print("No keys were pressed/registered")
    targetsDet = 0
else:
    # if buttons were pressed:
    for index, target in enumerate(targets):
        for TimeKeyPress in TimeKeyPressedArray:
            if (float(TimeKeyPress) >= float(target) and
                    float(TimeKeyPress) <= float(target)+1): #target TR and the follwing TR (have 2 TRs to respond)
                targetDetected[index] = 1

# logging.data('ArrayOfDetectedTargets' + unicode(targetDetected))
logging.data('ArrayOfDetectedTargets' + str(targetDetected))
print("Array Of Detected Targets:")
print(targetDetected)

# number of detected targets
targetsDet = sum(targetDetected)
# logging.data('NumberOfDetectedTargets' + unicode(targetsDet))
logging.data('NumberOfDetectedTargets' + str(targetsDet))
# detection ratio
DetectRatio = targetsDet/len(targetDetected)
# logging.data('RatioOfDetectedTargets' + unicode(DetectRatio))
logging.data('RatioOfDetectedTargets' + str(DetectRatio))

# display target detection results to participant
resultText = 'You have detected %i out of %i targets.' % (targetsDet,
                                                          len(targets))
print(resultText)
logging.data(resultText)
# also display a motivational slogan
if DetectRatio >= 0.95:
    feedbackText = 'Excellent! Keep up the good work'
elif DetectRatio < 0.95 and DetectRatio >= 0.85:
    feedbackText = 'Well done! Keep up the good work'
elif DetectRatio < 0.85 and DetectRatio >= 0.65:
    feedbackText = 'Please try to focus more'
else:
    feedbackText = 'You really need to focus more!'

targetText = visual.TextStim(
    win=myWin,
    color='white',
    height=0.5,
    pos=(0.0, 0.0),
    autoLog=False,
    )
targetText.setText(resultText+feedbackText)
logFile.write(str(resultText) + '\n')
logFile.write(str(feedbackText) + '\n')
targetText.draw()
myWin.flip()
core.wait(5)
myWin.close()

#%%"""SAVE DATA"""
# log important parameters
try:
    logFile.write('TargetDuration=' + str(targetDur) + '\n')
    logFile.write('TimeKeyPressedArray=' + str(TimeKeyPressedArray) + '\n')
except:
    print("(Important parameters could not be logged.)")

# create a pickle file with important arrays
try:
    os.chdir(outFolderName)
    # create python dictionary containing important arrays
    output = {'ExperimentName': expInfo['expName'],
              'Date': expInfo['date'],
              'SubjectID': expInfo['participant'],
              'Run_Number': expInfo['run'],
              'Conditions': conditions,
              'Durations': durations,
              'KeyPresses': TimeKeyPressedArray,
              'DetectedTargets': targetDetected,
              }
    # save dictionary as a pickle in output folder
    #misc.toFile(outFileName + '.pickle', output) # we don't need pickle, csv is enough
    #print('Pickle data saved as: ' + outFileName + '.pickle')
    # Save the same dictionary as a CSV file
    with open(outFileName + '.csv', mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        
        # Write the headers (keys of the dictionary)
        csv_writer.writerow(output.keys())
        
        # Write the values (values of the dictionary)
        # Ensure lists or arrays are joined into a string if needed
        csv_writer.writerow([", ".join(map(str, v)) if isinstance(v, (list, tuple)) else v for v in output.values()])
    
    print("***")
    os.chdir(_thisDir)
except:
    print('(OUTPUT folder could not be created.)')

# EYETRACKER CLOSE DISPLAY AND SAVE EDF
os.chdir(parentDir)
myWin.close()

#%%"""FINISH"""
core.quit()