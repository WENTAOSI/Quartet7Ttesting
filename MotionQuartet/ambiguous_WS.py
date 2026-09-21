"""
Created on Sat May 07 2022
It presents ambiguous motion quartet stimulus.
MotQuart lastes for 80s, FlickQuart lastes for 16s.
n repetitions depends on TRs
The stimulus begins and ends with a fixation condition of 12s.

Total triggers: depends on TR
Total time: depends on TR

Psychopy3 (v2020.2.4)
Based on https://github.com/MSchnei/motion_quartet_scripts (@author: Marian.Schneider)

adopted from @author: Alessandra Pizzuti adopted from @author Marian.Schneider

@author: siwentao 
"""
from psychopy import visual, event, core, monitors, logging, gui, data, misc, sound
import numpy as np
import pandas as pd
import os
import sys
import time
#%% SET PARAMS
###############################################################################
# GUI Store info about experiment and experimental run
expName = 'Amb_MotQuart'  # set experiment name here
expInfo = {
    'run': '1',
    'participant': 'sub-test',
    'display': ['Vanderbilt7T', 'dbic'],
    'aspect_ratio': '1.19',
    'TR': '2',
    'ho_dva': '-0.0981',
    'vo_dva': '1.7652'
    }
# Create GUI at the beginning of exp to get more expInfo
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName, sortKeys=False)
if dlg.OK == False: core.quit()  # user pressed cancel
TRIGGERKEY = 'quoteleft'
# specify vertical or horizontal switch buttom 
vertical_buttom = "1"; horizontal_buttom = "2"; ITI_buttom = "4"
# set global offset
ho_dva = float(expInfo["ho_dva"]); vo_dva = float(expInfo["vo_dva"])
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
###############################################################################
TR = float(expInfo['TR'])
print(f'TR = {TR}')
# if integer TR we can set percise timing 
if expInfo['TR'] == '2':
    fix_TR = 6; 
    #flicker_quartet = 8; 
    amb_quartet = 40
    #DurElem = np.array([fix_TR, flicker_quartet, amb_quartet])
    DurElem = np.array([fix_TR, amb_quartet]) # fix = 12s; flickerQuartet = 16s, AmbiguousQuartet = 80s + 16s = 96s
    NumQuartets = 3  # set number of repetitions of quartet blocks
    # NOTE: Fixation at the beginning and at the end lasts both for 6 triggers.
    # fixation = 0; flicker = 1; quartet = 2
else: # error has to be 2 sec TR
    raise ValueError("TR must be 2 seconds for this experiment.")
   
#total_TR = int((fix_TR * 2) + (amb_quartet + flicker_quartet) * NumQuartets)
total_TR = int(fix_TR + (amb_quartet + fix_TR) * NumQuartets)
print(f"total_TR: {total_TR}")

Conditions = np.zeros(int(NumQuartets * 2))
Conditions[::2] = 2
Conditions[1::2] = 0
Conditions = np.hstack(([0], Conditions))
Durations = np.zeros(len(Conditions))
Durations[Conditions == 0] = 6           # fixation
Durations[Conditions == 2] = amb_quartet # ambiguous quartet duration
print('Conditions:', Conditions)
print('Durations:', Durations)
print('Total number of triggers:', np.sum(Durations))

# Circle properties
circle_dva = 6  # Diameter in dva
circle_radius = circle_dva / 2  # Radius in dva
# ASPECT RATIO (Height/ width)
aspect_ratio = float(expInfo['aspect_ratio'])
# Define aspect ratio (Width / Height)
aspect_ratio = 1/aspect_ratio 
# Calculate the maximum width 
# this solved from np.sqr((width/2)**2 + ((width/aspect_ratio)/2)**2) <= radius
max_width = 2 * circle_radius / np.sqrt(1 + (1 / aspect_ratio**2))
# Calculate the corresponding height
max_height = max_width / aspect_ratio
# specify vertical distance for this participant determined psychophysically outside of the scanner
VertiDist = max_height / 2
HoriDist = max_width / 2
# specificy background color
backColor = [-0.5, -0.5, -0.5]  # from -1 (black) to 1 (white)
# specificy square color
squareColor = np.multiply(backColor, -1)  # from -1 (black) to 1 (white)
#%% SAVING and LOGGING
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
# get the path that this script is in and change dir to it
_thisDir = os.path.dirname(os.path.abspath(__file__))  # get current path
parentDir = os.path.dirname(_thisDir)
os.chdir(parentDir)  # change directory to this path
# Name and create specific subject folder
subjFolderName = '%s_SubjData' % (expInfo['participant'])
if not os.path.isdir(subjFolderName):
    os.makedirs(subjFolderName)
# Name and create data folder for the experiment
dataFolderName = subjFolderName + os.path.sep + '%s' % (expInfo['expName'])
if not os.path.isdir(dataFolderName):
    os.makedirs(dataFolderName)
# Name and create specific folder for logging results
logFolderName = dataFolderName + os.path.sep + 'Logging'
if not os.path.isdir(logFolderName):
    os.makedirs(logFolderName)
logFileName = logFolderName + os.path.sep + '%s_%s_Run%s_%s' % (
    expInfo['participant'], expInfo['expName'], expInfo['run'],
    expInfo['date'])
# Name and create specific folder for output and protocol files
outFolderName = dataFolderName + os.path.sep + 'Output'
if not os.path.isdir(outFolderName):
    os.makedirs(outFolderName)
prtFolderName = dataFolderName + os.path.sep + 'Protocol'
if not os.path.isdir(prtFolderName):
    os.makedirs(prtFolderName)
# save a log file and set level for msg to be received
logFile = logging.LogFile(logFileName+'.log', level=logging.INFO)
logging.console.setLevel(logging.WARNING)  # set console to receive warnings
#%% MONITOR AND WINDOW
# source:https://www.dartmouth.edu/dbic/research_infrastructure/peripherals.html
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
moni.setSizePix([PixW, PixH])  # [1920.0, 1080.0] in psychoph lab
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
myWin = visual.Window(size=(PixW, PixH),
                      screen = screen,
                      winType='pyglet',  # winType : None, ‘pyglet’, ‘pygame’
                      allowGUI=False,
                      allowStencil=False,
                      fullscr=True,  # for psychoph lab: fullscr = True
                      monitor=moni,
                      color=backColor,
                      colorSpace='rgb',
                      units='deg',
                      blendMode='avg',
                      )
myWin.mouseVisible = False
logFile.write('Conditions=' + str(Conditions) + '\n')
logFile.write('Durations (Triggers) =' + str(Durations) + '\n')
# create array to log key pressed events
KeyPressedArray = np.array(['KeyPressed', 'KeyPressedt'])

# %% STIMULI
SquareSize = 1.0  # 1.1 #1.8
SquareDur = 0.15  # in seconds # 9 frames
BlankDur = 0.067  # in seconds # 5 frames
logFile.write(f'Durations : {Durations}' + '/n')
logFile.write('SquareSize=' + str(SquareSize) + '\n')
logFile.write('SquareDur=' + str(SquareDur) + '\n')
logFile.write('BlankDur=' + str(BlankDur) + '\n')

message = visual.TextStim(myWin,text='Condition',pos=(-16, -8))
dotFix = visual.Circle(
    myWin,autoLog=False,name='dotFix',units='deg',radius=.15,
    pos=apply_global_offset((0,0), global_offset),fillColor='red', lineColor='red'
    )
Square = visual.GratingStim(
    myWin,autoLog=False,name='Square',tex=None,units='deg',
    size=(SquareSize, SquareSize),color=squareColor
    )
# Four Circles
circle_size = 1  # diameter of each circle
if expInfo["display"] == 'dbic':
    positions = [apply_global_offset((-3.5, 4), global_offset), apply_global_offset((-2, 4), global_offset),\
                 apply_global_offset((2, 4), global_offset), apply_global_offset((3.5, 4), global_offset)]  # Anchored positions 1, 2, 3, 4
elif expInfo["display"] == 'Vanderbilt7T':
    positions = [apply_global_offset((-3.5, 4), global_offset), apply_global_offset((-1.25, 4.5), global_offset),\
                 apply_global_offset((1.25, 4.5), global_offset), apply_global_offset((3.5, 4), global_offset)]  # Anchored positions 1, 2, 3, 4
# Generate circle objects at the specified positions
circles = []
for pos in positions:
    # Create a circle object at each specified position
    circle = visual.Circle(
        myWin, autoLog=False, units='deg',
        radius=circle_size/2,  # Diameter of the circle
        pos=pos,  # Position from the list
        lineColor=squareColor, lineWidth=2, fillColor=None
    )
    circles.append(circle)
# Generate H/V letters in the circle
Hs = []; Vs = []; Fs = []
for pos in positions:
    H = visual.TextStim(win=myWin, autoLog=False,color='white',height=circle_size-0.2,text='H',pos=pos)
    V = visual.TextStim(win=myWin,color='white',height=circle_size-0.2,text='V',pos=pos)
    #F = visual.TextStim(win=myWin,color='white', height=circle_size-0.2, text='F', pos=pos)
    Hs.append(H); Vs.append(V); 
    #Fs.append(F)
triggerText = visual.TextStim(
    win=myWin,color='white',height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text='Experiment will start soon. Waiting for scanner'
    )
instructText = visual.TextStim(
    win=myWin, color='white',height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text=f'Press {vertical_buttom} when you perceive VERTICAL\n\
        Press {horizontal_buttom} when you perceive HORIZONTAL\n\
            PRESS on the keys to practice'
    )
'''
instruct_ITI = visual.TextStim(
    win=myWin,color='white',height=0.5,
    pos=apply_global_offset(base_pos=(0,-3.5),global_offset=global_offset),
    text=f'Press {ITI_buttom} when you perceive Four Squares FLASHING\n\Press the Key to Practice')
'''
anykeyText = visual.TextStim(
    win=myWin, color='white',height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text='Good Job!\n\Press on any key to continue'
    )    
# %% TIME AND TIMING PARAMeTERS
# parameters
totalTrigger = np.sum(Durations)
print('Total number of triggers:', totalTrigger)
# get screen refresh rate
refr_rate = myWin.getActualFrameRate()  # get screen refresh rate
if refr_rate is not None:
    frameDur = 1.0/round(refr_rate)
else:
    frameDur = 1.0/60.0
    refr_rate = 60.0
logFile.write('RefreshRate=' + str(refr_rate) + '\n')
logFile.write('FrameDuration=' + str(frameDur) + '\n')
# define clock
clock = core.Clock()
logging.setDefaultClock(clock)
# %% FUNCTIONS
# create necessary functions for quartet and flicker
def quartetPart1(Hori, Verti):
    Square.setPos(apply_global_offset((-Hori, Verti), global_offset))
    Square.draw()
    Square.setPos(apply_global_offset((Hori, -Verti), global_offset))
    Square.draw()
    dotFix.draw()
def quartetPart2(Hori, Verti):
    Square.setPos(apply_global_offset((Hori, Verti), global_offset))
    Square.draw()
    Square.setPos(apply_global_offset((-Hori, -Verti), global_offset))
    Square.draw()
    dotFix.draw()
def fixation():
    dotFix.draw()
    myWin.flip()
def quartet(Hori, Verti):
    NumSquareFrames = int(round(SquareDur/frameDur))
    NumBlankFrames = int(round(BlankDur/frameDur))
    for frameN in range(NumSquareFrames):
        quartetPart1(Hori, Verti)
        myWin.flip()
    for frameN in range(NumBlankFrames):
        dotFix.draw()
        myWin.flip()
    for frameN in range(NumSquareFrames):
        quartetPart2(Hori, Verti)
        myWin.flip()
    for frameN in range(NumBlankFrames):
        dotFix.draw()
        myWin.flip()
'''
def flickerSl(Hori, Verti, instruct):
    NumSquareFrames = int(round(SquareDur/frameDur))
    NumBlankFrames = 2*int(round(BlankDur/frameDur)) + NumSquareFrames
    for frameN in range(NumSquareFrames):
        Square.setPos(apply_global_offset((-Hori, Verti), global_offset))
        Square.draw()
        Square.setPos(apply_global_offset((Hori, -Verti), global_offset))
        Square.draw()
        Square.setPos(apply_global_offset((Hori, Verti), global_offset))
        Square.draw()
        Square.setPos(apply_global_offset((-Hori, -Verti), global_offset))
        Square.draw()
        dotFix.draw()
        if instruct:
            instruct_ITI.draw()
            for circle in circles:
                circle.draw()
            Fs[int(ITI_buttom)-1].draw()
        myWin.flip()
    for frameN in range(NumBlankFrames):
        dotFix.draw()
        if instruct:
            instruct_ITI.draw()
            for circle in circles:
                circle.draw()  
            Fs[int(ITI_buttom)-1].draw()  
        myWin.flip()
'''
def buttom_instruct(win,vertical_buttom, horizontal_buttom, ITI_buttom):
    '''
    Displays buttom instructions at the report stage 
    Parameters
    win: (Psychopy object) window setting
    vertical_buttom: (str): "1","2","3",or"4"
    horizontal_buttom: (str): "1","2","3",or"4"
    ITI_buttom:(str): "1","2","3",or"4"
    Returns None
    '''
    # Draw and display the circles
    for circle in circles:
        circle.draw()
    Vs[int(vertical_buttom)-1].draw()   # Because python start counting from 0, draw the first one in the 0th on the list 
    Hs[int(horizontal_buttom)-1].draw() #    
    instructText.draw() # show instruct text
    win.flip()
    # show red color V,H after buttom press
    event.waitKeys(keyList=[vertical_buttom], timeStamped=False)
    for circle in circles:
        circle.draw()
    Vs[int(vertical_buttom)-1].setColor('red')
    Vs[int(vertical_buttom)-1].draw() 
    Hs[int(horizontal_buttom)-1].draw() #    
    instructText.draw() # show instruct text
    win.flip()
    event.waitKeys(keyList=[horizontal_buttom], timeStamped=False)
    for circle in circles:
        circle.draw()
    Hs[int(horizontal_buttom)-1].setColor('red')
    Hs[int(horizontal_buttom)-1].draw()
    Vs[int(vertical_buttom)-1].draw()
    instructText.draw()
    win.flip()
    core.wait(2)
    '''
    # Display instruction for flashing four dots
    while not event.getKeys(keyList=[ITI_buttom]):
        flickerSl(HoriDist, VertiDist, instruct=True) # this includes the instruction text
    '''
    '''
    local_clock = core.Clock()
    #continue to dispay the circles except turning F to red for 2 sec
    while local_clock.getTime() < 2:
        flickerSl(HoriDist, VertiDist, instruct=True)
        for circle in circles:
            circle.draw()
        Fs[int(ITI_buttom)-1].setColor('red')
        Fs[int(ITI_buttom)-1].draw()
        win.flip()
        '''
# %% RENDER_LOOP
# give the system time to settle
core.wait(1)
# Only shows instruction when 
if expInfo['run'] == '1':
    buttom_instruct(myWin,vertical_buttom, horizontal_buttom, None)
    anykeyText.draw()
    myWin.flip()
    event.waitKeys(timeStamped=False)
# wait for scanner trigger
triggerText.draw()
myWin.flip()
# # --------------------------------------------
# # launch: operator selects Scan or Test (emulate); see API documentation
# vol = launchScan(win, MR_settings, globalClock=clock, mode='Test')
# #----------------------------------------------
event.waitKeys(keyList=[TRIGGERKEY], timeStamped=False)
# Create Counters
i = 0             # counter for blocks
trigCount = 0     # counter triggers
# reset clocks
clock.reset()   # Comment out only for the simulation
logging.data('StartOfRun' + str(expInfo['run']))
logging.data(msg='Scanner trigger %i' % (trigCount))
# ============================================================
# RUN CONDITIONS
# ============================================================
last_trigger_time = None
while trigCount < totalTrigger:
    logging.data('StartOfCondition' + str(Conditions[i]))
    # --------------------------------------------------------
    # RUN CURRENT CONDITION
    # --------------------------------------------------------
    while trigCount < np.sum(Durations[0:i+1]):
        t = clock.getTime()
        # Draw stimulus
        if Conditions[i] == 0:
            fixation()
        elif Conditions[i] == 1:
            fixation()
        elif Conditions[i] == 2:
            quartet(HoriDist,VertiDist)
        # ----------------------------------------------------
        # CHECK KEYS
        # ----------------------------------------------------
        for key in event.getKeys():
            if key in ['escape', 'q']:
                logging.data(msg='User pressed quit')
                myWin.close()
                core.quit()
            elif key in ['1', 'num_1']:
                t = clock.getTime()
                KeyPressed = '1'
                KeyPressedNew = np.array([KeyPressed,t])
                KeyPressedArray = np.vstack(( KeyPressedArray,KeyPressedNew))
                logging.data(msg='Key1 pressed')
            elif key in ['2', 'num_2']:
                t = clock.getTime()
                KeyPressed = '2'
                KeyPressedNew = np.array([KeyPressed,t])
                KeyPressedArray = np.vstack((KeyPressedArray,KeyPressedNew))
                logging.data(msg='Key2 pressed')
            elif key in ['3', 'num_3']:
                t = clock.getTime()
                KeyPressed = '3'
                KeyPressedNew = np.array([KeyPressed,t])
                KeyPressedArray = np.vstack((KeyPressedArray,KeyPressedNew))
                logging.data(msg='Key3 pressed')
            elif key in ['4', 'num_4']:
                t = clock.getTime()
                KeyPressed = '4'
                KeyPressedNew = np.array([KeyPressed,t])
                KeyPressedArray = np.vstack((KeyPressedArray,KeyPressedNew))
                logging.data(msg='Key4 pressed')
            elif key == TRIGGERKEY:
                # Time of this scanner trigger
                last_trigger_time = clock.getTime()
                # Count this trigger
                trigCount += 1
                logging.data(msg='Scanner trigger %i' % trigCount)
    # Finished current condition
    i += 1
    # ============================================================
    # FINAL TR
    # ============================================================
    # The final scanner trigger starts the final TR.
    # The final condition is always fixation.
    # Since there is no next trigger, use the known TR duration
    # to allow the final TR to elapse.
    if last_trigger_time is not None:
        while clock.getTime() - last_trigger_time < TR:
            fixation()
            # Still check for quit
            for key in event.getKeys():
                if key in ['escape', 'q']:
                    logging.data(msg='User pressed quit')
                    myWin.close()
                    core.quit()
# END RUN
logging.data('EndOfRun' + str(expInfo['run']) + '\n')
# %% SAVE DATA
# calculate speed [degrees per frame]
TravelTime = 2*int(round(SquareDur/frameDur))+2*int(round(BlankDur/frameDur))
logFile.write('TravelTime=' + str(TravelTime) + '\n')
HoriSpeed = (HoriDist*4)/TravelTime
logFile.write('HoriSpeed=' + str(HoriSpeed) + '\n')
VertiSpeed = (VertiDist*4)/TravelTime
logFile.write('HoriSpeed=' + str(HoriSpeed) + '\n')
logFile.write('VertiSpeed=' + str(VertiSpeed) + '\n')
logFile.write('HoriDist=' + str(HoriDist) + '\n')
logFile.write('VertiDist=' + str(VertiDist) + '\n')

# Change into output folder
os.chdir(outFolderName)
# Define a mapping for conditions
condition_labels = {
            '0': 'fixation',
            vertical_buttom: 'vertiM',
            horizontal_buttom: 'horiM',
            ITI_buttom: 'flickerSl'
        }
# Skip the header row and map the labels
labels = ['Label'] + [condition_labels.get(row[0], 'Unknown') for row in KeyPressedArray[1:]]
# Add the labels as a new column
KeyPressedArray = np.column_stack((KeyPressedArray, labels))
# First row contains the column names
keyPressed_df = pd.DataFrame(KeyPressedArray[1:],columns=KeyPressedArray[0])
# Convert timestamp column back to numeric
keyPressed_df["KeyPressedt"] = pd.to_numeric(keyPressed_df["KeyPressedt"],errors="raise")
# Save np.array into output folder 
np.save(f"{expInfo['participant']}_amb_run{expInfo['run']}_key_presses.npy", KeyPressedArray)

# Save the array to a CSV file for visual inspection 
np.savetxt(
    f"{expInfo['participant']}_amb_run{expInfo['run']}_key_presses.csv",            # File name
    KeyPressedArray,              # Data to save
    fmt='%s',                     # Format: string for all columns
    delimiter=",",                # CSV delimiter
    comments=''                   # Prevent '#' before the header
)
################################### SAVE BIDS event files ########################################
os.chdir(parentDir)
amb_run_duration = total_TR * TR
amb_fixation_duration = fix_TR * TR
amb_final_fixation_onset = amb_run_duration - amb_fixation_duration
# Rename stimulus labels
condition_mapping = {
    "fixation": "fixation",
    "vertiM": "vertical_motion",
    "horiM": "horizontal_motion",
    "flickerSI": "flicker_static",
}
# Build BIDS events
amb_events = (keyPressed_df[["KeyPressedt", "Label"]].copy().rename(columns={"KeyPressedt": "onset", "Label": "trial_type"}))
amb_events["onset"] = pd.to_numeric(amb_events["onset"],errors="raise",)
amb_events = (amb_events.sort_values("onset").reset_index(drop=True))
amb_events["trial_type"] = (amb_events["trial_type"].astype(str).str.strip().replace(condition_mapping))
# Duration = next onset - current onset
amb_events["duration"] = (amb_events["onset"].shift(-1) - amb_events["onset"])
# Last logged event ends when final fixation begins
amb_events.loc[amb_events.index[-1],"duration",] = (amb_final_fixation_onset - amb_events.loc[amb_events.index[-1],"onset",])
# Add initial fixation
initial_fixation = pd.DataFrame({"onset": [0.0], "duration": [amb_fixation_duration], "trial_type": ["fixation"]})
# Add final fixation
final_fixation = pd.DataFrame({"onset": [amb_final_fixation_onset],"duration": [amb_fixation_duration],"trial_type": ["fixation"]})
amb_events = pd.concat([initial_fixation, amb_events,final_fixation,],ignore_index=True,)
# Sort once more
amb_events = (amb_events.sort_values("onset").reset_index(drop=True))
# Validate
if (amb_events["duration"] <= 0).any():
    bad_rows = amb_events.loc[amb_events["duration"] <= 0]
    raise ValueError(f"Negative durations found:\n{bad_rows}")
BIDS_dir = os.path.join('BIDS_events', expInfo['participant'], 'func')
os.makedirs(BIDS_dir, exist_ok=True)
amb_output_file = os.path.join(BIDS_dir, f"{expInfo['participant']}_task-ambiguous_run-{int(expInfo['run']):02d}_events.tsv")
# Save
amb_events.to_csv(amb_output_file, sep="\t", index=False, float_format="%.3f")
print(f"Saved {amb_output_file}")
'''
# Change into protocol folder
os.chdir(parentDir)
os.chdir(prtFolderName)
## construct protocol file for AMB

# set key KeyPressedArray to pd.DataFrame
KeyPressed_df = pd.DataFrame(KeyPressedArray[1:], columns=KeyPressedArray[0])
# create timestamp as stop time for events 
KeyPressed_df['Timestamp'] = KeyPressed_df['KeyPressedt'].shift(-1)
# define flicker block ends
if TR == 2:
    flicker_ends = [108, 204, 300]
elif TR == 4.217:
    flicker_ends = [168.68, 320.492]
# find all flickerSl rows
flicker_indices = KeyPressed_df[KeyPressed_df['Label'] == 'flickerSl'].index
# loop through and assign the fixed values
for i, flicker_end in zip(flicker_indices, flicker_ends):
    flicker_start = flicker_end - 16
    KeyPressed_df.at[i, 'Timestamp'] = flicker_end #assign flicker end
    KeyPressed_df.at[i-1, 'KeyPressedt'] = flicker_start # assign flicker start a block before
    KeyPressed_df.at[i, 'KeyPressedt'] = flicker_start
# Add fixation row at the top and bottom
if TR == 2:
    fixation_start = pd.DataFrame({
    'KeyPressed': [0],
    'KeyPressedt': [0.0],
    'Label': ['fixation'],
    'Timestamp': [12.0]
    })

    fixation_end = pd.DataFrame({
    'KeyPressed': [0],
    'KeyPressedt': [596.0],
    'Label': ['fixation'],
    'Timestamp': [312.0]
    })
elif TR == 4.217:
    fixation_start = pd.DataFrame({
        'KeyPressed': [0],
        'KeyPressedt': [0.0],
        'Label': ['fixation'],
        'Timestamp': [16.868]
    })

    fixation_end = pd.DataFrame({
    'KeyPressed': [0],
    'KeyPressedt': [0],
    'Label': ['fixation'],
    'Timestamp': [337.36]
    })
# Concat start row + main dataframe + end row
KeyPressed_df = pd.concat([fixation_start, KeyPressed_df, fixation_end], ignore_index=True)

# add Duration, Onset change Label to Stim
KeyPressed_df['Timestamp'] = KeyPressed_df['Timestamp'].astype(float)
KeyPressed_df['Duration'] = KeyPressed_df['Timestamp'].diff().fillna(KeyPressed_df['Timestamp'])
KeyPressed_df['Stim'] = KeyPressed_df['Label']
KeyPressed_df['Duration'] = KeyPressed_df['Duration'].astype(float)
KeyPressed_df['Onset'] = KeyPressed_df['Timestamp'] - KeyPressed_df['Duration']

KeyPressed_df.to_csv(f"{expInfo['participant']}_amb_run{expInfo['run']}_protocol.csv", index=False)
'''
os.chdir(parentDir)
myWin.close()
# %% FINISH#
core.quit()