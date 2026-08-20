"""
It presents physical motion quartet stimulus.

Psychopy3 (v2024.2.4)
Based on https://github.com/MSchnei/motion_quartet_scripts (@author: Marian.Schneider)
adopted from @author: Alessandra Pizzuti adopted from @author Marian.Schneider

@author: siwentao
modfied by: siwentao to create V2 of physical motion quartet stimulus for 7T fMRI testing
Key modifications:
- removed not needed eyetracker to make it cleaner 
- Flicker baseline becomes fixation baseline
- fixation task on fixation dot (press 1 when it turns white)
- switch H and V position psudorandomly while preserving the total number of trials for each block 

"""
from psychopy import visual, event, core, monitors, logging, gui, data, misc
from itertools import cycle
import numpy as np
import os
import sys
import time
import pandas as pd

# %% SET PARAMS
###############################################################################
# %% GUI
# Store info about experiment and experimental run
expName = 'Phy_MotQuart'  # set experiment name here
expInfo = {
    'run': '1',
    'participant': 'sub-test',
    'display': ['Vanderbilt7T', 'dbic'],
    'aspect_ratio': '1.12',
    'TR': '2',
    'ho_dva': '-0.0981',
    'vo_dva': '1.7652'
    }
# Create GUI at the beginning of exp to get more expInfo
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName, sortKeys=False)
if dlg.OK == False: core.quit()  # user pressed cancel
TRIGGERKEY = 'quoteleft'
# BLOCK DURATIONS [in TR]
# set durations of conditions and baseline
# specify vertical or horizontal condition number 
vertical_cond = "1"
horizontal_cond = "2"
fixation_cond = "4"
# set global offset
ho_dva = float(expInfo["ho_dva"])
vo_dva = float(expInfo["vo_dva"])
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
# %% BLOCK DURATIONS [Triggers]
# set durations of conditions and baseline
if TR == 4.217:
    MotionDur = 4
    BaseDur = 4
    Fixation = 4
    NumOf12PerBlock = 8 # number of (hor + ver) per block
    NumQuartets = 2

elif TR == 1.612:
    MotionDur = 8     # Vertical or Horizontal motion 8 TR 
    BaseDur = 10     # 4 sqares flickering 10 TR
    Fixation = 8     # Fixation (beginning and end) 8 TR
    NumOf12PerBlock = 6 # number of (hor + ver) per block
    NumQuartets = 4 # number of cycles

elif TR == 2:
    MotionDur = 5     # Vertical or Horizontal motion 5 TR 
    BaseDur = 6     # formerly 4 sqares flickering 8 TR, now fixation 
    Fixation = 6     # Fixation (beginning and end) 6 TR
    NumOf12PerBlock = 8 # number of (hor + ver) per block
    NumQuartets = 3 # number of cycles
    print(MotionDur)

else:
    MotionDur = 5     # Vertical or Horizontal motion 5 TR 
    BaseDur = 6     # formerly 4 sqares flickering 8 TR, now fixation
    Fixation = 6     # Fixation (beginning and end) 6 TR
    NumOf12PerBlock = 8 # number of (hor + ver) per block
    NumQuartets = 3 # number of cycles
    print(MotionDur)

# ============================================================
# RANDOMIZE H/V ORDER WITHIN EACH BIG BLOCK
# ============================================================
# set number of repetitions for each condition
# fixation = 0; horiM = 2; vertiM = 1; baseline fixation = 4
rng = np.random.default_rng()
all_motion_conditions = []
for block_i in range(NumQuartets):
    # Equal number of vertical and horizontal trials
    block_conditions = np.array(
        [1] * (NumOf12PerBlock // 2) +
        [2] * (NumOf12PerBlock // 2)
    )
    # Randomize H/V order within this big block
    rng.shuffle(block_conditions)
    all_motion_conditions.extend(block_conditions)
Conditions = np.array(all_motion_conditions, dtype=int)
print(f"randomized motion conditions: {Conditions}")

# Insert baseline fixation after each big motion block
pos_baseline = np.arange(
    NumOf12PerBlock,
    NumOf12PerBlock * NumQuartets,
    NumOf12PerBlock
)
Conditions = np.insert(Conditions,pos_baseline,4)

# Add beginning and ending fixation
Conditions = np.hstack(([0], Conditions, [0]))
Conditions = Conditions.astype(int)
print(f"final conditions: {Conditions}")

# ============================================================
# FIXED CONDITION DURATIONS
# ============================================================
Durations = np.ones(len(Conditions),dtype=int) * MotionDur
# Inter-block fixation
Durations[Conditions == 4] = BaseDur
# Beginning/end fixation
Durations[Conditions == 0] = Fixation
totalTrigger = np.sum(Durations)
print("Durations:")
print(Durations)
print("Total TRs:", totalTrigger)

# ============================================================
# FIXATION COLOR-CHANGE TIMELINE
# ============================================================
# Red fixation normally; occasionally turns white for exactly 1 TR.
FIXATION_NORMAL_COLOR = 'red'
FIXATION_TARGET_COLOR = 'white'
# Number of white-fixation events in one run
N_WHITE_FIXATIONS = 10
# Minimum distance between white events, in TRs
MIN_WHITE_GAP = 6
# Subject must respond within one TR of target onset
DETECTION_WINDOW = TR
rng_fix = np.random.default_rng()
# Avoid very beginning/end of run
candidate_TRs = np.arange(4, totalTrigger - 4)
white_fix_TRs = []
while len(white_fix_TRs) < N_WHITE_FIXATIONS:
    candidate = int(rng_fix.choice(candidate_TRs))
    # Prevent white events from clustering together
    if all(abs(candidate - existing) >= MIN_WHITE_GAP for existing in white_fix_TRs):
        white_fix_TRs.append(candidate)

white_fix_TRs = np.sort(np.array(white_fix_TRs))
print(f"White fixation TRs: {white_fix_TRs}")

# %%Circle properties
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
# specify vertical distance for this participant 
VertiDist = max_height / 2
HoriDist = max_width / 2
# specificy background color
backColor = [-0.5, -0.5, -0.5]  # from -1 (black) to 1 (white)
# specificy square color
squareColor = np.multiply(backColor, -1)  # from -1 (black) to 1 (white)

# %% SAVING and LOGGING
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
# Name and create specific folder for protocol files
prtFolderName = dataFolderName + os.path.sep + 'Protocols'
if not os.path.isdir(prtFolderName):
    os.makedirs(prtFolderName)
# Name and create folder for Output keyresponses
outFolderName = dataFolderName + os.path.sep + 'Output'
if not os.path.isdir(outFolderName):
    os.makedirs(outFolderName)
print(f"outFolderName: {outFolderName}")
# save a log file and set level for msg to be received
logFile = logging.LogFile(logFileName+'.log', level=logging.INFO)
logging.console.setLevel(logging.WARNING)  # set console to receive warnings
logFile.write('Conditions=' + str(Conditions) + '\n')
logFile.write('Durations (Triggers) =' + str(Durations) + '\n')
# %% MONITOR AND WINDOW
# set monitor information:
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

# %% STIMULI
# INITIALISE SOME STIMULI
SquareSize = 1.0  # 1.1 #1.8
SquareDur = 0.15  # in seconds # 9 frames
BlankDur = 0.067  # in seconds # 5 frames

logFile.write('SquareSize=' + str(SquareSize) + '\n')
logFile.write('SquareDur=' + str(SquareDur) + '\n')
logFile.write('BlankDur=' + str(BlankDur) + '\n')
logFile.write(f'Durations: {Durations} \n')

message = visual.TextStim(myWin,text='Condition',pos=apply_global_offset((-16, -8), global_offset))

dotFix = visual.Circle(myWin,autoLog=False,name='dotFix',units='deg',radius=.15,
                       pos=apply_global_offset((0,0), global_offset))

def update_fixation_color(current_TR):
    if current_TR in white_fix_TRs:
        color = FIXATION_TARGET_COLOR
    else:
        color = FIXATION_NORMAL_COLOR
    dotFix.fillColor = color
    dotFix.lineColor = color
    
# Keep the dot the same with lower resolution e.g. vanderbilt 7T screen
#if expInfo['display'] == 'Vanderbilt7T':
#    dotFix.radius = int(10/(1920/1024))

Square = visual.GratingStim(myWin,autoLog=False,name='Square',tex=None,units='deg',
                            size=(SquareSize, SquareSize),color= squareColor)
    
triggerText = visual.TextStim(
    win=myWin,color='white',height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text='Experiment will start soon. Waiting for scanner'
    )

instructText = visual.TextStim(
    win=myWin,color='white',height=0.5,
    pos=apply_global_offset(base_pos=(0, 0), global_offset=global_offset),
    text="Press 1 when the fixation dot turns white.\n\nPress '1' to start the experiment."
)

# %% TIME AND TIMING PARAMeTERS
# get screen refresh rate
refr_rate = myWin.getActualFrameRate()  # get screen refresh rate
if refr_rate is not None:
    frameDur = 1.0/round(refr_rate)
else:
    frameDur = 1.0/60.0  # couldn't get a reliable measure so guess
logFile.write('RefreshRate=' + str(refr_rate) + '\n')
logFile.write('FrameDuration=' + str(frameDur) + '\n')

# define clock
clock = core.Clock()
logging.setDefaultClock(clock)

# %% FUNCTIONS
# create necessary functions for quartet and flicker
NumSquareFrames = int(round(SquareDur/frameDur))  # num of square frames
NumBlankFrames = int(round(BlankDur/frameDur))  # num of blank frames
TravelTime = 2*NumSquareFrames+2*NumBlankFrames
TravelTimeArray = np.arange(TravelTime+1)/TravelTime
HoriTimeCycle = cycle(TravelTimeArray)  # iterate through the TravelTimeArray
VertiTimeCycle = cycle(TravelTimeArray)  # iterate through the TravelTimeArray

def HMotion_update(Hori, Verti):
    x = next(HoriTimeCycle)  # pick next element when cycling TravelTimeArray
    mHori = np.cos((2*np.pi)*x)*Hori
    Square.setPos(apply_global_offset((mHori, Verti), global_offset))  # square northeast
    Square.draw()
    Square.setPos(apply_global_offset((-mHori, -Verti), global_offset))  # square southwest
    Square.draw()
    dotFix.draw()
    myWin.flip()
    return mHori

def VMotion_update(Hori, Verti):
    x = next(VertiTimeCycle)
    mVerti = np.cos((2*np.pi)*x)*Verti
    Square.setPos(apply_global_offset((-Hori, mVerti), global_offset))  # square northwest
    Square.draw()
    Square.setPos(apply_global_offset((Hori, -mVerti), global_offset))  # square southeast
    Square.draw()
    dotFix.draw()
    myWin.flip()
    return mVerti

# %% RENDER_LOOP
# create array to log key pressed events
KeyPressedArray = np.array(['KeyPressedt', 't'])
ButtonPressTimes = []

# give the system time to settle
core.wait(1)
# instructions for the participant
instructText.draw()
myWin.flip()
event.waitKeys(keyList=['1','num_1'], timeStamped=False)

# wait for scanner trigger
triggerText.draw()
myWin.flip()
# # --------------------------------------------
# launch: operator selects Scan or Test (emulate); see API documentation
# vol = launchScan(win, MR_settings, globalClock=clock, mode='Scan')
# #----------------------------------------------
# reset clocks
event.waitKeys(keyList=[TRIGGERKEY], timeStamped=False)
clock.reset()
# Create Counters
i = 0           # counter for blocks
trigCount = 0   # counter triggers

# reset clocks
# clock.reset()
print(totalTrigger)
logging.data('StartOfRun' + str(expInfo['run']))
logging.data(msg='Scanner trigger %i' % (trigCount))

while trigCount < totalTrigger:    # 

    logging.data('StartOfCondition'+ str(Conditions[i]))

    while trigCount < np.sum(Durations[0:i+1]):
        t = clock.getTime()
        
        # update fixation color based on the current trigger count
        update_fixation_color(trigCount)
        
        if Conditions[i] == 0:
            dotFix.draw()
            myWin.flip()
        elif Conditions[i] == 2:
            mHori = HMotion_update(HoriDist, VertiDist)
        elif Conditions[i] == 1:
            mVerti = VMotion_update(HoriDist, VertiDist)
        elif Conditions[i] == 4:
            dotFix.draw()
            myWin.flip()

        for key in event.getKeys():
                if key in ['escape', 'q']:
                    logging.data(msg='User pressed quit')
                    myWin.close()
                    core.quit()
                elif key in ['1', 'num_1']:
                    t = clock.getTime()
                    KeyPressed = '1'
                    KeyPressedNew = np.array([KeyPressed, t])
                    KeyPressedArray = np.vstack((KeyPressedArray,KeyPressedNew))
                    ButtonPressTimes.append(t)
                    logging.data(msg=f'Fixation detection button pressed at {t:.3f} s')
                                   
                elif key == TRIGGERKEY:
                    t = clock.getTime()
                    trigCount = trigCount + 1
                    logging.data(msg='Scanner trigger %i' % (trigCount))
    i = i+1
    print('Block counter: %i' % i)
logging.data('EndOfRun' + str(expInfo['run']) + '\n')

# %% SAVE DATA
 # calculate speed [degrees per frame]
HoriSpeed = (HoriDist*4)/TravelTime
VertiSpeed = (VertiDist*4)/TravelTime
logFile.write('HoriSpeed=' + str(HoriSpeed) + '\n')
logFile.write('VertiSpeed=' + str(VertiSpeed) + '\n')
logFile.write('HoriDist=' + str(HoriDist) + '\n')
logFile.write('VertiDist=' + str(VertiDist) + '\n')
print('horizontalDistance: %f' % HoriDist)
print('verticalDistance: %f' % VertiDist)
print('Button Press Times:')
print(ButtonPressTimes)
print(KeyPressedArray)
print("Durations:")
print(Durations)

output_csv =  pd.DataFrame(KeyPressedArray)
# Save key pressed events to CSV
output_csv.to_csv(f"{outFolderName}/{expInfo['participant']}_Phy_keyPressed_run{expInfo['run']}.csv", index=False)
# Calculate timestamps (cumulative durations)
Timestamps = np.cumsum(Durations) 
print("Timestamps:")
print(Timestamps)
# Define a mapping for conditions
condition_labels = {
        '0': 'fixation',
        vertical_cond: 'vertiM',
        horizontal_cond: 'horiM',
        fixation_cond: 'fixation'
    }
# Map the labels to the Conditions array
Labels = np.array([condition_labels[str(cond)] for cond in Conditions]) # extra step to convert cond to str
print(f'Labels {Labels}')
# Initialize the base with column headers
protocol_array = np.array([['Conditions', 'Durations', 'Timestamp','Stim']]) 
print(f"protocol_array{protocol_array}")
# Combine all columns
protocol_data = np.column_stack((Conditions, Durations, Timestamps, Labels))   
print(f"protocol_data {protocol_data}")
# Stack the data below the header
protocol_array = np.vstack((protocol_array, protocol_data.astype(str)))  # Convert to string for uniformity
print(f"protocol_array {protocol_array}")

# Change into protocol folder
os.chdir(prtFolderName)
 # save protocol first 
#np.save(f"{expInfo['participant']}_Phy_protocol.npy", protocol_array)
# convert protocal_array to pd.dataframe
# Convert the protocol array to a DataFrame starting from the second row for the data
protocol_array_df = pd.DataFrame(protocol_array[1:], columns=protocol_array[0])
# Set Onset time
# Convert 'Timestamp' and 'Durations' to numeric types.
protocol_array_df['Timestamp'] = protocol_array_df['Timestamp'].astype(float)
protocol_array_df['Durations'] = protocol_array_df['Durations'].astype(float)
protocol_array_df['Onset'] = protocol_array_df['Timestamp'] - protocol_array_df['Durations']
# Save protocol_array to CSV for visualization Only one prot across runs 
protocol_array_df.to_csv(f"{expInfo['participant']}_phy_protocol.csv", index=False)
print(protocol_array)
# change into output folder
os.chdir(parentDir)

################################### SAVE BIDS EVENT FILE ########################################
# =============================================================================
# ORGANIZE MAIN MOTION/FIXATION EVENTS
# =============================================================================
BIDS_df = (
    protocol_array_df[["Onset", "Durations", "Stim"]]
    .rename(
        columns={
            "Onset": "onset",
            "Durations": "duration",
            "Stim": "trial_type"
        }
    )
)
condition_mapping = {
    "fixation": "fixation",
    "vertiM": "vertical_motion",
    "horiM": "horizontal_motion",
}
# Convert TR units -> seconds
BIDS_df["onset"] = (pd.to_numeric(BIDS_df["onset"], errors="raise") * TR)
BIDS_df["duration"] = (pd.to_numeric(BIDS_df["duration"], errors="raise") * TR)
BIDS_df["trial_type"] = (BIDS_df["trial_type"].astype(str).str.strip().replace(condition_mapping))

# =============================================================================
# INITIALIZE FIXATION-TASK COLUMNS
# =============================================================================
# target_present:
#   0 = no white fixation target during this condition
#   1 = white fixation target occurred during this condition
# target_onset:
#   exact onset of white fixation in seconds relative to scan onset
# detected:
#   1 = target detected
#   0 = target missed
#   n/a = no target occurred
# response_time:
#   seconds from target onset to button press
BIDS_df["target_present"] = 0
BIDS_df["target_onset"] = "n/a"
BIDS_df["detected"] = "n/a"
BIDS_df["response_time"] = "n/a"

# =============================================================================
# 3. ASSIGN FIXATION TARGETS TO THEIR CORRESPONDING CONDITION ROWS
# =============================================================================
# Track which button presses have already been assigned to a target.
# This prevents one button press from being counted for multiple targets.
used_button_presses = set()
for target_TR in white_fix_TRs:
    # white_fix_TRs is 1-indexed:
    # TR 1 -> 0 sec
    # TR 2 -> TR sec
    # TR 3 -> 2*TR sec
    target_onset = (target_TR - 1) * TR
    # -------------------------------------------------------------------------
    # Find the condition that was active when the white fixation appeared
    condition_mask = ((BIDS_df["onset"] <= target_onset) & (target_onset < BIDS_df["onset"] + BIDS_df["duration"]))
    matching_rows = BIDS_df.index[condition_mask]

    if len(matching_rows) != 1:
        print(
            f"WARNING: fixation target at {target_onset:.3f}s "
            f"matched {len(matching_rows)} condition rows."
        )
        continue
    row_idx = matching_rows[0]

    # -------------------------------------------------------------------------
    # Mark target occurrence
    BIDS_df.loc[row_idx, "target_present"] = 1
    BIDS_df.loc[row_idx, "target_onset"] = target_onset

    # -------------------------------------------------------------------------
    # Find a valid button press after this target
    valid_responses = []
    for press_idx, press_time in enumerate(ButtonPressTimes):
        # Don't use the same press twice
        if press_idx in used_button_presses:
            continue
        if (target_onset <= press_time <target_onset + DETECTION_WINDOW):
            valid_responses.append((press_idx, press_time))

    # -------------------------------------------------------------------------
    # Score detection
    if len(valid_responses) > 0:
        # First valid response after target onset
        press_idx, first_response = min(valid_responses,key=lambda x: x[1])
        used_button_presses.add(press_idx)
        BIDS_df.loc[row_idx, "detected"] = 1
        BIDS_df.loc[row_idx, "response_time"] = (first_response - target_onset)
    else:
        BIDS_df.loc[row_idx, "detected"] = 0

# =============================================================================
# 4. SAVE
# =============================================================================
BIDS_dir = os.path.join('BIDS_events', expInfo['participant'],'func')

if not os.path.isdir(BIDS_dir):
    os.makedirs(BIDS_dir)

BIDS_output_file = os.path.join(BIDS_dir, f"{expInfo['participant']}_task-physical_"
    f"run-{int(expInfo['run']):02d}_events.tsv")

BIDS_df.to_csv(BIDS_output_file,sep="\t",index=False,float_format="%.3f")

print("\n================ BIDS EVENTS ================")
print(BIDS_df.to_string(index=False))
print(f"\nSAVED {BIDS_output_file}")

os.chdir(parentDir)
myWin.close()

# %% FINISH
core.quit()