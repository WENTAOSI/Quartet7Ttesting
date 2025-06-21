#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 11:52:06 2025
@author: siwentao
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 11:50:16 2025
@author: siwentao
"""
from psychopy import visual, event, core, monitors, logging, gui, data, misc
import numpy as np
import os
import random
import pandas as pd
import pickle
import sys
import time
import pylink
from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy
# %% SET PARAMS

###############################################################################
TRIGGERKEY = 'quoteleft'
# BLOCK DURATIONS [in TR]
# set durations of conditions and baseline
TR = 4.217     # sec in int if whole number  or float  

# set global offset
global_offset = (0,0)
 
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

# Store info about experiment and experimental run
expName = 'Volitional_MotQuart'  # set experiment name here
expInfo = {
    'run': '1',
    'participant': 'test',
    'Eyelink':['False','True'],
    'display': ['Vanderbilt7T', 'dbic'],
    'aspect_ratio': '1.12'
    }

# Create GUI at the beginning of exp to get more expInfo
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False: core.quit()  # user pressed cancel

# Circle properties
circle_dva = 6  # Diameter in dva
circle_radius = circle_dva / 2  # Radius in dva

# Define aspect ratio (Width / Height)
aspect_ratio = float(expInfo['aspect_ratio'])
aspect_ratio = 1/aspect_ratio  # Example: 2 means Width is twice the height

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
squareColor = np.multiply(backColor, -1)  # from -1 (black) to 1 (white)    Back dark grey; square light grey


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
# Name and create specific folder for output
outFolderName = dataFolderName + os.path.sep + 'Output'
if not os.path.isdir(outFolderName):
    os.makedirs(outFolderName)
outFileName = outFolderName + os.path.sep + '%s_%s_Run%s_%s' % (
expInfo['participant'], expInfo['expName'], expInfo['run'],
    expInfo['date'])

# Name and create specific folder for protocol files
prtFolderName = dataFolderName + os.path.sep + 'Protocols'
if not os.path.isdir(prtFolderName):
    os.makedirs(prtFolderName)
prtFileName = prtFolderName + os.path.sep + f'{expInfo["participant"]}_Volitional_Run{expInfo["run"]}_protocol'


# save a log file and set level for msg to be received
logFile = logging.LogFile(logFileName+'.log', level=logging.INFO)
logging.console.setLevel(logging.WARNING)  # set console to receive warningVEs

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
moni.setSizePix([PixW, PixH]) 

# log monitor info
logFile.write('MonitorDistance=' + str(distanceMon) + 'cm' + '\n')
logFile.write('MonitorWidth=' + str(widthMon) + 'cm' + '\n')
logFile.write('PixelWidth=' + str(PixW) + '\n')
logFile.write('PixelHeight=' + str(PixH) + '\n')

# set screen:
myWin = visual.Window(size=(PixW, PixH),
                      screen = 0,
                      winType='pyglet',  # winType : None, ‘pyglet’, ‘pygame’
                      allowGUI=False,
                      allowStencil=False,
                      fullscr=True,  # for psychoph lab: fullscr = True
                      monitor=moni,
                      color=backColor,
                      colorSpace='rgb',
                      units='deg',
                      blendMode='avg',
                      waitBlanking=True
                      )

# %% TRIAL DURATIONS SETUP
# Load the pickle file containing catch trial distribution
with open(os.path.join("Volitional_MotQuart","catch_trials_distribution.pkl"), "rb") as file:
    catch_trials_distribution = pickle.load(file)

current_run_catch_trials = catch_trials_distribution[int(expInfo['run']) - 1]

num_trials = 10

# Initialize parameters 
if TR == 1 or TR == 2: # ALL number is even 
    total_time = 30
    report = 6
    precue = 10
    delay = 10
    switch = 4

elif TR == 4.217:
    total_time = 8*TR
    report = 2*TR
    precue = 2*TR
    delay = 3*TR
    switch = 1*TR


elif isinstance(TR, float) and TR != 4.217:
    total_time = 15*TR
    report = 3*TR
    precue = 5*TR
    delay = 5*TR
    switch = 2*TR

# Create a balanced list of QuartetOrder values
quartet_orders = ["quartetPart1, quartetPart2"] * (num_trials // 2) + \
                 ["quartetPart2, quartetPart1"] * (num_trials // 2)

# Create balanced list of tone V OR H
instruct_V_H = ["vertical"] * (num_trials // 2) + ["horizontal"] * (num_trials // 2)

# Handle catch trials
physical_catch_trials = [trial for trial in current_run_catch_trials if trial in ["V", "H"]]

# Pair one "physical" trial with "vertical" and one with "horizontal"
physical_pairs = [("vertical", "physical") if trial == "V" else ("horizontal", "physical") for trial in physical_catch_trials]

# Remove "vertical" and "horizontal" trials from instruct_V_H to pair with physical trials
for pair in physical_pairs:
    instruct_V_H.remove(pair[0])

# Create the remaining illusory trials
illusory_physical = ["illusory"] * (num_trials - len(physical_pairs))
combined_trials = physical_pairs + list(zip(instruct_V_H, illusory_physical))

# Shuffle the combined list to randomize positions
random.shuffle(combined_trials)

# Shuffle other lists
random.shuffle(quartet_orders)
instruct_V_H = [trial[0] for trial in combined_trials]
illusory_physical = [trial[1] for trial in combined_trials]

# Initialize button press instructions
button_seq = [(1, 2), (1, 3), (1, 4), (2, 1), (2, 3), (2, 4), (3, 1), (3, 2), (3, 4), (4, 1), (4, 2), (4, 3)]
random.shuffle(button_seq)

# Generate conditions and timing
conditions = []
for trial in range(1, num_trials + 1):

    # Assign values for the current trial
    this_quartet_order = quartet_orders.pop()
    this_instruct_V_H = instruct_V_H.pop()
    this_illusory_physical = illusory_physical.pop()
    this_button_seq = button_seq.pop()

    # Get V/H key for the trial
    this_V = str(this_button_seq[0])
    this_H = str(this_button_seq[1])

    # Store the trial data
    conditions.append({
        "Trial": trial,
        "PrecueTime": precue,
        "DelayTime": delay,
        "SwitchTime": switch,
        "ReportTime": report,
        "QuartetOrder": this_quartet_order,
        "Instruct_V_H": this_instruct_V_H,
        "illusory_physical": this_illusory_physical,
        "V_buttom": this_V,
        "H_buttom": this_H
    })

# Convert to a DataFrame for visualization or saving
conditions_df = pd.DataFrame(conditions)

# %% STIMULI
# INITIALISE SOME STIMULI
SquareSize = 1.0  # 1.1 #1.8

logFile.write('SquareSize=' + str(SquareSize) + '\n')

dotFix = visual.Circle(myWin,
                       autoLog=False,
                       name='dotFix',
                       units='deg',
                       radius=0.1,
                       pos=apply_global_offset((0,0), global_offset),
                       fillColor='red',
                       lineColor='red'
                       )

Square = visual.GratingStim(myWin,
                            autoLog=False,
                            name='Square',
                            tex=None,
                            units='deg',
                            size=(SquareSize, SquareSize),
                            color= squareColor,
                            )
# Four Circles
circle_size = 1  # width of each circle
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
        myWin,
        autoLog=False,
        units='deg',
        radius=circle_size/2,  # Diameter of the circle
        pos=pos,  # Position from the list
        lineColor=squareColor,
        lineWidth=2,
        fillColor=None
    )
    circles.append(circle)

# Generate H/V letters in the circle
Hs = []
Vs = []

for pos in positions:
    H = visual.TextStim(
        win=myWin,
        color='white',
        height=circle_size-0.2,
        text='H',
        pos=pos
        )
    V = visual.TextStim(
        win=myWin,
        color='white',
        height=circle_size-0.2,
        text='V',
        pos=pos
        )
    Hs.append(H)
    Vs.append(V)

blue_Square = visual.GratingStim(myWin,
                            autoLog=False,
                            name='Square',
                            tex=None,
                            units='deg',
                            size=(SquareSize*1.5, SquareSize*1.5),
                            color='blue',
                            pos=apply_global_offset((0,4), global_offset)
                            )
blue_Square.name = "BLUE"

red_Square = visual.GratingStim(myWin,
                            autoLog=False,
                            name='Square',
                            tex=None,
                            units='deg',
                            size=(SquareSize*1.5, SquareSize*1.5),
                            color='red',
                            pos=apply_global_offset((0,4), global_offset)
                            )
red_Square.name = "RED"

triggerText = visual.TextStim(
    win=myWin,
    color='white',
    height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text='Experiment will start soon. Waiting for scanner'
    )

anykeyText = visual.TextStim(
    win=myWin,
    color='white',
    height=0.5,
    text='Press any key to continue',
    pos=apply_global_offset(base_pos=(0,-2), global_offset=global_offset)
    )

confirm_report_V = visual.TextStim(
    win=myWin,
    color='white',
    height=0.5,
    text='You have pressed VERTICAL',
    pos=apply_global_offset(base_pos=(0,-4), global_offset=global_offset)
    )

confirm_report_H = visual.TextStim(
    win=myWin,
    color='white',
    height=0.5,
    text='You have pressed HORIZONTAL',
    pos=apply_global_offset(base_pos=(0,-4), global_offset=global_offset)
    )

endText = visual.TextStim(
    win=myWin,
    color="white",
    height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text="Please rest until further instructions"
    )

# %% VOLITIONAL INSTRUCTION COLOR MAPPING 

# Define the blue and red mappings as functions
if int(expInfo['run'])<= 8:
    color_mapping = {
        "vertical": blue_Square,
        "horizontal": red_Square,
    }
    for condition in conditions:
        condition["vertical"] = "blue_Square"
        condition["horizontal"] = "red_Square"
    
else:
    color_mapping = {
    "vertical": red_Square,
    "horizontal": blue_Square,
    }
    for condition in conditions:
        condition["vertical"] = "red_Square"
        condition["horizontal"] = "blue_Square" 

mapping_instruct_v = visual.TextStim(
    win=myWin,
    color='white',
    height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text=(f'In this run, be prepared to see VERTICAL motion\n'
        f'after {color_mapping["vertical"].name} onset')
    )

mapping_instruct_h = visual.TextStim(
    win=myWin,
    color="white",
    height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text=(f'In this run, be prepared to see HORIZONTAL motion\n'
        f'after {color_mapping["horizontal"].name} onset')
    )

# %% TIME AND TIMING PARAMeTERS
# parameters
refr_rate = myWin.getActualFrameRate()  # get screen refresh rate

print(f"refr_rate{refr_rate}")
if refr_rate is None:
    refr_rate = 120.0 # if could not get reliable refresh rate

if refr_rate is not None:
    frameDur = 1.0/round(refr_rate)
else:
    frameDur = 1.0/round(refr_rate)  # couldn't get a reliable measure so guess

# physical quartet motion setup 
physical_duration = 0.5  # Total duration of the motion (seconds)
num_frames_physical = int(physical_duration * refr_rate)  # Number of frames for the motion
frame_interval = physical_duration / num_frames_physical  # Time per frame
phases = np.linspace(0, 1, num_frames_physical)  # Phase values from 0 to 1

logFile.write('RefreshRate=' + str(refr_rate) + '\n')
logFile.write('FrameDuration=' + str(frameDur) + '\n')

# define clock
clock = core.Clock()
logging.setDefaultClock(clock)

# %% FUNCTIONS
# create necessary functions for quartet
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

def quartetIntermedian(Hori, Verti, V_or_H):
    if V_or_H == "vertical":
        Square.setPos(apply_global_offset((0, Verti), global_offset))
        Square.draw()
        Square.setPos(apply_global_offset((0, -Verti), global_offset))
        Square.draw()
        dotFix.draw()
    elif V_or_H == "horizontal":
        Square.setPos(apply_global_offset((Hori, 0), global_offset))
        Square.draw()
        Square.setPos(apply_global_offset((-Hori, 0), global_offset))
        Square.draw()
        dotFix.draw()

# For catch trials physcial motion HORIZONTAL
def HMotion_update(Hori, Verti, sequence, progress=0):
    """
    Perform one-shot horizontal motion based on progress.
    :param Hori: Maximum horizontal distance.
    :param Verti: Fixed vertical distance.
    :param progress: Progress of the motion (0 to 1).
    :return: Horizontal position (mHori) of the square.
    """
    mHori = progress * Hori  # Move from center to maximum Hori
    
    if sequence == "quartetPart1, quartetPart2":
        Square.setPos(apply_global_offset((mHori, Verti), global_offset))  # Square northwest
        Square.draw()
        Square.setPos(apply_global_offset((-mHori, -Verti), global_offset))  # Square southeast
        Square.draw()
        dotFix.draw()
        myWin.flip()
        return mHori
    
    elif sequence == "quartetPart2, quartetPart1":
        Square.setPos(apply_global_offset((-mHori, Verti), global_offset))  # Square northeast
        Square.draw()
        Square.setPos(apply_global_offset((mHori, -Verti), global_offset))  # Square southwest
        Square.draw()
        dotFix.draw()
        myWin.flip()
        return mHori

def VMotion_update(Hori, Verti, sequence, progress=0):
    """
    Perform one-shot vertical motion based on progress.
    :param Hori: Fixed horizontal distance.
    :param Verti: Maximum vertical distance.
    :param progress: Progress of the motion (0 to 1).
    :return: Vertical position (mVerti) of the square.
    """
    mVerti = progress * Verti  # Move from center to maximum Verti
    
    if sequence == "quartetPart2, quartetPart1":
        Square.setPos(apply_global_offset((-Hori, mVerti), global_offset))  # Square northeast
        Square.draw()
        Square.setPos(apply_global_offset((Hori, -mVerti), global_offset))  # Square southwest
        Square.draw()
        dotFix.draw()
        myWin.flip()
        return mVerti
    
    elif sequence == "quartetPart1, quartetPart2":
        Square.setPos(apply_global_offset((Hori, mVerti), global_offset))  # Square northwest
        Square.draw()
        Square.setPos(apply_global_offset((-Hori, -mVerti), global_offset))  # Square southeast
        Square.draw()
        dotFix.draw()
        myWin.flip()
        return mVerti

def buttom_instruct(vertical_buttom, horizontal_buttom):
    '''
    Displays buttom instructions at the report stage 
    Parameters
    vertical_buttom: (str): "1","2","3",or"4"
    horizontal_buttom: (str): "1","2","3",or"4"

    Returns
    None
    '''
    # Draw and display the circles
    for circle in circles:
        circle.draw()
    Vs[int(vertical_buttom)-1].draw()   # Because python start counting from 0, draw the first one in the 0th on the list 
    Hs[int(horizontal_buttom)-1].draw() # 

def check_for_escape():
    """Check if the Escape key is pressed and exit the program."""
    keys = event.getKeys(keyList=['escape'])
    if 'escape' in keys:
        core.quit()

# %% EYELINK SETUP

# Set this variable to True if you use the built-in retina screen as your
# primary display device on macOS. If have an external monitor, set this
# variable True if you choose to "Optimize for Built-in Retina Display"
# in the Displays preference settings.
use_retina = False
# Set this variable to True to run the script in "Dummy Mode"
if expInfo['Eyelink'] == 'True':
    dummy_mode = False
elif expInfo['Eyelink'] == 'False':
    dummy_mode = True
# Set this variable to True to run the task in full screen mode
full_screen = True
# Set up EDF data file name and local data folder
# The EDF data filename should not exceed 8 alphanumeric characters
# use ONLY number 0-9, letters, & _ (underscore) in the filename
edf_fname = 'EYE'
# We download EDF data file from the EyeLink Host PC to the local hard
# drive at the end of each testing session, here we rename the EDF to
# include session start date/time
#time_str = time.strftime("_%Y_%m_%d_%H_%M", time.localtime())
session_identifier = edf_fname
# create a folder for the current testing session in the data folder for the subject 
session_folder = os.path.join(dataFolderName, session_identifier)
if not os.path.exists(session_folder):
    os.makedirs(session_folder)

print(f"session_folder = {session_folder}")
'''
----------------# Connect to the EyeLink Host PC-----------------------
'''
# The Host IP address, by default, is "100.1.1.1".
# the "el_tracker" objected created here can be accessed through the Pylink
# Set the Host PC address to "None" (without quotes) to run the script
# in "Dummy Mode"
if dummy_mode:
    el_tracker = pylink.EyeLink(None)
else:
    try:
        el_tracker = pylink.EyeLink("100.1.1.1")
    except RuntimeError as error:
        print('ERROR:', error)
        core.quit()
        sys.exit()

# Step 2: Open an EDF data file on the Host PC
edf_file =  f"vo_run{expInfo['run']}.EDF"
try:
    el_tracker.openDataFile(edf_file)
except RuntimeError as err:
    print('ERROR:', err)
    # close the link if we have one open
    if el_tracker.isConnected():
        el_tracker.close()
    core.quit()
    sys.exit()

# Add a header text to the EDF file to identify the current experiment name
# This is OPTIONAL. If your text starts with "RECORDED BY " it will be
# available in DataViewer's Inspector window by clicking
# the EDF session node in the top panel and looking for the "Recorded By:"
# field in the bottom panel of the Inspector.
preamble_text = 'RECORDED BY %s' % os.path.basename(__file__)
el_tracker.sendCommand("add_file_preamble_text '%s'" % preamble_text)

# Put the tracker in offline mode before we change tracking parameters
el_tracker.setOfflineMode()

# Get the software version:  1-EyeLink I, 2-EyeLink II, 3/4-EyeLink 1000,
# 5-EyeLink 1000 Plus, 6-Portable DUO
eyelink_ver = 5  # set version to 0, in case running in Dummy mode
if not dummy_mode:
    vstr = el_tracker.getTrackerVersionString()
    eyelink_ver = int(vstr.split()[-1].split('.')[0])
    # print out some version info in the shell
    print('Running experiment on %s, version %d' % (vstr, eyelink_ver))

'''
----------------# Setting up for CALIBRATION and other stuff-----------------------
'''
# File and Link data control
# what eye events to save in the EDF file, include everything by default
file_event_flags = 'LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT'
# what eye events to make available over the link, include everything by default
link_event_flags = 'LEFT,RIGHT,FIXATION,SACCADE,BLINK,BUTTON,FIXUPDATE,INPUT'
# what sample data to save in the EDF data file and to make available
# over the link, include the 'HTARGET' flag to save head target sticker
# data for supported eye trackers
if eyelink_ver > 3:
    file_sample_flags = 'LEFT,RIGHT,GAZE,HREF,RAW,AREA,HTARGET,GAZERES,BUTTON,STATUS,INPUT'
    link_sample_flags = 'LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT'
else:
    file_sample_flags = 'LEFT,RIGHT,GAZE,HREF,RAW,AREA,GAZERES,BUTTON,STATUS,INPUT'
    link_sample_flags = 'LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT'
el_tracker.sendCommand("file_event_filter = %s" % file_event_flags)
el_tracker.sendCommand("file_sample_data = %s" % file_sample_flags)
el_tracker.sendCommand("link_event_filter = %s" % link_event_flags)
el_tracker.sendCommand("link_sample_data = %s" % link_sample_flags)

# Optional tracking parameters
# Sample rate, 250, 500, 1000, or 2000, check your tracker specification
# if eyelink_ver > 2:
#     el_tracker.sendCommand("sample_rate 1000")
# Choose a calibration type, H3, HV3, HV5, HV13 (HV = horizontal/vertical),
el_tracker.sendCommand("calibration_type = HV9")

# Set a gamepad button to accept calibration/drift check target
# You need a supported gamepad/button box that is connected to the Host PC
el_tracker.sendCommand("button_function 5 'accept_target_fixation'")

# Get screen resolution in pixels
scn_width, scn_height = myWin.size

# EyeLink coordinate setup
el_coords = "screen_pixel_coords = 0 0 %d %d" % (scn_width - 1, scn_height - 1)
el_tracker.sendCommand(el_coords)

# Data Viewer setup
dv_coords = "DISPLAY_COORDS  0 0 %d %d" % (scn_width - 1, scn_height - 1)
el_tracker.sendMessage(dv_coords)

os.chdir(_thisDir) # make sure that EyeLinkCoreGraphicsPsychoPy can be read
# Graphics environment setup
#sound.AudioDeviceInfo(deviceIndex = -1, deviceName =u'Microsoft Sound Mapper - Output', inputChannels = 0, outputLatency = (0.09, 0.18), inputLatency = (0.09, 0.18), defaultSampleRate = 44100)
genv = EyeLinkCoreGraphicsPsychoPy(el_tracker, myWin)
print(genv)  # print out the version number of the CoreGraphics library

# Set background and foreground colors for the calibration target
# in PsychoPy, (-1, -1, -1)=black, (1, 1, 1)=white, (0, 0, 0)=mid-gray
foreground_color = (-1, -1, -1)
background_color = myWin.color
genv.setCalibrationColors(foreground_color, background_color)

# Set up the calibration target
# The target could be a "circle" (default), a "picture", a "movie" clip,
# or a rotating "spiral". To configure the type of calibration target, set
# genv.setTargetType to "circle", "picture", "movie", or "spiral", e.g.,
# genv.setTargetType('picture')
#
# Use gen.setPictureTarget() to set a "picture" target
# genv.setPictureTarget(os.path.join('images', 'fixTarget.bmp'))
#
# Use genv.setMovieTarget() to set a "movie" target
# genv.setMovieTarget(os.path.join('videos', 'calibVid.mov'))

# Use a picture as the calibration target
genv.setTargetType('picture')

genv.setPictureTarget(os.path.join(_thisDir,'images', 'fixTarget.bmp'))

# Configure the size of the calibration target (in pixels)
# this option applies only to "circle" and "spiral" targets
# genv.setTargetSize(24)

# Beeps to play during calibration, validation and drift correction
# parameters: target, good, error
#     target -- sound to play when target moves
#     good -- sound to play on successful operation
#     error -- sound to play on failure or interruption
# Each parameter could be ''--default sound, 'off'--no sound, or a wav file
genv.setCalibrationSounds('', '', '')

# resolution fix for macOS retina display issues
if use_retina:
    genv.fixMacRetinaDisplay()

# Request Pylink to use the PsychoPy window we opened above for calibration
pylink.openGraphicsEx(genv)

os.chdir(parentDir) # change dirctory back
'''
 -------------# define a few helper functions for EYELINK handling ---------------------
'''
def clear_screen(win):
    """ clear up the PsychoPy window"""
    win.fillColor = genv.getBackgroundColor()
    win.flip()
    
def show_msg(win, text, wait_for_keypress=True):
    """ Show task instructions on screen"""
    msg = visual.TextStim(win, text,
                          color=genv.getForegroundColor(),
                          wrapWidth=scn_width/2)
    clear_screen(win)
    msg.draw()
    win.flip()
    # wait indefinitely, terminates upon any key press
    if wait_for_keypress:
        event.waitKeys()
        clear_screen(win)

def abort_trial():
    """Ends recording """
    el_tracker = pylink.getEYELINK()
    # Stop recording
    if el_tracker.isRecording():
        # add 100 ms to catch final trial events
        pylink.pumpDelay(100)
        el_tracker.stopRecording()
    # clear the screen
    clear_screen(myWin)
    # Send a message to clear the Data Viewer screen
    bgcolor_RGB = (116, 116, 116)
    el_tracker.sendMessage('!V CLEAR %d %d %d' % bgcolor_RGB)
    # send a message to mark trial end
    el_tracker.sendMessage('TRIAL_RESULT %d' % pylink.TRIAL_ERROR)
    return pylink.TRIAL_ERROR

def terminate_task():
    """ Terminate the task gracefully and retrieve the EDF data file

    file_to_retrieve: The EDF on the Host that we would like to download
    win: the current window used by the experimental script
    """
    el_tracker = pylink.getEYELINK()
    if el_tracker.isConnected():
        # Terminate the current trial first if the task terminated prematurely
        error = el_tracker.isRecording()
        if error == pylink.TRIAL_OK:
            abort_trial()
        # Put tracker in Offline mode
        el_tracker.setOfflineMode()
        # Clear the Host PC screen and wait for 500 ms
        el_tracker.sendCommand('clear_screen 0')
        pylink.msecDelay(500)
        # Close the edf data file on the Host
        el_tracker.closeDataFile()
        # Show a file transfer message on the screen
        #msg = 'EDF data is transferring from EyeLink Host PC...'
        #show_msg(myWin, msg, wait_for_keypress=False)
        # Download the EDF data file from the Host PC to a local data folder
        # parameters: source_file_on_the_host, destination_file_on_local_drive
        local_edf = os.path.join(session_folder, edf_file)
        print("LOCAL_EDF = " + local_edf)
        try:
            el_tracker.receiveDataFile(edf_file, local_edf)
        except RuntimeError as error:
            print('ERROR:', error)
        # Close the link to the tracker.
        el_tracker.close()
    # close the PsychoPy window
    myWin.close()
    # quit PsychoPy
    core.quit()
    sys.exit()

'''
---------------# Set up the camera and calibrate the tracker ---------------------
'''
# Show the task instructions
task_msg = 'In the task, you may press the SPACEBAR to end a trial\n' + \
    '\nPress Ctrl-C to if you need to quit the task early\n'
if dummy_mode:
    task_msg = task_msg + '\nNow, press ENTER to start the task'
else:
    task_msg = task_msg + '\nNow, press ENTER twice to calibrate tracker'
show_msg(myWin, task_msg)

# Calibration 
# skip this step if running the script in Dummy Mode
if not dummy_mode:
    try:
        el_tracker.doTrackerSetup()        
    except RuntimeError as err:
        print('ERROR:', err)
        el_tracker.exitCalibration()

# put tracker in idle/offline mode before recording
el_tracker.setOfflineMode()
# Start recording
# arguments: sample_to_file, events_to_file, sample_over_link,
# event_over_link (1-yes, 0-no)
try:
    el_tracker.startRecording(1, 1, 1, 1)
except RuntimeError as error:
    print("ERROR:", error)
    abort_trial()
        
# %% RENDER_LOOp
# Mapping instruction 
mapping_instruct_v.draw()
color_mapping["vertical"].draw()
anykeyText.draw()
myWin.flip()
# Wait for any key press to continue
event.waitKeys()

mapping_instruct_h.draw()
color_mapping["horizontal"].draw()
anykeyText.draw()
myWin.flip()
# Wait for any key press to continue
event.waitKeys()

# wait for scanner trigger
triggerText.draw()
myWin.flip()
event.waitKeys(keyList=[TRIGGERKEY], timeStamped=False)
# reset clocks
clock.reset()
logFile.write('StartOfRun' + str(expInfo['run']))

el_tracker.sendMessage(f"EXPERIMENT_START {expName}")

num_trial = 0

# Main trial loop 
for trial in conditions:
     
    num_trial += 1
    logFile.write(f'Time at start of trial {num_trial} is {clock.getTime()}\n')
    
    # PRECUE PERIOD ----------------------------------------------------------------------------------------------------------
    PrecueDur = trial["PrecueTime"]  # Duration in seconds
    if num_trial == 1:
        precue_start_time = clock.getTime()  # Get the current time at the start of the precue period
    else:
        precue_start_time = report_start + ReportDur
    # Signal tracker
    el_tracker.sendMessage(f"Trial: {num_trial} ; PRECUE START")
    
    # Run a while loop until the total PrecueDur is completed
    while clock.getTime() - precue_start_time < PrecueDur:
        
        check_for_escape()
        
        if trial["QuartetOrder"] == "quartetPart1, quartetPart2":
            quartetPart1(HoriDist, VertiDist)  # Draw quartetPart1
        elif trial["QuartetOrder"] == "quartetPart2, quartetPart1":
            quartetPart2(HoriDist, VertiDist)  # Draw quartetPart2
        
        myWin.flip()  # Flip the window to update the display
    
    # DELAY PERIOD -----------------------------------------------------------------------------------------------------------
    DelayDur = trial["DelayTime"]  # Delay duration in seconds
    
    delay_start_time = precue_start_time + PrecueDur  # Record the start time of the delay period
    
    el_tracker.sendMessage(f"Trial: {num_trial} ; DELAY START")
    # Run a while loop until the total DelayDur is completed
    while clock.getTime() - delay_start_time < DelayDur:
        
        check_for_escape()
        
        # Maintain the same visual display as the precue period
        if trial["QuartetOrder"] == "quartetPart1, quartetPart2":
            while clock.getTime() - delay_start_time < 2:
                
                check_for_escape()
                
                # Dynamically call the function based on the mapping
                color_mapping[trial["Instruct_V_H"]].draw()
                # Draw quartetPart1
                quartetPart1(HoriDist, VertiDist)  
                myWin.flip()
            # Draw quartetPart1    
            quartetPart1(HoriDist, VertiDist)  
            
        elif trial["QuartetOrder"] == "quartetPart2, quartetPart1":
            while clock.getTime() - delay_start_time < 2:
                
                check_for_escape()
                
                # Dynamically call the function based on the mapping
                color_mapping[trial["Instruct_V_H"]].draw()
                # Draw quartetPart2
                quartetPart2(HoriDist, VertiDist)  
                myWin.flip()
            # Draw quartetPart2    
            quartetPart2(HoriDist, VertiDist)   
        # Flip the window to display visuals
        myWin.flip()
 
    # SWITCH PERIOD ----------------------------------------------------------------------------------------------------------
    # Get the duration of the SWITCH period in frames
    SwitchDur = trial["SwitchTime"]
    # Initialize variables to track invalid responses
    invalid_key = None
    invalid_timestamp = None
    switch_start_time = delay_start_time + DelayDur
    
    el_tracker.sendMessage(f"Trial: {num_trial} ; SWITCH START")
    
    # Illusory Motion: Static switch between parts
    if trial["illusory_physical"] == "illusory":
        
        while clock.getTime() - switch_start_time < SwitchDur:
            
            check_for_escape()
            
            # Check for button presses during each frame
            keys = event.getKeys(keyList=['1', '2', '3', '4'], timeStamped=clock)
            if keys and invalid_key is None:  # Record only the first invalid key press
                invalid_key, invalid_timestamp = keys[0]
    
            # Display static visual based on QuartetOrder
            if trial["QuartetOrder"] == "quartetPart1, quartetPart2":
                quartetPart2(HoriDist, VertiDist)
            elif trial["QuartetOrder"] == "quartetPart2, quartetPart1":
                quartetPart1(HoriDist, VertiDist)
            
            myWin.flip()  # Flip window every frame
    
    # Physical Motion: Motion followed by static display
    elif trial["illusory_physical"] == "physical":
        # Get the sequence (QuartetOrder)
        sequence = trial["QuartetOrder"]
    
        if trial["Instruct_V_H"] == "vertical":
            # Vertical motion across frames
            physical_start_time = clock.getTime()
              
            # IF USING INTERMEDIATE STEP TO FACILITATE V/H
            quartetIntermedian(HoriDist, VertiDist, trial["Instruct_V_H"])
            core.wait(.25) # let it be short
            # record physical end time 
            physical_end_time = clock.getTime()
            
            # Static display after motion (remaining time)
            while clock.getTime() - physical_end_time <  SwitchDur - (physical_end_time - physical_start_time):
                
                check_for_escape()
                
                # Check for button presses during each frame
                keys = event.getKeys(keyList=['1', '2', '3', '4'], timeStamped=clock)
                if keys and invalid_key is None:  # Record only the first invalid key press
                    invalid_key, invalid_timestamp = keys[0]
    
                if sequence == "quartetPart1, quartetPart2":
                    quartetPart2(HoriDist, VertiDist)
                elif sequence == "quartetPart2, quartetPart1":
                    quartetPart1(HoriDist, VertiDist)
                myWin.flip()  # Flip the screen for static frames
    
        elif trial["Instruct_V_H"] == "horizontal":
            
            # Horizontal motion across frames
            physical_start_time = clock.getTime()
            
            # IF USING INTERMEDIATE STEP TO FACILITATE V/H
            quartetIntermedian(HoriDist, VertiDist, trial["Instruct_V_H"])
            core.wait(.25) # let it be short         
            # record physical end time    
            physical_end_time = clock.getTime()
            
            # Static display after motion (remaining time)
            while clock.getTime() - physical_end_time <  SwitchDur - (physical_end_time - physical_start_time):
                
                check_for_escape()
            
                # Check for button presses during each frame
                keys = event.getKeys(keyList=['1', '2', '3', '4'], timeStamped=clock)
                if keys and invalid_key is None:  # Record only the first invalid key press
                    invalid_key, invalid_timestamp = keys[0]
                    
                if sequence == "quartetPart1, quartetPart2":
                    quartetPart2(HoriDist, VertiDist)
                elif sequence == "quartetPart2, quartetPart1":
                    quartetPart1(HoriDist, VertiDist)
                myWin.flip()  # Flip the screen for static frames
    
    # Record invalid response
    if invalid_key is not None:
        # Save invalid button press to trial data
        trial["invalid_ResponseKey"] = invalid_key
        trial["invalid_ResponseTime"] = invalid_timestamp
    else:
        trial["invalid_ResponseKey"] = "None"
        trial["invalid_ResponseTime"] = "None"
    
    # REPORT PERIOD ----------------------------------------------------------------------------------------------------------
    ReportDur = trial["ReportTime"]  # Report duration in seconds
    #logFile.write(f"ReportDur == {ReportDur} sec\n")
    
    V_buttom = trial["V_buttom"]
    H_buttom = trial["H_buttom"]
    
    response_recorded = False  # Track if a response is recorded
    response_key = None
    response_time = None
    
    # Start the report period
    report_start = switch_start_time + SwitchDur  # Record the start time
    
    el_tracker.sendMessage(f"Trial: {num_trial} ; REPORT START")
    
    while clock.getTime() - report_start < ReportDur:  # Loop for the duration of the report period
        # Draw the report text and button image
        buttom_instruct(V_buttom, H_buttom)
        
        # Check for button presses
        keys = event.getKeys(keyList=['1', '2', '3', '4'], timeStamped=clock)
        if keys and not response_recorded:  # Process only the first response
            response_key, response_time = keys[0]  # Extract the key and timestamp
            response_recorded = True  # Mark the response as recorded
    
            # Record response
            trial["ResponseKey"] = response_key
            trial["ResponseTime"] = response_time
    
            # Draw confirmation text based on the response
            if response_key == V_buttom:
                confirm_report_V.draw()  # Vertical confirmation text
            elif response_key == H_buttom:
                confirm_report_H.draw()  # Horizontal confirmation text
        
        if response_recorded:
            # Keep confirmation text 
            if response_key == V_buttom:
                confirm_report_V.draw()  # Vertical confirmation text
            elif response_key == H_buttom:
                confirm_report_H.draw()  # Horizontal confirmation tex
    
        # Flip the window to display stimuli
        myWin.flip()
    
    logFile.write(f'Time at the end of trial {num_trial} is {clock.getTime()}\n')
    
    # Handle timeout case (no response)
    if not response_recorded:
        trial["ResponseKey"] = "None"
        trial["ResponseTime"] = "None"
    
    # Log response
    logFile.write(f"Trial {trial['Trial']} INVALID Response: {trial['invalid_ResponseKey']} at {trial['invalid_ResponseTime']} sec\n")
    logFile.write(f"Trial {trial['Trial']} Response: {trial['ResponseKey']} at {trial['ResponseTime']} sec\n")

# SIGNAL TRACKER EXPERIMENT END
el_tracker.sendMessage("EXPERIMENT_END")
    
# Convert conditions to a DataFrame
conditions_df = pd.DataFrame(conditions)

# Adding expected press and success True/False
conditions_df['expected_key'] = conditions_df.apply(
    lambda row: row['H_buttom'] if (row['Instruct_V_H'] == 'horizontal' and row['illusory_physical'] == 'illusory') else
                row['V_buttom'] if (row['Instruct_V_H'] == 'vertical' and row['illusory_physical'] == 'illusory') else
                row['V_buttom'] if (row['Instruct_V_H'] == 'horizontal' and row['illusory_physical'] == 'physical') else
                row['H_buttom'], axis=1
)
conditions_df['success'] = conditions_df.apply(
    lambda row: True if (row['illusory_physical'] == 'illusory' and row['ResponseKey'] == row['expected_key']) else '', axis=1
)

# Save responses DataFrame to the Output folder as a CSV file
conditions_df.to_csv(outFileName + '.csv', index=False)

# Log the saving process 
logFile.write(f"Responses saved to {outFileName}.csv")

# Construct protocol file save into protocol folder 
protocol_df = conditions_df.set_index([col for col in conditions_df.columns if col not in ['PrecueTime', 'DelayTime', 'SwitchTime', 'ReportTime']])
protocol_df = protocol_df.stack().reset_index()
protocol_df.columns = [*protocol_df.columns[:-2], 'Condition', 'Duration']
# Reorder columns to place Condition and Duration as the second and third columns
cols = list(protocol_df.columns)
cols.insert(1, cols.pop(cols.index('Condition')))
cols.insert(2, cols.pop(cols.index('Duration')))
protocol_df = protocol_df[cols]
# Duration should be modifed as TRs by dividing the time/TR
protocol_df['Duration'] = (protocol_df['Duration'] / TR).astype(int)
# Create Timestamp column with cumulative time per Trial
protocol_df['Timestamp'] = protocol_df['Duration'].cumsum()
# Create Onset time column
protocol_df['Onset'] = protocol_df['Timestamp'] - protocol_df['Duration']
# change 'Condition' column label to 'Stim'
protocol_df.rename(columns={'Condition': 'Stim'}, inplace=True)
# Save protocol DataFrame to the protocol folder as a CSV file
protocol_df.to_csv(prtFileName + '.csv', index=False)

# End of experiment 
endText.draw()
myWin.flip()
core.wait(5) # wait for 5 sec

# EYETRACKER CLOSE DISPLAY AND SAVE EDF
os.chdir(parentDir)
el_tracker.stopRecording()
terminate_task()
myWin.close()
core.quit()    

