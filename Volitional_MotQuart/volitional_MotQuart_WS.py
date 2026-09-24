#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 11:52:06 2025
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
from pathlib import Path
# %% SET PARAMS
# GUI Store info about experiment and experimental run
expName = 'Volitional_MotQuart'  # set experiment name here
expInfo = {
    'run': '1',
    'participant': 'test',
    'display': ['Vanderbilt7T', 'dbic'],
    'aspect_ratio': '1.12',
    'ho_dva': '-0.0981',
    'vo_dva': '1.7652',
    'TR': '2',
    }
# Create GUI at the beginning of exp to get more expInfo
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName, sortKeys=False)
if dlg.OK == False: core.quit()  # user pressed cancel
TRIGGERKEY = 'quoteleft'
# set global offset
ho_dva = expInfo['ho_dva']; vo_dva = expInfo['vo_dva']
global_offset = (float(ho_dva), float(-float(vo_dva)))
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
TR = float(expInfo['TR'])
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
VertiDist = max_height / 2; HoriDist = max_width / 2
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
# BIDS output directory
BIDSoutput_dir = os.path.join('BIDS_events', f"sub-{expInfo['participant']}", 'func')
if not os.path.isdir(BIDSoutput_dir):
    os.makedirs(BIDSoutput_dir)
# Name and create data folder for the experiment
dataFolderName = subjFolderName + os.path.sep + '%s' % (expInfo['expName'])
if not os.path.isdir(dataFolderName):
    os.makedirs(dataFolderName)
# Name and create specific folder for logging results
logFolderName = dataFolderName + os.path.sep + 'Logging'
if not os.path.isdir(logFolderName):
    os.makedirs(logFolderName)
logFileName = logFolderName + os.path.sep + '%s_%s_Run%s_%s' % (expInfo['participant'], expInfo['expName'], expInfo['run'],expInfo['date'])
# Name and create specific folder for output
outFolderName = dataFolderName + os.path.sep + 'Output'
if not os.path.isdir(outFolderName):
    os.makedirs(outFolderName)
outFileName = outFolderName + os.path.sep + '%s_%s_Run%s_%s' % (expInfo['participant'], expInfo['expName'], expInfo['run'], expInfo['date'])
# Name and create specific folder for protocol files
prtFolderName = dataFolderName + os.path.sep + 'Protocols'
if not os.path.isdir(prtFolderName):
    os.makedirs(prtFolderName)
prtFileName = prtFolderName + os.path.sep + f'{expInfo["participant"]}_Volitional_Run{expInfo["run"]}_protocol'
# save a log file and set level for msg to be received
logFile = logging.LogFile(logFileName+'.log', level=logging.INFO)
logging.console.setLevel(logging.WARNING)  # set console to receive warningVEs
# %% MONITOR AND WINDOW
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
if expInfo['display'] == 'Vanderbilt7T':
    screen=1
elif expInfo['display'] == 'dbic':
    screen=0
myWin = visual.Window(size=(PixW, PixH), screen = screen, winType='pyglet', allowGUI=False, allowStencil=False,fullscr=True, 
                      monitor=moni, color=backColor, colorSpace='rgb', units='deg', blendMode='avg', waitBlanking=True)
# %% TRIAL DURATIONS SETUP
use_catch_trials = True
# %% TRIAL DURATIONS SETUP
num_trials = 12
report = 2
precue = 3
delay = [3,4,5]
switch = 1
use_catch_trials = True

# Counterbalanced REAL trials: 2 V + 2 H at each delay
timing_instruct = [
    (3,"vertical"),(3,"vertical"),(3,"horizontal"),(3,"horizontal"),
    (4,"vertical"),(4,"vertical"),(4,"horizontal"),(4,"horizontal"),
    (5,"vertical"),(5,"vertical"),(5,"horizontal"),(5,"horizontal")
]

button_seq = [(1,2),(1,3),(1,4),(2,1),(2,3),(2,4),(3,1),(3,2),(3,4),(4,1),(4,2),(4,3)]
random.shuffle(button_seq)

# Generate 12 counterbalanced real trials
conditions = []
for this_delay,this_instruct in timing_instruct:
    this_button = button_seq.pop()
    conditions.append({
        "PrecueTime":precue,
        "DelayTime":this_delay,
        "SwitchTime":switch,
        "ReportTime":report,
        "QuartetOrder":"quartetPart1, quartetPart2",
        "Instruct_V_H":this_instruct,
        "illusory_physical":"illusory",
        "V_buttom":str(this_button[0]),
        "H_buttom":str(this_button[1])
    })

# Add ONE additional catch trial
if use_catch_trials:
    catch_instruct = random.choice(["vertical","horizontal"])
    catch_button = random.choice([(1,2),(1,3),(1,4),(2,1),(2,3),(2,4),(3,1),(3,2),(3,4),(4,1),(4,2),(4,3)])

    conditions.append({
        "PrecueTime":precue,
        "DelayTime":3,
        "SwitchTime":switch,
        "ReportTime":report,
        "QuartetOrder":"quartetPart1, quartetPart2",
        "Instruct_V_H":catch_instruct,
        "illusory_physical":"physical",
        "V_buttom":str(catch_button[0]),
        "H_buttom":str(catch_button[1])
    })

# Shuffle real + catch trials together
random.shuffle(conditions)

# Assign trial numbers AFTER shuffling
for trial,c in enumerate(conditions,start=1):
    c["Trial"] = trial

conditions = [{"Trial":c.pop("Trial"),**c} for c in conditions]
conditions_df = pd.DataFrame(conditions)

# Print
for c in conditions:
    print(f'Trial {c["Trial"]}: Instruct {c["Instruct_V_H"]}, Type {c["illusory_physical"]}, Delay {c["DelayTime"]}')

# %% STIMULI# %% STIMULI
# INITIALISE SOME STIMULI
SquareSize = 1.0  # 1.1 #1.8
logFile.write('SquareSize=' + str(SquareSize) + '\n')
dotFix = visual.Circle(myWin, autoLog=False, name='dotFix', units='deg',radius=0.1, pos=apply_global_offset((0,0), global_offset), fillColor='white', lineColor='white' )
# fixation changes color to blue or red to indicate cue (initialize it here )
dotFix_blue = visual.Circle(myWin, autoLog=False, name='dotFix_blue', units='deg',radius=0.1, pos=apply_global_offset((0,0), global_offset), fillColor='blue', lineColor='blue' )
dotFix_blue.name = "BLUE"
dotFix_red = visual.Circle(myWin, autoLog=False, name='dotFix_red', units='deg',radius=0.1, pos=apply_global_offset((0,0), global_offset), fillColor='red', lineColor='red' )
dotFix_red.name = "RED"
# quartet square 
Square = visual.GratingStim(myWin, autoLog=False, name='Square', tex=None, units='deg', size=(SquareSize, SquareSize), color= squareColor)
# Four Circles
circle_size = 1  # width of each circle
if expInfo["display"] == 'dbic':
    positions = [apply_global_offset((-3.5, 4), global_offset), apply_global_offset((-2, 4), global_offset),\
                 apply_global_offset((2, 4), global_offset), apply_global_offset((3.5, 4), global_offset)]  # Anchored positions 1, 2, 3, 4
elif expInfo["display"] == 'Vanderbilt7T':
    positions = [apply_global_offset((-2.5, 0), global_offset), apply_global_offset((-1, 0.5), global_offset),\
                 apply_global_offset((1, 0.5), global_offset), apply_global_offset((2.5, 0), global_offset)]  # Anchored positions 1, 2, 3, 4
# Generate circle objects at the specified positions
circles = []
for pos in positions:
    # Create a circle object at each specified position
    circle = visual.Circle(myWin, autoLog=False, units='deg', radius=circle_size/2, pos=pos, lineColor=squareColor, lineWidth=2, fillColor=None)
    circles.append(circle)
# Generate H/V letters in the circle
Hs = []; Vs = []
for pos in positions:
    H = visual.TextStim(win=myWin, color='white', height=circle_size-0.2,text='H', pos=pos)
    V = visual.TextStim(win=myWin, color='white', height=circle_size-0.2,text='V', pos=pos)
    Hs.append(H); Vs.append(V)
# generate text
triggerText = visual.TextStim(
    win=myWin, color='white', height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text='Experiment will start soon. Waiting for scanner'
    )
anykeyText = visual.TextStim(
    win=myWin, color='white', height=0.5,
    text='Press any key to continue',
    pos=apply_global_offset(base_pos=(0,-2), global_offset=global_offset)
    )
confirm_report_V = visual.TextStim(
    win=myWin, color='white', height=0.5,
    text='You have pressed VERTICAL',
    pos=apply_global_offset(base_pos=(0,-1), global_offset=global_offset)
    )
confirm_report_H = visual.TextStim(
    win=myWin, color='white', height=0.5,
    text='You have pressed HORIZONTAL',
    pos=apply_global_offset(base_pos=(0,-1), global_offset=global_offset)
    )
endText = visual.TextStim(
    win=myWin, color="white", height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text="Please rest until further instructions"
    )
# %% VOLITIONAL INSTRUCTION COLOR MAPPING 
# Define the blue and red mappings as functions
run_number = int(''.join(filter(str.isdigit, expInfo['run'])))
print(run_number)
if run_number == 1 or run_number == 2 or run_number == 5 or run_number == 6 or run_number == 9 or run_number == 10:
    color_mapping = {"vertical": dotFix_blue, "horizontal": dotFix_red}
    for condition in conditions:
        condition["vertical"] = "dotFix_blue"
        condition["horizontal"] = "dotFix_red"
elif run_number == 3 or run_number == 4 or run_number == 7 or run_number == 8 or run_number == 11 or run_number == 12:
    color_mapping = { "vertical": dotFix_red, "horizontal": dotFix_blue}
    for condition in conditions:
        condition["vertical"] = "dotFix_red"
        condition["horizontal"] = "dotFix_blue" 
        
mapping_instruct_v = visual.TextStim(
    win=myWin,color='white', height=0.5,
    pos=apply_global_offset(base_pos=(0,-1), global_offset=global_offset),
    text=(f'In this run, be prepared to see VERTICAL motion\n' f'after {color_mapping["vertical"].name} onset')
    )
mapping_instruct_h = visual.TextStim(
    win=myWin, color="white", height=0.5,
    pos=apply_global_offset(base_pos=(0,-1), global_offset=global_offset),
    text=(f'In this run, be prepared to see HORIZONTAL motion\n' f'after {color_mapping["horizontal"].name} onset')
    )
# %% TIME AND TIMING PARAMeTERS
# parameters
refr_rate = myWin.getActualFrameRate()  # get screen refresh rate
print(f"refr_rate{refr_rate}")
if refr_rate is None:
    refr_rate = 60.0 # if could not get reliable refresh rate
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
    Hs[int(horizontal_buttom)-1].draw()  
    dotFix.draw()
    
def check_for_escape():
    keys = event.getKeys(keyList=['escape'])
    if 'escape' in keys:
        core.quit()
        
tr_count = 0
last_trigger_time = None
def check_TR_trigger():
    """
    Check for scanner trigger without blocking. Call this every frame.
    """
    global tr_count, last_trigger_time
    keys = event.getKeys(keyList=[TRIGGERKEY],timeStamped=clock)
    for key, timestamp in keys:
        tr_count += 1
        last_trigger_time = timestamp
        logFile.write(f"TR {tr_count}: {timestamp:.6f} sec\n")
    return tr_count

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
# start a test clock
test_clock = core.Clock()
# reset clocks
clock.reset()
logFile.write('StartOfRun' + str(expInfo['run']))
#=============================================================
# initial fixation
initial_fix_TRs = final_fix_TR =  4
while tr_count < initial_fix_TRs:
    check_TR_trigger()
    dotFix.draw()
    myWin.flip()

num_trial = 0
print(conditions)
# Main trial loop 
for trial in conditions:
    num_trial += 1
    logFile.write(f'Time at start of trial {num_trial} is {clock.getTime()}\n')
    # Establish TR boundaries for this trial
    trial_start_TR = tr_count
    precue_end_TR = trial_start_TR + trial['PrecueTime']
    delay_end_TR  = precue_end_TR + trial['DelayTime']
    switch_end_TR = delay_end_TR + trial['SwitchTime']
    report_end_TR = switch_end_TR + trial['ReportTime']
    logFile.write(
        f"Trial {num_trial}: "
        f"start TR={trial_start_TR}, "
        f"precue end={precue_end_TR}, "
        f"delay end={delay_end_TR}, "
        f"switch end={switch_end_TR}, "
        f"report end={report_end_TR}\n"
    )
    #INITIALIZE TRIAL RESPONSE VARIABLES
    trial["invalid_ResponseKey"] = "None"; trial["invalid_ResponseTime"] = "None"; trial["ResponseKey"] = "None"
    trial["ResponseTime"] = "None"; trial["ResponseRT"] = "None"
    #========================================================
    #PRECUE
    while tr_count < precue_end_TR:
        check_for_escape()
        check_TR_trigger()
        if trial["QuartetOrder"] == "quartetPart1, quartetPart2":
            quartetPart1(HoriDist, VertiDist)
        elif trial["QuartetOrder"] == "quartetPart2, quartetPart1":
            quartetPart2(HoriDist, VertiDist)
        myWin.flip()
    #========================================================
    #DELAY
    delay_start_clock = clock.getTime()
    while tr_count < delay_end_TR:
        check_for_escape()
        check_TR_trigger()
        if trial["QuartetOrder"] == "quartetPart1, quartetPart2":
            quartetPart1(HoriDist, VertiDist)
        elif trial["QuartetOrder"] == "quartetPart2, quartetPart1":
            quartetPart2(HoriDist, VertiDist)
        # Instruction cue for first 1 seconds
        show_instruction = (clock.getTime() - delay_start_clock < 1)
        if show_instruction:
            color_mapping[trial["Instruct_V_H"]].draw()
        myWin.flip()
    #========================================================
    #SWITCH
    invalid_key = None
    invalid_timestamp = None

    if trial["illusory_physical"]=="physical":
        quartetIntermedian(HoriDist,VertiDist, trial["Instruct_V_H"])
        myWin.flip()

    while tr_count < switch_end_TR:
        check_for_escape()
        check_TR_trigger()
        keys = event.getKeys(keyList=['1','2','3','4'],timeStamped=clock)
        if keys and invalid_key is None:
            invalid_key,invalid_timestamp = keys[0]

        if trial["QuartetOrder"]=="quartetPart1, quartetPart2":
            quartetPart2(HoriDist,VertiDist)
        elif trial["QuartetOrder"]=="quartetPart2, quartetPart1":
            quartetPart1(HoriDist,VertiDist)

        myWin.flip()
    # ========================================================
    # REPORT
    ReportDur = trial["ReportTime"]
    V_buttom = trial["V_buttom"]
    H_buttom = trial["H_buttom"]
    response_recorded = False
    response_key = None
    response_time = None
    # Record report onset
    report_start_TR = tr_count
    report_start_time = clock.getTime()
    logFile.write(
        f"Trial {num_trial} REPORT started at "
        f"TR {report_start_TR}, time {report_start_time:.6f} sec\n"
    )
    report_trigger_end_TR = report_end_TR

    while tr_count < report_trigger_end_TR:
        check_for_escape()
        # Check scanner trigger
        check_TR_trigger()
        # Draw report instruction
        buttom_instruct(V_buttom, H_buttom)
        # Check participant response
        keys = event.getKeys(keyList=['1', '2', '3', '4'],timeStamped=clock)
        if keys and not response_recorded:
            response_key, response_time = keys[0]
            response_recorded = True
            trial["ResponseKey"] = response_key
            trial["ResponseTime"] = response_time
            trial["ResponseRT"] = response_time - report_start_time
        # Keep confirmation on screen
        if response_recorded:
            if response_key == V_buttom:
                confirm_report_V.draw()
            elif response_key == H_buttom:
                confirm_report_H.draw()
        myWin.flip()
    # ========================================================
    # REPORT FINISHED
    report_end_time = clock.getTime()
    logFile.write(f"Trial {num_trial} REPORT ended at " f"TR {tr_count}, time {report_end_time:.6f} sec\n")
    logFile.write(f"Time at the end of trial {num_trial} is " f"{report_end_time:.6f} sec\n")
    logFile.write(
        f"Trial {trial['Trial']} INVALID Response: "
        f"{trial['invalid_ResponseKey']} at "
        f"{trial['invalid_ResponseTime']} sec\n"
    )
    logFile.write(
        f"Trial {trial['Trial']} Response: "
        f"{trial['ResponseKey']} at "
        f"{trial['ResponseTime']} sec\n"
    )
# ============================================================
# final 3 TR fixation using TR counting 
final_fix_duration_TR = 3
final_fix_TR = tr_count + final_fix_duration_TR
# First 2 TRs using scanner triggers
while tr_count < final_fix_TR - 1:
    check_for_escape()
    check_TR_trigger()
    dotFix.draw()
    myWin.flip()
# ============================================================
# FINAL TR: time-based
# The trigger that ended the loop above marks the beginning of the final fixation TR.
final_TR_start_time = last_trigger_time
while clock.getTime() - final_TR_start_time < TR:
    check_for_escape()
    dotFix.draw()
    myWin.flip()
#========================================================
# Convert conditions to DataFrame
conditions_df = pd.DataFrame(conditions)

# Adding expected press and success True/False
conditions_df['expected_key'] = conditions_df.apply(
    lambda row: row['H_buttom'] if (row['Instruct_V_H']=='horizontal' and row['illusory_physical']=='illusory') else
                row['V_buttom'] if (row['Instruct_V_H']=='vertical' and row['illusory_physical']=='illusory') else
                row['V_buttom'] if (row['Instruct_V_H']=='horizontal' and row['illusory_physical']=='physical') else
                row['H_buttom'],axis=1)

conditions_df['success'] = conditions_df['ResponseKey']==conditions_df['expected_key']

# Save responses
conditions_df.to_csv(outFileName+'.csv',index=False)
logFile.write(f"Responses saved to {outFileName}.csv")

#========================================================
# Construct protocol file
protocol_df = conditions_df.set_index([col for col in conditions_df.columns if col not in ['PrecueTime','DelayTime','SwitchTime','ReportTime']])
protocol_df = protocol_df.stack().reset_index()
protocol_df.columns = [*protocol_df.columns[:-2],'Condition','Duration']

cols = list(protocol_df.columns)
cols.insert(1,cols.pop(cols.index('Condition')))
cols.insert(2,cols.pop(cols.index('Duration')))
protocol_df = protocol_df[cols]

# Add initial and final 4-TR fixation
initial_fix = {col:'n/a' for col in protocol_df.columns}
initial_fix.update({'Trial':'n/a','Condition':'Fixation','Duration':4})

final_fix = {col:'n/a' for col in protocol_df.columns}
final_fix.update({'Trial':'n/a','Condition':'Fixation','Duration':4})

protocol_df = pd.concat([pd.DataFrame([initial_fix]),protocol_df,pd.DataFrame([final_fix])],ignore_index=True)

# Timing in TRs
protocol_df['Duration'] = pd.to_numeric(protocol_df['Duration'],errors='raise')
protocol_df['Timestamp'] = protocol_df['Duration'].cumsum()
protocol_df['Onset'] = protocol_df['Timestamp']-protocol_df['Duration']
protocol_df.rename(columns={'Condition':'Stim'},inplace=True)

# Save protocol
protocol_df.to_csv(prtFileName+'.csv',index=False)

#========================================================
# Convert protocol into BIDS events.tsv
vol_BIDS_output_file = Path(BIDSoutput_dir)/f"sub-{expInfo['participant']}_task-volitional_run-{int(expInfo['run']):02d}_events.tsv"

# Select BIDS events
vol_events = protocol_df.loc[protocol_df['Stim'].isin(['Fixation','PrecueTime','DelayTime','SwitchTime','ReportTime'])].copy()
if vol_events.empty:
    raise ValueError('No task events found')

# Convert timing from TRs to seconds
vol_events['Onset'] = pd.to_numeric(vol_events['Onset'],errors='raise')
vol_events['Duration'] = pd.to_numeric(vol_events['Duration'],errors='raise')
vol_events['onset'] = vol_events['Onset']*TR
vol_events['duration'] = vol_events['Duration']*TR

# Event labels
stim_mapping = {'Fixation':'fixation','PrecueTime':'precue','DelayTime':'delay','SwitchTime':'switch','ReportTime':'report'}
vol_events['phase'] = vol_events['Stim'].astype(str).str.strip().replace(stim_mapping)
vol_events['instructed_axis'] = vol_events['Instruct_V_H'].astype(str).str.strip().str.lower()
vol_events['trial_type'] = np.where(vol_events['phase']=='fixation','fixation',vol_events['phase']+'_'+vol_events['instructed_axis'])

# Trial metadata
vol_events['trial'] = pd.to_numeric(vol_events['Trial'],errors='coerce').astype('Int64')
vol_events['quartet_order'] = vol_events['QuartetOrder'].astype(str).str.strip()
vol_events['stimulus_type'] = vol_events['illusory_physical'].astype(str).str.strip().str.lower()
vol_events.loc[vol_events['phase']=='fixation','stimulus_type'] = 'n/a'

# Success
vol_events['success'] = vol_events['success'].astype(str).str.strip().str.lower().replace({'true':'1','false':'0','':'n/a'})
vol_events.loc[vol_events['phase']=='fixation','success'] = 'n/a'

# Response information
vol_events['response_key'] = vol_events['ResponseKey'].astype(str).str.strip()
vol_events['expected_key'] = vol_events['expected_key'].astype(str).str.strip()
vol_events['response_time'] = pd.to_numeric(vol_events['ResponseTime'],errors='coerce')
vol_events['response_onset'] = vol_events['response_time']

# Optional metadata
optional_column_mapping = {
    'V_buttom':'vertical_button',
    'H_buttom':'horizontal_button',
    'vertical':'vertical_cue',
    'horizontal':'horizontal_cue',
    'invalid_ResponseKey':'invalid_response_key',
    'invalid_ResponseTime':'invalid_response_time'
}

for original_column,bids_column in optional_column_mapping.items():
    if original_column in vol_events.columns:
        vol_events[bids_column] = vol_events[original_column]

# Final BIDS columns
output_columns = ['onset','duration','trial_type','stimulus_type','success']
vol_events = vol_events[output_columns].sort_values(['onset','trial_type']).reset_index(drop=True)

# Validate
if vol_events['onset'].isna().any():
    raise ValueError('Missing onset values')
if vol_events['duration'].isna().any():
    raise ValueError('Missing duration values')
if (vol_events['onset']<0).any():
    raise ValueError('Negative onset found')

# Save BIDS events.tsv
vol_events.to_csv(vol_BIDS_output_file,sep='\t',index=False,na_rep='n/a',float_format='%.3f')
print(f'Saved {vol_BIDS_output_file}')

#========================================================
# End experiment
endText.draw()
myWin.flip()
core.wait(2)
os.chdir(parentDir)
myWin.close()
core.quit()