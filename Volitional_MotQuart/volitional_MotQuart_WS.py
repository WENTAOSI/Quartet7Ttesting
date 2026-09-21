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
BIDSoutput_dir = os.path.join('BIDS_events', expInfo['participant'], 'func')
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
num_trials = 12
# Initialize parameters these are time as seconds () 
total_time = 13
report = 2
precue = [4, 6]
delay = [4, 6]
switch = [1]
total_TRs = int(total_time * num_trials)
# valid combinations of precue, delay, switch that sum to total_time - report 
# is a list of dictionaries
valid_combinations = [
    {"precue": p, "delay": d, "switch": s}
    for p in precue
    for d in delay
    for s in switch
    if p + d + s == total_time - report
]
# Create a balanced list of QuartetOrder 
quartet_orders = ["quartetPart1, quartetPart2"] * num_trials # to avoid complication we only do one order for now 
# Create balanced list of tone V OR H
instruct_V_H = ["vertical"] * (num_trials // 2) + ["horizontal"] * (num_trials // 2)
#############################################################################################################
# For the current stage all illusory trial NO catch trials are used
'''
# Load the pickle file containing catch trial distribution
with open(os.path.join("Volitional_MotQuart","catch_trials_distribution.pkl"), "rb") as file:
    catch_trials_distribution = pickle.load(file)
current_run_catch_trials = catch_trials_distribution[int(expInfo['run']) - 1]
# Handle catch trials
physical_catch_trials = [trial for trial in current_run_catch_trials if trial in ["V", "H"]]
# Pair one "physical" trial with "vertical" and one with "horizontal"
physical_pairs = [("vertical", "physical") if trial == "V" else ("horizontal", "physical") for trial in physical_catch_trials]
# Remove "vertical" and "horizontal" trials from instruct_V_H to pair with physical trials
for pair in physical_pairs:
    instruct_V_H.remove(pair[0])
# Create the remaining illusory trials
#illusory_physical = ["illusory"] * (num_trials - len(physical_pairs))
#combined_trials = physical_pairs + list(zip(instruct_V_H, illusory_physical))
'''
#####################################################################################################################
combined_trials =  list(zip(instruct_V_H, ['illusory']*num_trials))
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
    # choose a valid combination randomly from the pre-defined valid combinations
    chosen_combo = random.choice(valid_combinations)
    precue_choice = chosen_combo["precue"]
    delay_choice = chosen_combo["delay"]
    switch_choice = chosen_combo["switch"]
    # Store the trial data
    conditions.append({
        "Trial": trial,
        "PrecueTime": precue_choice,
        "DelayTime": delay_choice,
        "SwitchTime": switch_choice,
        "ReportTime": report,
        "QuartetOrder": this_quartet_order,
        "Instruct_V_H": this_instruct_V_H,
        "illusory_physical": this_illusory_physical,
        "V_buttom": this_V,
        "H_buttom": this_H
    })
# Convert to a DataFrame for visualization or saving
conditions_df = pd.DataFrame(conditions)
print(f'Trial {trial}: Precue {precue_choice}, Delay {delay_choice}, Switch {switch_choice}, Report {report}, Total {precue_choice + delay_choice + switch_choice + report}')
# %% STIMULI
# INITIALISE SOME STIMULI
SquareSize = 1.0  # 1.1 #1.8
logFile.write('SquareSize=' + str(SquareSize) + '\n')
dotFix = visual.Circle(myWin, autoLog=False, name='dotFix', units='deg',radius=0.1, pos=apply_global_offset((0,0), global_offset), fillColor='red', lineColor='red' )
Square = visual.GratingStim(myWin, autoLog=False, name='Square', tex=None, units='deg', size=(SquareSize, SquareSize), color= squareColor)
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
    circle = visual.Circle(myWin, autoLog=False, units='deg', radius=circle_size/2, pos=pos, lineColor=squareColor, lineWidth=2, fillColor=None)
    circles.append(circle)
# Generate H/V letters in the circle
Hs = []; Vs = []
for pos in positions:
    H = visual.TextStim(win=myWin, color='white', height=circle_size-0.2,text='H', pos=pos)
    V = visual.TextStim(win=myWin, color='white', height=circle_size-0.2,text='V', pos=pos)
    Hs.append(H); Vs.append(V)
blue_Square = visual.GratingStim(myWin,autoLog=False,name='Square',tex=None,units='deg',size=(SquareSize*1.5, SquareSize*1.5),color='blue',pos=apply_global_offset((0,4), global_offset))
blue_Square.name = "BLUE"
red_Square = visual.GratingStim(myWin,autoLog=False,name='Square',tex=None,units='deg',size=(SquareSize*1.5, SquareSize*1.5), color='red',pos=apply_global_offset((0,4), global_offset))
red_Square.name = "RED"
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
    pos=apply_global_offset(base_pos=(0,-4), global_offset=global_offset)
    )
confirm_report_H = visual.TextStim(
    win=myWin, color='white', height=0.5,
    text='You have pressed HORIZONTAL',
    pos=apply_global_offset(base_pos=(0,-4), global_offset=global_offset)
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
if run_number <= 6:
    color_mapping = {"vertical": blue_Square, "horizontal": red_Square}
    for condition in conditions:
        condition["vertical"] = "blue_Square"
        condition["horizontal"] = "red_Square"
else:
    color_mapping = { "vertical": red_Square, "horizontal": blue_Square}
    for condition in conditions:
        condition["vertical"] = "red_Square"
        condition["horizontal"] = "blue_Square" 
        
mapping_instruct_v = visual.TextStim(
    win=myWin,color='white', height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text=(f'In this run, be prepared to see VERTICAL motion\n' f'after {color_mapping["vertical"].name} onset')
    )
mapping_instruct_h = visual.TextStim(
    win=myWin, color="white", height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text=(f'In this run, be prepared to see HORIZONTAL motion\n' f'after {color_mapping["horizontal"].name} onset')
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
    Hs[int(horizontal_buttom)-1].draw()  
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
    #========================================================
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
        # Instruction cue for first 2 seconds
        show_instruction = (clock.getTime() - delay_start_clock < 2)
        if show_instruction:
            color_mapping[trial["Instruct_V_H"]].draw()
        if trial["QuartetOrder"] == "quartetPart1, quartetPart2":
            quartetPart1(HoriDist, VertiDist)
        elif trial["QuartetOrder"] == "quartetPart2, quartetPart1":
            quartetPart2(HoriDist, VertiDist)
        myWin.flip()
    #========================================================
    #SWITCH
    invalid_key = None
    invalid_timestamp = None
    if trial["illusory_physical"] == "illusory":
        while tr_count < switch_end_TR:
            check_for_escape()
            check_TR_trigger()
            keys = event.getKeys(keyList=['1', '2', '3', '4'],timeStamped=clock)
            if keys and invalid_key is None:
                invalid_key, invalid_timestamp = keys[0]
            if trial["QuartetOrder"] == "quartetPart1, quartetPart2":
                quartetPart2(HoriDist, VertiDist)
            elif trial["QuartetOrder"] == "quartetPart2, quartetPart1":
                quartetPart1(HoriDist, VertiDist)
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
    # ========================================================
    # SET REPORT TR BOUNDARY
    # For all normal trials, use report_end_TR normally.
    # For the LAST trial, stop the trigger-based loop one TR
    # earlier. The final TR will then be displayed using time.
    if num_trial == len(conditions):
        report_trigger_end_TR = report_end_TR - 1
    else:
        report_trigger_end_TR = report_end_TR
    # ========================================================
    # NORMAL TR-BASED REPORT
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
    # FINAL TR OF THE FINAL TRIAL
    # ========================================================
    # The trigger that caused the loop above to finish marks
    # the beginning of the final TR.
    # Do not wait for another trigger. Instead, display the
    # report for one full TR measured from that trigger.
    # ========================================================
    if num_trial == len(conditions):
        final_TR_start_time = last_trigger_time
        logFile.write(
            f"Final clock-timed TR started at "
            f"TR {tr_count}, time {final_TR_start_time:.6f} sec\n"
        )
        while clock.getTime() - final_TR_start_time < TR:
            check_for_escape()
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
    # =======================================================
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
#========================================================    
# Convert conditions to a DataFrame
conditions_df = pd.DataFrame(conditions)
# Adding expected press and success True/False
conditions_df['expected_key'] = conditions_df.apply(
    lambda row: row['H_buttom'] if (row['Instruct_V_H'] == 'horizontal' and row['illusory_physical'] == 'illusory') else
                row['V_buttom'] if (row['Instruct_V_H'] == 'vertical' and row['illusory_physical'] == 'illusory') else
                row['V_buttom'] if (row['Instruct_V_H'] == 'horizontal' and row['illusory_physical'] == 'physical') else
                row['H_buttom'], axis=1
)
conditions_df['success'] = conditions_df.apply(lambda row: True if (row['illusory_physical'] == 'illusory' and row['ResponseKey'] == row['expected_key']) else '', axis=1)
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
protocol_df['Duration'] = protocol_df['Duration'] 
# Create Timestamp column with cumulative time per Trial
protocol_df['Timestamp'] = protocol_df['Duration'].cumsum()
# Create Onset time column
protocol_df['Onset'] = protocol_df['Timestamp'] - protocol_df['Duration']
# change 'Condition' column label to 'Stim'
protocol_df.rename(columns={'Condition': 'Stim'}, inplace=True)
# Save protocol DataFrame to the protocol folder as a CSV file
protocol_df.to_csv(prtFileName + '.csv', index=False)

#========================================================
# converts protocol file into BIDS csv 
vol_BIDS_output_file = (Path(BIDSoutput_dir) / f"{expInfo['participant']}_task-volitional_run-{int(expInfo['run']):02d}_events.tsv")

vol_events = protocol_df.loc[protocol_df["Stim"].isin(["PrecueTime", "DelayTime",])].copy()
if vol_events.empty:
    raise ValueError(f"No PrecueTime or DelayTime rows found in ")

# Convert timing from TRs to seconds
vol_events["Onset"] = pd.to_numeric(vol_events["Onset"], errors="raise",)
vol_events["Duration"] = pd.to_numeric(vol_events["Duration"],errors="raise",)
vol_events["onset"] = vol_events["Onset"] * TR
vol_events["duration"] = vol_events["Duration"] * TR

 # Clean phase and instruction labels
stim_mapping = {"PrecueTime": "precue", "DelayTime": "delay",}
vol_events["phase"] = (vol_events["Stim"].astype(str).str.strip().replace(stim_mapping))
vol_events["instructed_axis"] = (vol_events["Instruct_V_H"].astype(str).str.strip().str.lower())

# trial_type examples: # precue_horizontal # delay_horizontal # precue_vertical   # delay_vertical
vol_events["trial_type"] = (vol_events["phase"]+ "_" + vol_events["instructed_axis"])

# Preserve trial metadata
vol_events["trial"] = pd.to_numeric(vol_events["Trial"], errors="raise",).astype(int)
vol_events["quartet_order"] = (vol_events["QuartetOrder"].astype(str).str.strip())
vol_events["stimulus_type"] = (vol_events["illusory_physical"].astype(str).str.strip().str.lower())
# Convert success to consistent lowercase text BIDS permits extra columns containing strings.
vol_events["success"] = (vol_events["success"].astype(str).str.strip().str.lower().replace({ "true": "1","false": "0"}))
# Preserve response information
vol_events["response_key"] = (vol_events["ResponseKey"].astype(str).str.strip())
vol_events["expected_key"] = (vol_events["expected_key"].astype(str).str.strip())
vol_events["response_time"] = pd.to_numeric(vol_events["ResponseTime"],errors="coerce",)
# The logged ResponseTime appears to be measured from the beginning of the run. Keep it as response_onset rather than calling it reaction time.
vol_events["response_onset"] = vol_events["response_time"]
# Optional additional metadata
optional_column_mapping = {
            "V_buttom": "vertical_button",
            "H_buttom": "horizontal_button",
            "vertical": "vertical_cue",
            "horizontal": "horizontal_cue",
            "invalid_ResponseKey": "invalid_response_key",
            "invalid_ResponseTime": "invalid_response_time",
        }
for original_column, bids_column in optional_column_mapping.items():
    if original_column in vol_events.columns:
        vol_events[bids_column] = vol_events[original_column]
    # Select and order output columns
    output_columns = ["onset","duration","trial_type","success"]
    '''
    # Add optional columns if needed
    optional_output_columns = [
            "vertical_button",
            "horizontal_button",
            "vertical_cue",
            "horizontal_cue",
            "invalid_response_key",
            "invalid_response_time",
            "trial",
            "phase",
            "instructed_axis",
            "quartet_order",
            "stimulus_type",
            "response_key",
            "response_onset",
            "expected_key"
        ]
    output_columns.extend(
            column
            for column in optional_output_columns
            if column in vol_events.columns
        )
    '''
    vol_events = vol_events[output_columns]
    # Sort and validate
    vol_events = (vol_events.sort_values(["onset","trial_type"]).reset_index(drop=True))
    if vol_events["onset"].isna().any():
        raise ValueError(f"Missing onset values in df")
    if vol_events["duration"].isna().any():
        raise ValueError(f"Missing duration values in df")
    if (vol_events["onset"] < 0).any():
            raise ValueError(f"Negative onset found in df")

event_ends = (vol_events["onset"] + vol_events["duration"])
# Save BIDS events.tsv
vol_events.to_csv(vol_BIDS_output_file,sep="\t",index=False,na_rep="n/a",float_format="%.3f")
print(f"Saved {vol_BIDS_output_file}")

# End of experiment 
endText.draw()
myWin.flip()
core.wait(2) # wait for 2 sec

os.chdir(parentDir)
myWin.close()
core.quit()    