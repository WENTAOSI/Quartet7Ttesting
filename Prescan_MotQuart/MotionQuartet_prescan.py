#!/Applications/PsychoPy.app/Contents/Resources/bin/python
##!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
TOP-DOWN MODULATION IN ILLUSORY QUARTET MOTION
Pre-scan Psychophysics session (Genc et al., 2011)

Lastly updated: June 25 2025
"""

# %% IMPORT LIBRARIES
# ==============================================================================
from psychopy import visual, event, core, monitors, logging, gui, data, misc
import numpy as np
import os
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
import markdown


# %% SCREEN AND SYSTEM CONFIG
# ==============================================================================
# Eccentricity properties
circle_dva = 6  # Diameter in dva
circle_radius = circle_dva / 2  # Radius in dva

# specificy background color
backColor = [-0.5, -0.5, -0.5]  # from -1 (black) to 1 (white)

# specificy square color
squareColor = np.multiply(backColor, -1)  # from -1 (black) to 1 (white)    Back dark grey; square light grey


# %% SAVING and LOGGING
# ==============================================================================
# Store info about experiment and experimental run
expName = 'Prescan_MotQuart'  # set experiment name here
expInfo = {
    'participant': 'test',
    'sub_id': '99'
    }

# Monitor Info
my_monitors = {
    "Beaver": {"mon_dist": 80, "size_cm": (52.3, 0), "size_px": (1920, 1080), "refresh_rate": 60, "screen": 1},
    "Scanner": {"mon_dist": 128.7, "size_cm": (42.8, 0), "size_px": (1920, 1080), "refresh_rate": 60, "screen": 1},
    "TseLab": {"mon_dist": 63, "size_cm": (47.2, 0), "size_px": (1920, 1080), "refresh_rate": 60, "screen": 1},
    "TongLab": {"mon_dist": 39, "size_cm": (38, 0), "size_px": (1600, 1200), "refresh_rate": 85, "screen": 0}
}

# Create GUI
dlg = gui.Dlg(title=expName)
dlg.addField("Subject ID:", expInfo['sub_id'])
dlg.addField("Initials:", expInfo['participant'])
dlg.addField("Monitor:", choices=list(my_monitors.keys()), initial="TongLab")
dlg.addField("Debug:", choices=["Yes", "No"], initial="No")
dlg.show()  # Show GUI
if dlg.OK == False: core.quit()  # user pressed cancel

# Update expInfo with selected values
expInfo.update({
    'sub_id': dlg.data[0],
    'participant': dlg.data[1],
    'monitor': dlg.data[2],
    'debug': dlg.data[3],
})
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
DEBUG = True if expInfo['debug'] == "Yes" else False

# get the path that this script is in and change dir to it
_thisDir = os.path.dirname(os.path.abspath(__file__))  # get current path
os.chdir(os.path.dirname(_thisDir))  # change directory

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
logFileName = logFolderName + os.path.sep + '%s_%s_%s' % (
    expInfo['participant'], expInfo['expName'], expInfo['date'])
# Name and create specific folder for output
outFolderName = dataFolderName + os.path.sep + 'Output'
if not os.path.isdir(outFolderName):
    os.makedirs(outFolderName)
outFileName = outFolderName + os.path.sep + '%s_%s_%s' % (
    expInfo['participant'], expInfo['expName'], expInfo['date'])
# Name and create specific folder for protocol files
prtFolderName = dataFolderName + os.path.sep + 'Protocols'
if not os.path.isdir(prtFolderName):
    os.makedirs(prtFolderName)

# save a log file and set level for msg to be received
logFile = logging.LogFile(logFileName+'.log', level=logging.INFO)
logging.console.setLevel(logging.WARNING)  # set console to receive warningVEs

# %% MONITOR AND WINDOW
# ==============================================================================
# Set Monitor Info
monName = expInfo['monitor']
monInfo = my_monitors[monName]

# Initialize Monitor
moni = monitors.Monitor(name=monName,
                           width=monInfo["size_cm"][0],
                           distance=monInfo["mon_dist"])
moni.setSizePix(monInfo["size_px"])

# Log Monitor Info
logFile.write(f"Monitor: {monName}\n")
logFile.write(f"Distance = {monInfo['mon_dist']} cm\n, Width = {monInfo['size_cm'][0]} cm\n")
logFile.write(f"Pixel Width = {monInfo['size_px'][0]}, Pixel Height = {monInfo['size_px'][1]}\n")
logFile.write(f"Refresh Rate = {monInfo['refresh_rate']} Hz\n")

myWin = visual.Window(monitor=moni,
                        size=monInfo["size_px"],
                        screen=monInfo["screen"],
                        winType='pyglet',
                        allowGUI=False,
                        allowStencil=False,
                        fullscr=True,
                        color=backColor,
                        colorSpace='rgb',
                        units='deg',
                        blendMode='avg',
                        waitBlanking=True
                        )


# %% EXPERIMENT DESIGN
# ==============================================================================
# Parameters for experiments

phase1 = {
    "num_runs": 2,
    "num_trials": 10,
    "duration": 250 / 1000, # ms
    "report_time": 0, # s
    "feedback_time": 1,
    "response_delay": 150 / 1000,
    "ITI": 1,
    "num_practice": 1
}
phase2 = {
    "num_runs": 4,
    "num_trials": 80,
    "duration_frame1": 500 / 1000, # ms
    "duration": 500 / 1000, # ms
    "report_time": 5, # s
    "ITI": 1,
    "feedback_time": 1,
    "response_delay": 150 / 1000,
    "cycle": 1,
    "condition_labels": ["PR-3","PR-2","PR-1","PR","PR+1","PR+2","PR+3","PR+4"]
}
# phase3 = {
#     "num_runs": 1,
#     "num_trials": 80,
#     "duration_frame1": 500 / 1000, # ms
#     "duration": 500 / 1000, # ms
#     "cue_time": 1, 
#     "report_time": 5, # s
#     "ITI": 1,
#     "feedback_time": 1,
#     "response_delay": 150 / 1000,
#     "cycle": 1, 
#     "condition_labels": ["vertical","horizontal"]
# }

# for debugging
if DEBUG is True:
    phase1["num_trials"] = 4
    phase1["num_practice"] = 2
    phase2["num_runs"] = 2
    phase2["num_trials"] = 8
    # phase3["num_trials"] = 4


# %% PHASE 1 (Method of Limits)
# ==============================================================================
phase = phase1
conditions = []
for run in range(1, phase["num_runs"] + 1):
    # counterbalancing factor 1: starting dots
    quartet_orders = [
        "left_tilted" if num == 0 else "right_tilted"
        for num in np.random.permutation(phase["num_trials"]) % 2
    ]
    # counterbalancing factor 2: starting ratio
    ratio_directions = [
        "ascending" if num == 0 else "descending"
        for num in np.random.permutation(phase["num_trials"]) % 2
    ]
    for trial in range(1, phase["num_trials"] + 1):
        conditions.append({
            "Trial": trial + phase["num_trials"] * (run - 1),
            "Run": run,
            "Duration": phase["duration"],
            "RespDelay": phase["response_delay"],
            "FeedbackTime": phase["feedback_time"],
            "ITI": phase["ITI"],
            "QuartetOrder": quartet_orders.pop(),
            "RatioDir": ratio_directions.pop()
        })
phase1_conditions = conditions
phase1_df = pd.DataFrame(conditions)

# %% PHASE 2 (Method of Constant Stimuli)
# ==============================================================================
phase = phase2
conditions = []
for run in range(1, phase["num_runs"] + 1):
    # condition: aspect ratio (8 personalized ratio)
    condition_ratio = [
        phase["condition_labels"][num]
        for num in np.random.permutation(phase["num_trials"]) % len(phase["condition_labels"])
    ]
    # counterbalancing factor: starting dots
    quartet_orders = [
        "left_tilted" if num == 0 else "right_tilted"
        for num in np.random.permutation(phase["num_trials"]) % 2
    ]
    for trial in range(1, phase["num_trials"] + 1):
        conditions.append({
            "Run": run,
            "Trial": trial + phase["num_trials"] * (run - 1),
            "Duration_F1": phase["duration_frame1"],
            "Duration": phase["duration"],
            "ReportTime": phase["report_time"],
            "RespDelay": phase["response_delay"],
            "FeedbackTime": phase["feedback_time"],
            "ITI": phase["ITI"],
            "Cycle": phase["cycle"],
            "ConditionRatio": condition_ratio.pop(),
            "QuartetOrder": quartet_orders.pop()
        })
phase2_conditions = conditions
phase2_df = pd.DataFrame(conditions)

# %% PHASE 3 (Assessment for Volitional Control)
# # ==============================================================================
# phase = phase3
# conditions = []
# for run in range(1, phase["num_runs"] + 1):
#     # condition: color cue (vertical or horizontal)
#     condition_cue = [
#         phase["condition_labels"][num]
#         for num in np.random.permutation(phase["num_trials"]) % len(phase["condition_labels"])
#     ]
#     # counterbalancing factor: starting dots
#     quartet_orders = [
#         "left_tilted" if num == 0 else "right_tilted"
#         for num in np.random.permutation(phase["num_trials"]) % 2
#     ]
#     for trial in range(1, phase["num_trials"] + 1):
#         conditions.append({
#             "Run": run,
#             "Trial": trial + phase["num_trials"] * (run - 1),
#             "Duration_F1": phase["duration_frame1"],
#             "Duration": phase["duration"],
#             "CueTime": phase["cue_time"],
#             "ReportTime": phase["report_time"],
#             "RespDelay": phase["response_delay"],
#             "FeedbackTime": phase["feedback_time"],
#             "ITI": phase["ITI"],
#             "Cycle": phase["cycle"],
#             "ConditionCue": condition_cue.pop(),
#             "QuartetOrder": quartet_orders.pop()
#         })
# phase3_conditions = conditions
# phase3_df = pd.DataFrame(conditions)

# %%  STIMULUS PARAMETERS
# ==============================================================================

# RATIO PRESET
max_rad = np.arctan(3)        # naximum ratio in radians
min_rad = np.arctan(1/3)      # minimum ratio in radians
range_rad = np.linspace(min_rad, max_rad, 154)  # 154 evenly-spaced angles in rad
range_ratio = np.tan(range_rad)  # rad to ratio
list_ratio = {
    "ascending": range_ratio,
    "descending": range_ratio[::-1]
}

# SQUARE (quartets)
SquareSize = 1.0  # in dva
logFile.write('SquareSize=' + str(SquareSize) + '\n')
Square = visual.GratingStim(myWin,
                            autoLog=False,
                            name='Square',
                            tex=None,
                            units='deg',
                            size=(SquareSize, SquareSize),
                            color= squareColor,
                            )

# DOT (fixation)
dotFix = visual.Circle(myWin,
                       autoLog=False,
                       name='dotFix',
                       units='deg',
                       radius=0.1,
                       fillColor='red',
                       lineColor='red'
                       )

# DOT (green)
dotFix_green = visual.Circle(myWin,
                       autoLog=False,
                       name='dotFix',
                       units='deg',
                       radius=0.1,
                       fillColor='green',
                       lineColor='green'
                       )

# %%  INSTRUCTION
# ==============================================================================

def show_text(win, text, color='white', height=0.5, pos=(0, 0)):
    return visual.TextStim(
        win=win,
        color=color,
        height=height,
        text=text,
        pos=pos
    )

triggerText = show_text(myWin, "Experiment will start soon. Press p to continue.")
phase1Text = show_text(myWin, "Please keep your eyes on the red dot at all times.\n\nPress the space bar as soon as you notice\n\na change in the perceived direction.")
pracText = show_text(myWin, "Press any key to begin the practice trials.")
endofpracText = show_text(myWin, "Great job!\nWhen you're ready, press the space bar to begin the experiment.")
bridgeText = show_text(myWin, "You've finished phase 1! Now, move on to phase 2.\n\nPress any key to proceed.")
phase2Text = show_text(myWin, "Please keep your eyes on the red dot at all times.\n\nOnce the presentation ends,\n\npress V or H key to report the direction you perceived.")
phase3Text = show_text(myWin, "Please keep your eyes on the red dot at all times.\n\nTry to perceive the direction indicated in the instructions.\n\nReport the direction you actually perceived.")
reportText = show_text(myWin, "?")
norespText = show_text(myWin, "No response detected.\n\nPlease stay focused for the next trial.")
analysisText = show_text(myWin, "Now analyzing the data ...")
endText = show_text(myWin, "You're all set! Thank you for your participation.")

# Mapping response key to label
key_to_label = {'v': 'vertical', 'h': 'horizontal'}

def show_break(runs_left):
    text = f"{runs_left} run(s) to go!\n\nPlease take a short break and press any key to continue."
    return show_text(myWin, text)

def show_cue(trial_cue):
    text = f"Try to see {trial_cue.upper()}" 
    return show_text(myWin, text, pos = (0,0.5))

# %% TIME PARAMETERS
# ==============================================================================

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
# ==============================================================================

# Calculate the Hori/Verti distances
def ratio2dist(ratio, radius):
    """
    Calculate the horizontal (Hori) and vertical (Verti) distances
    for a point on a circle given the radius and a ratio.
    """
    angle = np.arctan(ratio)  # Ratio to Angle in radians
    Hori = radius * np.cos(angle)  # Horizontal distance
    Verti = radius * np.sin(angle)  # Vertical distance
    return Hori, Verti

# Generate quartets on the screen
def show_quartets(Hori, Verti, pair, is_green=False):
    if pair == 0:  # Left-tilted pair
        Square.setPos((-Hori, Verti))
        Square.draw()
        Square.setPos((Hori, -Verti))
        Square.draw()
    elif pair == 1:  # Right-tilted pair
        Square.setPos((Hori, Verti))
        Square.draw()
        Square.setPos((-Hori, -Verti))
        Square.draw()
    else:
        raise ValueError("Invalid value for pair. Use 0 for left-tilted or 1 for right-tilted.")
    dotFix.draw()
    if is_green:
        dotFix_green.draw()

# Show Fixation
def show_fixation(dotFix, win, duration, is_green=False):
    dotFix.draw()
    if is_green:
        dotFix_green.draw()
    win.flip()
    core.wait(duration)

# # Quartet: left-tilted dot pair
# def quartetPart1(Hori, Verti):
#     Square.setPos((-Hori, Verti))
#     Square.draw()
#     Square.setPos((Hori, -Verti))
#     Square.draw()
#     dotFix.draw()

# # Quartet: right-tilted dot pair
# def quartetPart2(Hori, Verti):
#     Square.setPos((Hori, Verti))
#     Square.draw()
#     Square.setPos((-Hori, -Verti))
#     Square.draw()
#     dotFix.draw()

def check_for_escape():
    """Check if the Escape key is pressed and exit the program."""
    keys = event.getKeys(keyList=['escape'])
    if 'escape' in keys:
        core.quit()

        
# %% INSTRUCTIONS
# ==============================================================================

event.Mouse(visible=False)
logFile.write(f"Start of Experiment {expInfo['expName']}\n")

triggerText.draw()
myWin.flip()
event.waitKeys(keyList=['p'], timeStamped=False)

# Phase 1 instruction
phase1Text.draw()
myWin.flip()
event.waitKeys(timeStamped=False)

# Practice 1 start instruction
pracText.draw()
myWin.flip()
event.waitKeys(timeStamped=False)

# reset clocks
clock.reset()
practice_start_time = clock.getTime()

# %% PRACTICE TRIALS
# ==============================================================================

logFile.write(f'Time at start of practice trials is {clock.getTime()}\n')

for trial in range(phase1["num_practice"]):

    trial_start_time = clock.getTime()

    trial_pair = np.random.randint(0, 2)
    trial_list_ratio = list_ratio[np.random.choice(["ascending", "descending"])]

    response_recorded = False
    response_key = None
    response_time = None

    n_flip = 0
    n_step = 0
    HoriDist, VertiDist = ratio2dist(trial_list_ratio[n_step], circle_radius)

    # Stimuli presentation until response
    while response_recorded is False:
        check_for_escape()
        # switch pair (per 250 ms)
        if clock.getTime() - trial_start_time > phase1["duration"] * (n_flip+1):
            trial_pair = 1 - trial_pair
            n_flip += 1
            # switch ratio (per 1 cycle of quartet)
            if n_flip % 2 == 0: 
                if n_step >= len(trial_list_ratio): # exit if the list is over
                    break
                trial_ratio = trial_list_ratio[n_step]
                HoriDist, VertiDist = ratio2dist(trial_ratio, circle_radius)
                n_step += 1
        show_quartets(HoriDist, VertiDist, trial_pair)
        myWin.flip()

        # Check for button presses during each frame
        keys = event.getKeys(keyList=['space'], timeStamped=clock)
        if keys and not response_recorded and clock.getTime() - trial_start_time > phase1["response_delay"]:
            response_key, response_time = keys[0]  # Extract the key and timestamp
            response_recorded = True  # Mark the response as recorded

    # Response Feedback
    if response_recorded:
        dotFix_green.draw() #show_quartets(HoriDist, VertiDist, trial_pair, is_green=True)
    else:
        norespText.draw()
    myWin.flip()
    core.wait(phase1["feedback_time"])

    # Intertrial interval
    dotFix.draw()
    myWin.flip()
    core.wait(phase1["ITI"])
            
logFile.write(f'Time at the end of practice trials is {clock.getTime()}\n')    

# End of practice instructions
endofpracText.draw()
myWin.flip()
event.waitKeys(timeStamped=False)

# %% RUN PHASE 1
# ==============================================================================

logFile.write(f'Phase 1: Method of Limits\n')

for trial in phase1_conditions:

    logFile.write(f'Phase 1: Time at start of trial {trial["Trial"]} is {clock.getTime()}\n')
    trial_start_time = clock.getTime()

    trial_pair = 0 if trial["QuartetOrder"] == "left_tilted" else 1
    trial_list_ratio = list_ratio[trial["RatioDir"]]

    response_recorded = False
    response_key = None
    response_time = None

    n_flip = 0
    n_step = 0
    HoriDist, VertiDist = ratio2dist(trial_list_ratio[n_step], circle_radius)

    # Stimuli presentation until response
    while response_recorded is False:
        check_for_escape()
        # switch pair (per 250 ms)
        if clock.getTime() - trial_start_time > trial["Duration"] * (n_flip+1):
            trial_pair = 1 - trial_pair
            n_flip += 1
            # switch ratio (per 1 cycle of quartet)
            if n_flip % 2 == 0: 
                if n_step >= len(trial_list_ratio): # exit if the list is over
                    break
                trial_ratio = trial_list_ratio[n_step]
                HoriDist, VertiDist = ratio2dist(trial_ratio, circle_radius)
                n_step += 1
        show_quartets(HoriDist, VertiDist, trial_pair)
        myWin.flip()

        # Check for button presses during each frame
        keys = event.getKeys(keyList=['space'], timeStamped=clock)
        if keys and not response_recorded and clock.getTime() - trial_start_time > trial["RespDelay"]:
            response_key, response_time = keys[0]  # Extract the key and timestamp
            response_recorded = True  # Mark the response as recorded    
            # Record response
            trial["ResponseKey"] = response_key
            trial["ResponseTime"] = response_time - trial_start_time

    # Record stimulus on response
    trial["ResponseFlip"] = n_flip 
    trial["ResponseRatio"] = trial_ratio

    logFile.write(f'Phase 1: Time at the end of trial {trial["Trial"]} is {clock.getTime()}\n')
    
    # Response Feedback
    if response_recorded is False:
        norespText.draw()
        trial["ResponseKey"] = "None"
        trial["ResponseTime"] = "None"

    show_fixation(dotFix, myWin, trial["FeedbackTime"], is_green=True)
    show_fixation(dotFix, myWin, trial["ITI"])
    
    # Log response
    logFile.write(f"Phase 1: Trial {trial['Trial']} Response: {trial['ResponseKey']} at {trial['ResponseTime']} sec\n")

    # Interblock break
    if trial["Trial"] % phase1["num_trials"] == 0 and trial["Trial"] < phase1["num_trials"] * phase1["num_runs"]:
        show_break(phase1["num_runs"] - trial["Run"]).draw()
        myWin.flip()
        event.waitKeys()
    
# Convert conditions to a DataFrame
phase1_df = pd.DataFrame(phase1_conditions)
# Save responses DataFrame to the Output folder as a CSV file
phase1_df.to_csv(outFileName + '_p1.csv', index=False)

# Log the saving process 
logFile.write(f"Phase 1 Responses saved to {outFileName}.csv\n")

# End of phase 1 instrunctions
bridgeText.draw()
myWin.flip()
event.waitKeys()

# %% CALCULATE PERSONALIZED ASPECT RATIO
# ==============================================================================

# Filter by 'ascending' and 'descending' bins
ascending_bin = phase1_df[phase1_df["RatioDir"] == "ascending"]
descending_bin = phase1_df[phase1_df["RatioDir"] == "descending"]

# Calculate means in angles
ascending_mean_rad = np.arctan(
    ascending_bin[ascending_bin["ResponseKey"].notna()]["ResponseRatio"].astype(float)
).mean()
descending_mean_rad = np.arctan(
    descending_bin[descending_bin["ResponseKey"].notna()]["ResponseRatio"].astype(float)
).mean()
overall_mean_rad = (ascending_mean_rad + descending_mean_rad) / 2

# Calculate means in ratio
ascending_mean = np.tan(ascending_mean_rad)
descending_mean = np.tan(descending_mean_rad)
overall_mean_ratio = np.tan(overall_mean_rad)

# Log the estimated parity ratio
logFile.write(f"Ascending Mean Ratio: {ascending_mean:.4f}\n")
logFile.write(f"Descending Mean Ratio: {descending_mean:.4f}\n")
logFile.write(f"Overall Mean Ratio: {overall_mean_ratio:.4f}\n")

# Generate personalized aspect ratio
mystep_rad = 0.075
myrange_rad = [overall_mean_rad + i * mystep_rad for i in range(-3, 5)]
myrange_ratio = np.tan(myrange_rad) 
subject_ratio = {label: value for label, value in zip(phase2["condition_labels"], myrange_ratio)}

# Log the personalized aspect ratio
logFile.write("Personalized Aspect Ratios:\n")
for label, ratio in subject_ratio.items():
    logFile.write(f"{label}: {ratio:.4f}\n")


# %% RUN PHASE 2
# ==============================================================================

# Phase 2 instruction
phase2Text.draw()
myWin.flip()
event.waitKeys()
logFile.write(f'Phase 2: Method of Constant Stimuli\n')

for trial in phase2_conditions:
     
    logFile.write(f'Phase 2: Time at start of trial {trial["Trial"]} is {clock.getTime()}\n')

    trial_pair = 0 if trial["QuartetOrder"] == "left_tilted" else 1
    trial_ratio = subject_ratio[trial["ConditionRatio"]]
    trial["trial_ratio"] = trial_ratio # record in the df what the ratio is 
    HoriDist, VertiDist = ratio2dist(trial_ratio, circle_radius)

    response_recorded = False
    response_key = None
    response_time = None

    trial_start_time = clock.getTime()
    n_flip = 0
    next_flip_time = trial["Duration_F1"] # duration for 1st frame

    # Stimuli presentation (for predetermined cycles)
    while 1:
        check_for_escape()
        # switch pair
        if clock.getTime() - trial_start_time > next_flip_time:
            trial_pair = 1 - trial_pair
            n_flip += 1
            next_flip_time += trial["Duration"]
        if n_flip >= trial["Cycle"] * 2:
            break
        show_quartets(HoriDist, VertiDist, trial_pair)
        myWin.flip()
    trial_resp_window = clock.getTime()

    # Check for button presses (until response or max report time)
    while response_recorded is False:
        check_for_escape()
        keys = event.getKeys(keyList=['v', 'h'], timeStamped=clock)
        if keys and not response_recorded and clock.getTime() - trial_resp_window > trial["RespDelay"]:  # Process only the first response
            response_key, response_time = keys[0]  # Extract the key and timestamp
            response_recorded = True  # Mark the response as recorded    
            response_label = key_to_label[response_key]
            # Record response
            trial["ResponseKey"] = response_key
            trial["ResponseTime"] = response_time - trial_resp_window
            trial["ResponseLabel"] = response_label
        reportText.draw()
        myWin.flip()

    logFile.write(f'Phase 2: Time at the end of trial {trial["Trial"]} is {clock.getTime()}\n')
    
    # Response Feedback
    if response_recorded is False:
        norespText.draw()
        trial["ResponseKey"] = "None"
        trial["ResponseTime"] = "None"
    show_fixation(dotFix, myWin, trial["FeedbackTime"], is_green=True)
    show_fixation(dotFix, myWin, trial["ITI"])
        
    # Log response
    logFile.write(f"Phase 2: Trial {trial['Trial']} Response: {trial['ResponseKey']} at {trial['ResponseTime']} sec\n")

    # Interblock break
    if trial["Trial"] % phase2["num_trials"] == 0 and trial["Trial"] < phase2["num_trials"] * phase2["num_runs"]:
        show_break(phase2["num_runs"] - trial["Run"]).draw()
        myWin.flip()
        event.waitKeys()
    
# Convert conditions to a DataFrame
phase2_df = pd.DataFrame(phase2_conditions)
# Save responses DataFrame to the Output folder as a CSV file
phase2_df.to_csv(outFileName + '_p2.csv', index=False)

# Log the saving process 
logFile.write(f"Phase 2 Responses saved to {outFileName}.csv\n")


# %% Wrap up the data ...
# ==============================================================================

analysisText.draw()
myWin.flip()

baseFileName = outFolderName + os.path.sep + '%s_%s' % (
    expInfo['participant'], expInfo['expName'])

# Phase 1 data: Density Plot
# ==============================================================================

ascending_bin["RatioDir"] = "Ascending"
descending_bin["RatioDir"] = "Descending"
combined_data = pd.concat([ascending_bin, descending_bin], ignore_index=True)

p1_data = combined_data[combined_data["ResponseKey"].notna()].copy()
p1_data["ResponseRad"] = np.arctan(p1_data["ResponseRatio"])

summary_table = (
    p1_data.groupby("RatioDir")
    .agg(
        MeanRad=("ResponseRad", "mean"),
        SDRad=("ResponseRad", "std"),
        MinRad=("ResponseRad", "min"),
        MaxRad=("ResponseRad", "max"),
        n=("ResponseRad", "count")
    )
    .reset_index()
)

# Add ratio versions
summary_table["MeanRatio"] = np.tan(summary_table["MeanRad"])
summary_table["SDRatio"] = np.tan(summary_table["SDRad"])
summary_table["MinRatio"] = np.tan(summary_table["MinRad"])
summary_table["MaxRatio"] = np.tan(summary_table["MaxRad"])


# %% Density Plot in Radian

plt.figure(figsize=(6, 2))
sns.kdeplot(
    data=p1_data,
    x="ResponseRad",
    hue="RatioDir",
    fill=True,
    common_norm=False,
    alpha=0.5,
    palette={"Ascending": "blue", "Descending": "red"}
)

for _, row in summary_table.iterrows():
    plt.axvline(row["MeanRad"], color="blue" if row["RatioDir"] == "Ascending" else "red",
                linestyle="--", linewidth=2)

plt.title("Density Plot of Response by Direction")
plt.xlabel("Response (Radian)")
plt.ylabel("Density")
plt.tight_layout()
plt.savefig(baseFileName + '_Phase1_RadDensity.png', dpi=300)
plt.close()

# %% Density Plot in Ratio
plt.figure(figsize=(6, 2))
sns.kdeplot(
    data=p1_data,
    x="ResponseRatio",
    hue="RatioDir",
    fill=True,
    common_norm=False,
    alpha=0.5,
    palette={"Ascending": "blue", "Descending": "red"}
)

for _, row in summary_table.iterrows():
    plt.axvline(row["MeanRatio"], color="blue" if row["RatioDir"] == "Ascending" else "red",
                linestyle="--", linewidth=2)

plt.title("Density Plot of Response by Direction")
plt.xlabel("Response (Ratio)")
plt.ylabel("Density")
plt.tight_layout()
plt.savefig(baseFileName + '_Phase1_RatDensity.png', dpi=300)
plt.close()

# %% Phase 2 data: Estimate Psychometric Curve
# ==============================================================================

# Filter out invalid responses
p2_data = phase2_df[phase2_df["ResponseKey"] != "None"].copy()

p2_data["ConditionRat"] = p2_data["ConditionRatio"].map(subject_ratio)
p2_data["ConditionRad"] = np.arctan(p2_data["ConditionRat"].astype(float))
p2_data["ResponseBinary"] = p2_data["ResponseLabel"].apply(lambda x: 1 if x == "vertical" else 0)

# Model fitting 
X = p2_data[["ConditionRad"]].values  # (n, 1) shape
y = p2_data["ResponseBinary"].values  # (n,) shape

model = LogisticRegression(solver="lbfgs", fit_intercept=True, max_iter=1000)
model.fit(X, y)

# Probability prediction for each ratio
p2_data["PredictedProb"] = model.predict_proba(X)[:, 1]  # P(y=1)

# Continuous predictions for visualization
x_pred = np.linspace(X.min(), X.max(), 100).reshape(-1, 1)
y_pred = model.predict_proba(x_pred)[:, 1]*100

# calculate the PSE
intercept = model.intercept_[0]
slope = model.coef_[0][0]
x_intercept = -intercept / slope

# Summarize data for visualization
data_summary = (
    p2_data.groupby("ConditionRad")
    .agg(
        vertical_percent=("ResponseBinary", lambda x: x.mean() * 100),
        n=("ResponseBinary", "count"),
        se=("ResponseBinary", lambda x: x.std(ddof=1) / np.sqrt(len(x)) * 100)
    )
    .reset_index()
)
data_summary["ConditionRatio"] = np.tan(data_summary["ConditionRad"])

# %% Plotting the results
plt.figure(figsize=(6, 6))

# fitted curve
plt.plot(x_pred, y_pred, color='red', label='Fitted curve')

# probability prediction + error bars
plt.errorbar(
    data_summary["ConditionRad"], 
    data_summary["vertical_percent"],
    yerr=data_summary["se"],
    fmt='o', 
    color='blue', 
    ecolor='blue', 
    capsize=3, 
    label='Observed mean ± SE'
)

# PSE line
plt.axvline(x=x_intercept, color='black', linestyle='--')
plt.axhline(y=50, color='black', linestyle='--')
plt.text(x_intercept+0.05, 3, f"PSE: {x_intercept:.2f}", fontweight='bold', ha='center', va='bottom')

# Labels
plt.xlabel("Angle in Radian")
plt.ylabel("Prediction for Vertical Response (%)")
plt.title("Psychometric Curve with Quartet Angle")
plt.legend()
plt.tight_layout()
plt.savefig(baseFileName + '_Phase2_RadCurve.png', dpi=300)
plt.close()

# %% Plotting in ratio
x_pred_rad = np.linspace(p2_data["ConditionRad"].min(), p2_data["ConditionRad"].max(), 100).reshape(-1, 1)
y_pred_prob = model.predict_proba(x_pred_rad)[:, 1] * 100
x_pred_ratio = np.tan(x_pred_rad.flatten())

intercept = model.intercept_[0]
slope = model.coef_[0][0]
x_intercept_rad = -intercept / slope
x_intercept_ratio = np.tan(x_intercept_rad)

plt.figure(figsize=(6, 6))

# Fitted curve
plt.plot(x_pred_ratio, y_pred_prob, color='red', label='Fitted curve')

# predictions + error bars
plt.errorbar(
    data_summary["ConditionRatio"], 
    data_summary["vertical_percent"],
    yerr=data_summary["se"],
    fmt='o', 
    color='blue', 
    ecolor='blue', 
    capsize=3, 
    label='Observed mean ± SE'
)

# PSE
plt.axvline(x=x_intercept_ratio, color='black', linestyle='--')
plt.axhline(y=50, color='black', linestyle='--')
plt.text(x_intercept_ratio+0.12, 3, f"PSE: {x_intercept_ratio:.2f}", fontweight='bold', ha='center', va='top')

# Labels
plt.xlabel("Aspect Ratio")
plt.ylabel("Prediction for Vertical Response (%)")
plt.title("Psychometric Curve with Quartet Ratio")
plt.tight_layout()
plt.legend()
plt.savefig(baseFileName + '_Phase2_RatCurve.png', dpi=300)
plt.close()


# %% Save Markdown
# ==============================================================================

pse_rad = x_intercept
pse_ratio = x_intercept_ratio
subject_ratio_df = pd.DataFrame(list(subject_ratio.items()), columns=["Condition Label", "Aspect Ratio"])
subject_ratio_df["Aspect Ratio"] = subject_ratio_df["Aspect Ratio"].map(lambda x: f"{x:.4f}")

# Markdown
mdFileName = baseFileName + '_summary.md'

def img_md(filename, width=400):
    return f'<img src="{filename}" width="{width}">'

md_lines = [
    f"# Summary Report for Participant {expInfo['participant']}\n",
    f"**Experiment:** {expInfo['expName']}\n",
    f"**Date:** {expInfo['date']}\n",

    "## Phase 1 - Method of Limits\n",
    "**Radian Density Plot:**\n",
    img_md(os.path.basename(baseFileName) + "_Phase1_RadDensity.png") + "\n\n",

    "**Ratio Density Plot:**\n",
    img_md(os.path.basename(baseFileName) + "_Phase1_RatDensity.png") + "\n\n",

    "**Summary Table:**\n",
    summary_table.to_markdown(index=False), 
    "\n\n",

    "## Personalized Aspect Ratios\n",
    subject_ratio_df.to_markdown(index=False),
    "\n\n",

    "## Phase 2 - Psychometric Curve\n",
    "**Radian Version:**\n",
    img_md(os.path.basename(baseFileName) + "_Phase2_RadCurve.png") + "\n\n",

    "**Ratio Version:**\n",
    img_md(os.path.basename(baseFileName) + "_Phase2_RatCurve.png") + "\n\n",

    f"**Estimated PSE (in radian):** `{pse_rad:.4f}`\n",
    f"**Estimated PSE (in ratio):** `{pse_ratio:.4f}`\n",
]
# Save as MD
with open(mdFileName, 'w', encoding='utf-8') as f:
    for line in md_lines:
        f.write(line + '\n')

with open(mdFileName, 'r', encoding='utf-8') as f:
    md_text = f.read()

html = markdown.markdown(md_text, extensions=['tables'])
htmlFileName = mdFileName.replace('.md', '.html')

with open(htmlFileName, 'w', encoding='utf-8') as f:
    f.write(html)

# %% RUN PHASE 3
# ==============================================================================

# # Phase 3 instruction
# phase3Text.draw()
# myWin.flip()
# event.waitKeys()
# logFile.write(f'Phase 3: Testing Volitional Control \n')

# # preset for personalized PR
# trial_ratio = estimated_PR
# HoriDist, VertiDist = ratio2dist(trial_ratio, circle_radius)

# for trial in phase3_conditions:
     
#     logFile.write(f'Phase 3: Time at start of trial {trial["Trial"]} is {clock.getTime()}\n')

#     trial_pair = 0 if trial["QuartetOrder"] == "left_tilted" else 1
#     trial_cue = trial["ConditionCue"]

#     response_recorded = False
#     response_key = None
#     response_time = None

#     trial_start_time = clock.getTime()
#     n_flip = 0
#     next_flip_time = trial["Duration_F1"] # duration for 1st frame

#     # Cue presentation
#     show_cue(trial_cue).draw()
#     dotFix.draw()
#     myWin.flip()
#     core.wait(trial["CueTime"])
#     trial_cue_time = clock.getTime()
    
#     # Stimuli presentation (for predetermined cycles)
#     while 1:
#         check_for_escape()
#         # switch pair
#         if clock.getTime() - trial_cue_time > next_flip_time:
#             trial_pair = 1 - trial_pair
#             n_flip += 1
#             next_flip_time += trial["Duration"]
#         if n_flip >= trial["Cycle"] * 2:
#             break
#         show_quartets(HoriDist, VertiDist, trial_pair)
#         myWin.flip()
#     trial_resp_window = clock.getTime()

#     # Check for button presses (until response or max report time)
#     while response_recorded is False:
#         check_for_escape()
#         keys = event.getKeys(keyList=['v', 'h'], timeStamped=clock)
#         if keys and not response_recorded and clock.getTime() - trial_resp_window > trial["RespDelay"]:  # Process only the first response
#             response_key, response_time = keys[0]  # Extract the key and timestamp
#             response_recorded = True  # Mark the response as recorded    
#             # Control success?
#             response_label = key_to_label[response_key]
#             trial_success = response_label == trial_cue
#             # Record response
#             trial["ResponseKey"] = response_key
#             trial["ResponseTime"] = response_time - trial_resp_window
#             trial["ResponseLabel"] = response_label
#             trial["TrialSuccess"] = trial_success
#         reportText.draw()
#         myWin.flip()

#     logFile.write(f'Phase 3: Time at the end of trial {trial["Trial"]} is {clock.getTime()}\n')
    
#     # Response Feedback
#     if response_recorded is False:
#         norespText.draw()
#         trial["ResponseKey"] = "None"
#         trial["ResponseTime"] = "None"
#         trial["ResponseLabel"] = "None"
#         trial["TrialSuccess"] = "None"
        
#     show_fixation(dotFix, myWin, trial["FeedbackTime"], is_green=True)
#     show_fixation(dotFix, myWin, trial["ITI"])
        
#     # Log response
#     logFile.write(f"Phase 3: Trial {trial['Trial']} Response: {trial['ResponseKey']} at {trial['ResponseTime']} sec\n")
    
# # Convert conditions to a DataFrame
# phase3_df = pd.DataFrame(phase3_conditions)
# # Save responses DataFrame to the Output folder as a CSV file
# phase3_df.to_csv(outFileName + '_p3.csv', index=False)

# # Log the saving process 
# logFile.write(f"Phase 3 Responses saved to {outFileName}.csv\n")

# # %% CALCULATE SUCCESS RATE

# # Filter by 'horizontal' and 'vertical' bins
# horizontal_bin = phase3_df[phase3_df["ConditionCue"] == "horizontal"]
# vertical_bin = phase3_df[phase3_df["ConditionCue"] == "vertical"]

# # Calculate success rate
# horizontal_success_count = horizontal_bin["TrialSuccess"].apply(lambda x: x == True).sum()
# vertical_success_count = vertical_bin["TrialSuccess"].apply(lambda x: x == True).sum()

# horizontal_success_rate = horizontal_success_count / len(horizontal_bin)
# vertical_success_rate = vertical_success_count / len(vertical_bin)

# # Log success rate
# logFile.write(f"Phase 3: Success rate for horizontal cue: {horizontal_success_rate * 100:.2f}%\n")
# logFile.write(f"Phase 3: Success rate for vertical cue: {vertical_success_rate * 100:.2f}%\n")


# %% End of experiment
logFile.write("End of Experiment")
endText.draw()
myWin.flip()
core.wait(5) # wait for 5 s
myWin.close()
event.Mouse(visible=False)

try:
    core.quit()
except SystemExit:
    os._exit(0)

# %%
