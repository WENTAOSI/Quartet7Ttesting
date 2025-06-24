"""
Created on Sat May 07 2022
It presents ambiguous motion quartet stimulus.
MotQuart lastes for 80s, FlickQuart lastes for 16s.
6 repetitions.
The stimulus begins and ends with a fixation condition of 20s.

Total triggers: 616
Total time: 616 = [10min20s]

Psychopy3 (v2020.2.4)
Based on https://github.com/MSchnei/motion_quartet_scripts (@author: Marian.Schneider)

adopted from @author: Alessandra Pizzuti adopted from @author Marian.Schneider

@author: siwentao to integrate EYELINK
"""

from psychopy import visual, event, core, monitors, logging, gui, data, misc, sound
import numpy as np
import pandas as pd
import os
import sys
import time
import pylink
from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy

#%% SET PARAMS
###############################################################################
TRIGGERKEY = 'quoteleft'
# BLOCK DURATIONS [in TR]
# set durations of conditions and baseline
TR = 4.217     # sec in int if whole number  or float 
# specify vertical or horizontal switch buttom 
vertical_buttom = "1"
horizontal_buttom = "2"
ITI_buttom = "4"

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
###############################################################################
# BLOCKS


# if integer TR we can set percise timing 
if TR == 2:
    DurElem = np.array([int(12/TR), int(16/TR), int(96/TR)])  # fix = 12s; flickerQuartet = 16s, AmbiguousQuartet = 96s
    NumQuartets = 3  # set number of repetitions of quartet blocks
    # NOTE: Fixation at the beginning and at the end lasts both for 10 triggers.

elif TR == 4.217:
    DurElem = np.array([4, 4, 32]) # fix = 4 TR; flickerQuartet = 4 TR, AmbiguousQuartet = 36 TR
    NumQuartets = 2  # set number of repetitions of quartet blocks
# if float TR, we have to enforce number of TR rather than time 
'''
else:
    DurElem = np.array([10, 8, 40]) # fix 10 TR; flickerQuartet = 8 TR; AmbiguousQuartet = 40 TR
'''


# fixation = 0; flicker = 1; quartet = 2

Conditions = np.zeros(int(NumQuartets*2))
Conditions[::2] = np.tile([2], NumQuartets)  # every 2nd element
Conditions[1::2] = np.random.permutation(np.tile(np.array([1]), NumQuartets))
Conditions = np.hstack(([0], Conditions, [0]))
Durations = np.zeros(int(len(Conditions)))
for ind in range(0, len(DurElem)):
    Durations[Conditions == ind] = DurElem[ind]

print('Conditions:', Conditions)
print('Durations:', Durations)
# Store info about experiment and experimental run
expName = 'Amb_MotQuart'  # set experiment name here
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

message = visual.TextStim(
    myWin,
    text='Condition',
    pos=(-16, -8)
    )

dotFix = visual.Circle(
    myWin,
    autoLog=False,
    name='dotFix',
    units='deg',
    radius=.1,
    pos=apply_global_offset((0,0), global_offset),
    fillColor='red',
    lineColor='red'
    )

Square = visual.GratingStim(
    myWin,
    autoLog=False,
    name='Square',
    tex=None,
    units='deg',
    size=(SquareSize, SquareSize),
    color=squareColor,
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
Fs = []
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
    F = visual.TextStim(
        win=myWin,
        color='white',
        height=circle_size-0.2,
        text='F',
        pos=pos
        )
    Hs.append(H)
    Vs.append(V)
    Fs.append(F)
    
triggerText = visual.TextStim(
    win=myWin,
    color='white',
    height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text='Experiment will start soon. Waiting for scanner'
    )

instructText = visual.TextStim(
    win=myWin, 
    color='white',
    height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text=f'Press {vertical_buttom} when you perceive VERTICAL\n\
        Press {horizontal_buttom} when you perceive HORIZONTAL\n\
            PRESS on the keys to practice'
    )
instruct_ITI = visual.TextStim(
    win=myWin,
    color='white',
    height=0.5,
    pos=apply_global_offset(base_pos=(0,-3.5),global_offset=global_offset),
    text=f'Press {ITI_buttom} when you perceive Four Squares FLASHING\n\
        Press the Key to Practice'
    )
anykeyText = visual.TextStim(
    win=myWin, 
    color='white',
    height=0.5,
    pos=apply_global_offset(base_pos=(0,0), global_offset=global_offset),
    text='Good Job!\n\
        Press on any key to continue'
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
    
    # Display instruction for flashing four dots
    while not event.getKeys(keyList=[ITI_buttom]):
        flickerSl(HoriDist, VertiDist, instruct=True) # this includes the instruction text
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
edf_file =  f"am_run{expInfo['run']}.EDF"
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

# %% RENDER_LOOP
# give the system time to settle
core.wait(1)

# Only shows instruction when 
if expInfo['run'] == '1':

    buttom_instruct(myWin,vertical_buttom, horizontal_buttom, ITI_buttom)
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
trigCount = 1     # counter triggers

# reset clocks
clock.reset()   # Comment out only for the simulation
logging.data('StartOfRun' + str(expInfo['run']))
logging.data(msg='Scanner trigger %i' % (trigCount))

el_tracker.sendMessage(f"EXPERIMENT_START {expInfo['expName']}")

while trigCount < totalTrigger:   # 

    logging.data('StartOfCondition'+ str(Conditions[i]))

    while trigCount < np.sum(Durations[0:i+1]):
        t = clock.getTime()
        
        # Log start of each TR and conditions to EYETRACKER
        #el_tracker.sendMessage(f"CONDITION_START Condition {Conditions[i]} Trigger {trigCount}")

        if Conditions[i] == 0:
            fixation()
        elif Conditions[i] == 1:
            flickerSl(HoriDist, VertiDist, instruct=False) #turn instruction mode off
        elif Conditions[i] == 2:
            quartet(HoriDist, VertiDist)

        for key in event.getKeys():   # e vuoto al primo trigger
                if key in ['escape', 'q']:
                    logging.data(msg='User pressed quit')
                    el_tracker.sendMessage("EXPERIMENT_ABORTED") # eyetracker log quit
                    el_tracker.stopRecording()
                    terminate_task()
                    myWin.close()
                    core.quit()
                elif key in ['1', 'num_1']:
                    t = clock.getTime()
                    KeyPressed = '1'
                    KeyPressedNew = np.array([KeyPressed, t])
                    KeyPressedArray = np.vstack((KeyPressedArray,
                                                 KeyPressedNew))
                    logging.data(msg='Key1 pressed')
                elif key in ['2', 'num_2']:
                    t = clock.getTime()
                    KeyPressed = '2'
                    KeyPressedNew = np.array([KeyPressed, t])
                    KeyPressedArray = np.vstack((KeyPressedArray,
                                                 KeyPressedNew))
                    logging.data(msg='Key2 pressed')
                elif key in ['3', 'num_3']:
                    t = clock.getTime()
                    KeyPressed = '3'
                    KeyPressedNew = np.array([KeyPressed, t])
                    KeyPressedArray = np.vstack((KeyPressedArray,
                                                 KeyPressedNew))
                    logging.data(msg='Key3 pressed')
                    
                elif key in ['4', 'num_4']:
                    t = clock.getTime()
                    KeyPressed = '4'
                    KeyPressedNew = np.array([KeyPressed, t])
                    KeyPressedArray = np.vstack((KeyPressedArray,
                                                 KeyPressedNew))
                    logging.data(msg='Key4 pressed')
                
                    
                elif key == TRIGGERKEY:
                    t = clock.getTime()
                    trigCount = trigCount+1
                    el_tracker.sendMessage(f"TRIGGER {trigCount}") # log trigger to eyetracker 
                    logging.data(msg='Scanner trigger %i' % (trigCount))
    i = i+1

logging.data('EndOfRun' + str(expInfo['run']) + '\n')

# %% time stamp ending EYETRACKER
el_tracker.sendMessage("EXPERIMENT_END")
#el_tracker.stopRecording()
#myWin.close()

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

# Change into protocol folder
os.chdir(parentDir)
os.chdir(prtFolderName)
## construct protocol file for AMB
# set key KeyPressedArray to pd.DataFrame
KeyPressed_df = pd.DataFrame(KeyPressedArray[1:], columns=KeyPressedArray[0])
# create timestamp as stop time for events 
KeyPressed_df['Timestamp'] = KeyPressed_df['KeyPressedt'].shift(-1)
# define flicker block ends
flicker_ends = [116, 212, 308, 404, 500, 596]
# find all flickerSl rows
flicker_indices = KeyPressed_df[KeyPressed_df['Label'] == 'flickerSl'].index
# loop through and assign the fixed values
for i, flicker_end in zip(flicker_indices, flicker_ends):
    flicker_start = flicker_end - 16
    KeyPressed_df.at[i, 'Timestamp'] = flicker_end #assign flicker end
    KeyPressed_df.at[i-1, 'KeyPressedt'] = flicker_start # assign flicker start a block before
    KeyPressed_df.at[i, 'KeyPressedt'] = flicker_start

# Add fixation row at the top and bottom
fixation_start = pd.DataFrame({
    'KeyPressed': [0],
    'KeyPressedt': [0.0],
    'Label': ['fixation'],
    'Timestamp': [20.0]
})

fixation_end = pd.DataFrame({
    'KeyPressed': [0],
    'KeyPressedt': [596.0],
    'Label': ['fixation'],
    'Timestamp': [616.0]
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

# EYETRACKER CLOSE DISPLAY AND SAVE EDF
os.chdir(parentDir)
el_tracker.stopRecording()
terminate_task()
myWin.close()
# %% FINISH#
core.quit()