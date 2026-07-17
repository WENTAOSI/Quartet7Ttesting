from psychopy import visual, core, event, monitors, gui
import numpy as np
import random
import csv
from pathlib import Path
import pandas as pd

##########################For Vandy 7T##########################################
# set global offset
ho_dva = -0.0981
vo_dva = 1.7652
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
# GUI Store info about experiment and experimental run
expName = 'MT_IPS_loc'  # set experiment name here
expInfo = {
    'run': '1',
    'sub-num': 'sub-00',
    'type': ['stationarybaseline', 'movingbaseline'],
    'Eyelink':['False'],
    'display': ['Vanderbilt7T'],
    'TR': ['2']
    }
# Create GUI at the beginning of exp to get more expInfo
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False: core.quit()  # user pressed cancel

# =====================================================
# PARAMETERS
# =====================================================
TR = float(expInfo['TR'][0])  # Use the selected TR value
trigger_key = "quoteleft"
quit_key = 'escape'

n_blocks = 14           # MOT blocks
cue_TRs = 1            # 2 s target cue
track_TRs = 5          # 10 s tracking
response_max = 2       # 4 s response period
base_TRs = 5           # 10 s baseline

n_objects = 10         # 5 left, 5 right
dot_radius = 0.1
speed = 2.0

# Foveal circular aperture
aperture_radius = 4.0   # dva radius
midline_gap = 0.3       # invisible gap around vertical meridian
# Keep the whole dot inside the circular aperture
effective_radius = aperture_radius - dot_radius

# =====================================================
# DEFINE SCANNER DISPLAY
# =====================================================
if expInfo['display'] == 'Vanderbilt7T':
    scanner_monitor = monitors.Monitor("scanner")
    # Physical width of display in cm
    scanner_monitor.setWidth(17.0)
    # Viewing distance from eyes to screen via mirror, in cm
    scanner_monitor.setDistance(48.0)
    # Projector resolution
    scanner_monitor.setSizePix([1024, 768])
    scanner_monitor.save()
else:
    raise ValueError(f"Invalid display input: {expInfo['display']}. ")

# =====================================================
# WINDOW
# =====================================================
win = visual.Window(
    monitor=scanner_monitor,
    fullscr=True,
    units="deg",
    color="black",
    checkTiming=False
)
fix = visual.TextStim(
    win,
    text="+",
    pos=apply_global_offset((0, 0)),
    color="white",
    height=0.35
)
dots = [
    visual.Circle(
        win,
        radius=dot_radius,
        fillColor="white",
        lineColor="white"
    )
    for _ in range(n_objects)
]
probe_ring = visual.Circle(
    win,
    radius=dot_radius * 1.8,
    lineColor="red",
    fillColor=None,
    lineWidth=3
)
text_stim = visual.TextStim(
    win,
    text="",
    color="white",
    height=0.4
)

# =====================================================
# FUNCTIONS
# =====================================================
def wait_trs(n_trs):
    count = 0
    while count < n_trs:
        keys = event.waitKeys(keyList=[trigger_key, quit_key])
        if quit_key in keys:
            win.close()
            core.quit()
        if trigger_key in keys:
            count += 1

def is_left(i):
    """
    Dots 0-4 remain in the left visual field.
    Dots 5-9 remain in the right visual field.
    """
    return i < n_objects // 2

def generate_positions():
    """
    Generate starting positions uniformly within the circular aperture.

    The first five dots are restricted to the left hemifield.
    The last five dots are restricted to the right hemifield.
    """
    positions = []

    for i in range(n_objects):
        while True:
            # sqrt gives approximately uniform spatial density
            # over the area of the circle
            radius = effective_radius * np.sqrt(np.random.random())
            angle = np.random.uniform(0, 2 * np.pi)
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            if is_left(i):
                if x <= -midline_gap / 2:
                    break
            else:
                if x >= midline_gap / 2:
                    break
        positions.append([x, y])
    return np.array(positions, dtype=float)

def generate_velocities():
    angles = np.random.uniform(0, 2 * np.pi, n_objects)
    return np.column_stack([np.cos(angles) * speed, np.sin(angles) * speed])

def update_positions(pos, vel, dt):
    """
    Update dot positions while enforcing:
    1. A circular aperture with radius 4 dva.
    2. An invisible vertical boundary at fixation.
    3. Five dots in each hemifield.
    """
    new_pos = pos + vel * dt
    
    for i in range(n_objects):
        # -------------------------------------------------
        # Vertical meridian boundary
        # -------------------------------------------------
        if is_left(i):
            left_boundary = -midline_gap / 2
            if new_pos[i, 0] > left_boundary:
                new_pos[i, 0] = left_boundary
                vel[i, 0] *= -1
        else:
            right_boundary = midline_gap / 2
            if new_pos[i, 0] < right_boundary:
                new_pos[i, 0] = right_boundary
                vel[i, 0] *= -1

        # -------------------------------------------------
        # Circular aperture boundary
        # -------------------------------------------------
        distance_from_center = np.linalg.norm(new_pos[i])
        
        if distance_from_center > effective_radius:
            # Normal vector pointing outward from the circle
            normal = (new_pos[i] / distance_from_center)
            # Reflect velocity across the circular boundary
            vel[i] = (vel[i] - 2 * np.dot(vel[i], normal) * normal)
            # Place dot back on the valid boundary
            new_pos[i] = (normal * effective_radius)

        # Recheck hemifield after circular reflection
        if is_left(i):
            left_boundary = -midline_gap / 2
            if new_pos[i, 0] > left_boundary:
                new_pos[i, 0] = left_boundary
                vel[i, 0] = -abs(vel[i, 0])
        else:
            right_boundary = midline_gap / 2
            if new_pos[i, 0] < right_boundary:
                new_pos[i, 0] = right_boundary
                vel[i, 0] = abs(vel[i, 0])
    return new_pos, vel

def draw_fixation():
    fix.draw()

def draw_dots(pos, target_indices=None,show_targets=False):
    for i, dot in enumerate(dots):
        dot.pos = apply_global_offset(pos[i])
        if (
            show_targets
            and target_indices is not None
            and i in target_indices
        ):
            dot.fillColor = "red"
            dot.lineColor = "red"
        else:
            dot.fillColor = "white"
            dot.lineColor = "white"
        dot.draw()
    draw_fixation()


def choose_bilateral_targets():
    """
    Choose exactly one target from each hemifield.
    """
    left_indices = list(range(n_objects // 2))
    right_indices = list(range(n_objects // 2, n_objects))
    return [random.choice(left_indices), random.choice(right_indices)]

def stationary_period(pos, n_trs):

    trigger_count = 0
    while trigger_count < n_trs:
        draw_dots(pos)
        win.flip()
        keys = event.getKeys(keyList=[trigger_key, quit_key])
        if quit_key in keys:
            win.close()
            core.quit()
        if trigger_key in keys:
            trigger_count += 1

def run_motion_period(pos, vel, n_trs):
    clock = core.Clock()
    last_t = clock.getTime()

    for tr in range(n_trs):
        got_trigger = False
        while not got_trigger:
            now = clock.getTime()
            dt = now - last_t
            last_t = now
            pos, vel = update_positions(pos, vel, dt)
            draw_dots(pos)
            win.flip()
            keys = event.getKeys(keyList=[trigger_key, quit_key])

            if quit_key in keys:
                win.close()
                core.quit()
            if trigger_key in keys:
                got_trigger = True
    return pos, vel

def response_period(pos, target_indices, probe_condition):

    if probe_condition == "target":
        probe_index = random.choice(target_indices)
    else:
        distractors = [
            i for i in range(n_objects)
            if i not in target_indices
        ]
        probe_index = random.choice(distractors)

    correct_answer = (
        "1"
        if probe_index in target_indices
        else "2"
    )

    resp_clock = core.Clock()
    response = None
    rt = None

    while resp_clock.getTime() < response_max:
        draw_dots(pos)
        probe_ring.pos = apply_global_offset(pos[probe_index])
        probe_ring.draw()
        text_stim.text = (
            "Was the highlighted dot a target?\n\n"
            "1 = yes    2 = no"
        )
        text_stim.pos = apply_global_offset((0, -4.5))
        text_stim.draw()
        win.flip()
        keys = event.getKeys(keyList=["1", "2", quit_key], timeStamped=resp_clock)

        for key, t in keys:
            if key == quit_key:
                win.close()
                core.quit()
            if response is None:
                response = key
                rt = t
    accuracy = (
        int(response == correct_answer)
        if response is not None
        else 0
    )
    return probe_index, correct_answer, response, rt, accuracy

# =====================================================
# EXPERIMENT
# =====================================================
results = []

probe_conditions = (["target"] * 9 +["distractor"] * 9)
random.shuffle(probe_conditions)

text_stim.text = (
    "Foveal bilateral MOT localizer\n\n"
    "Track the two red dots.\n"
    "One target stays left of fixation.\n"
    "One target stays right of fixation.\n\n"
    "Waiting for scanner trigger..."
)

text_stim.pos = apply_global_offset((0,0))
text_stim.draw()
win.flip()

# First scanner pulse starts the task
wait_trs(1)

# initialization 
bids_events = []
current_onset = 0.0

for block_num in range(1, n_blocks + 1):

    if block_num == 1:
        pos = generate_positions()
        vel = generate_velocities()
    # -----------------------------
    # FIXATION BASELINE: 5 TRs Stationary or Moving dots 
    # -----------------------------
    if expInfo['type'] == 'stationary_baseline':
        stationary_period(pos, base_TRs)
    elif expInfo['type'] == 'moving_baseline':
        pos, vel = run_motion_period(pos, vel, base_TRs)
    
    target_indices = (choose_bilateral_targets())
    fixation_duration = base_TRs * TR
    
    bids_events.append({
        "onset": current_onset,
        "duration": fixation_duration,
        "trial_type": expInfo['type'],
        "block": block_num,
        "probe_condition": "n/a",
        "response": "n/a",
        "response_time": "n/a",
        "accuracy": "n/a"
    })
    current_onset += fixation_duration
    
    # -----------------------------
    # TARGET CUE: 1 TR
    # -----------------------------
    draw_dots(pos, target_indices, show_targets=True)
    win.flip()
    wait_trs(cue_TRs)
    cue_duration = cue_TRs * TR
    
    bids_events.append({
        "onset": current_onset,
        "duration": cue_duration,
        "trial_type": "target_cue",
        "block": block_num,
        "probe_condition": "n/a",
        "response": "n/a",
        "response_time": "n/a",
        "accuracy": "n/a"
    })
    current_onset += cue_duration
    
    # -----------------------------
    # TRACKING: 5 TRs
    # -----------------------------
    pos, vel = run_motion_period(pos, vel, track_TRs)
    tracking_duration = track_TRs * TR
    bids_events.append({
        "onset": current_onset,
        "duration": tracking_duration,
        "trial_type": "tracking",
        "block": block_num,
        "probe_condition": "n/a",
        "response": "n/a",
        "response_time": "n/a",
        "accuracy": "n/a"
    })
    current_onset += tracking_duration

    # -----------------------------
    # RESPONSE: maximum 4 seconds
    # -----------------------------
    probe_condition = probe_conditions[block_num - 1]
    probe_index, correct_answer, response, rt, accuracy = response_period(
        pos,
        target_indices,
        probe_condition
    )
    response_max_duration = response_max * TR

    bids_events.append({
        "onset": current_onset,
        "duration": response_max_duration,
        "trial_type": "response",
        "block": block_num,
        "probe_condition": probe_condition,
        "response": response if response is not None else "n/a",
        "response_time": rt if rt is not None else "n/a",
        "accuracy": accuracy
    })
    current_onset += response_max_duration

# =====================================================
# SAVE
# =====================================================
sub_num = str(expInfo["sub-num"]).strip()
run = int(expInfo["run"])
output_dir = Path("BIDS_events", sub_num, "func")
output_dir.mkdir(parents=True, exist_ok=True)
output_file = (output_dir / f"{sub_num}_task-mtipsloc{expInfo['type']}_run-{run:02d}_events.tsv")

events_df = pd.DataFrame(bids_events)
events_df.to_csv(
    output_file,
    sep="\t",
    index=False,
    float_format="%.3f",
    na_rep="n/a"
)
print(f"Saved {output_file}")
print(f"Protocol duration: {current_onset:.1f} seconds")

# =====================================================
# END
# =====================================================
text_stim.text = "Done."
text_stim.pos = apply_global_offset((0,0))
text_stim.draw()
win.flip()
core.wait(2)
win.close()
core.quit()