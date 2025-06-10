"""Prepare condition order, targets and noise texture for stim presentation."""

import numpy as np
import os
import csv

runs = [1]
##############################3
TR = 2
period = 10 #sec 
###############################

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


# Save the conditions as .csv
filename_csv = os.path.join(f'MtLoc_MQ_Condition.csv')

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
