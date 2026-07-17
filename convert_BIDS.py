from pathlib import Path
import pandas as pd
import re
import numpy as np
# ============================================================
# Root directory
# ============================================================
root = Path("/Users/siwentao/Documents/Github/7Ttesting ")
# Output BIDS directory
bids_root = f"{root}/BIDS_events"
# Mapping from original folder names to BIDS subject IDs
subject_map = {
    "NH_SubjData": "sub-01",
    "EC_SubjData": "sub-02",
}
# Rename stimulus labels
condition_mapping = {
    "fixation": "fixation",
    "vertiM": "vertical_motion",
    "horiM": "horizontal_motion",
    "flickerSI": "flicker_static",
}

# ============================================================
# Convert each subject
# ============================================================
for subj_folder, bids_sub in subject_map.items():
    # --------------------------------------------------------
    # hMT loc
    # --------------------------------------------------------
    TR = 2.0 # hMT loc 2 sec TR
    loc_condition_mapping = {0: "stationary", 1: "moving"}
    loc_protocol_dir = Path(root, subj_folder, "hMT+localizer", "Protocols")
    output_dir = Path(bids_root, bids_sub, "func")
    output_dir.mkdir(parents=True, exist_ok=True)
    loc_protocol_files = list(loc_protocol_dir.glob("*.csv"))
    loc_files_by_run = {}
    
    for loc_input_file in loc_protocol_files:
        match = re.search(r"run(\d+)\.csv$", loc_input_file.name)
        if match is None:
            continue
        run_num = int(match.group(1))
        loc_files_by_run[run_num] = loc_input_file
    
        # Convert each localizer run
    for run_num, loc_input_file in sorted(loc_files_by_run.items()):
        loc_df = pd.read_csv(loc_input_file)
        # Remove accidental unnamed index columns
        loc_df = loc_df.loc[:, ~loc_df.columns.str.startswith("Unnamed")].copy()
        # Make sure required columns exist
        required_columns = {
            "Condition",
            "Duration",
            "Target",
        }
        missing_columns = required_columns - set(loc_df.columns)
        if missing_columns:
            raise ValueError(
                f"{loc_input_file.name} is missing columns: "
                f"{sorted(missing_columns)}"
            )
        # Convert values to numeric
        loc_df["Condition"] = pd.to_numeric(loc_df["Condition"], errors="raise")
        loc_df["Duration"] = pd.to_numeric(loc_df["Duration"], errors="raise")
        
        # ----------------------------------------------------
        # Block events
        # ----------------------------------------------------
        block_events = pd.DataFrame()
        # Recalculate onset from cumulative durations.
        # This also fixes the final row where Onset incorrectly returns to 0.
        block_events["onset"] = (loc_df["Duration"].cumsum().shift(fill_value=0))
        block_events["duration"] = loc_df["Duration"]
        block_events["trial_type"] = (loc_df["Condition"].map(loc_condition_mapping))
        # Check for unmapped condition values
        if block_events["trial_type"].isna().any():
            unknown_conditions = sorted(
                loc_df.loc[
                    block_events["trial_type"].isna(), "Condition",
                ].unique()
            )
            raise ValueError(
                f"Unknown localizer conditions in "
                f"{loc_input_file.name}: {unknown_conditions}"
            )

        # ----------------------------------------------------
        # Combine and sort all events
        # ----------------------------------------------------
        loc_events = (block_events.sort_values(["onset", "duration"]).reset_index(drop=True))
        # Ensure numeric column
        loc_events["onset"] = loc_events["onset"].astype(float)
        loc_events["duration"] = loc_events["duration"].astype(float)
        # multiply TR by TR LENGTH 
        loc_events["onset"] *= TR
        loc_events["duration"] *= TR

        # ----------------------------------------------------
        # Save BIDS events.tsv
        # ----------------------------------------------------
        loc_output_file = (output_dir/ f"{bids_sub}_task-hMTlocalizer_run-{run_num:02d}_events.tsv")
        loc_events.to_csv(
            loc_output_file,
            sep="\t",
            index=False,
            float_format="%.3f",
            na_rep="n/a",
        )
        print(f"Saved {loc_output_file}")
    
    # --------------------------------------------------------
    # PHY
    # --------------------------------------------------------
    TR = 1.612 # the rest of high res has 1.612 TR
    phy_protocol_dir = Path(root, subj_folder, "Phy_MotQuart", "Protocols")
    phy_protocol_files = list(phy_protocol_dir.glob("*.csv"))

    if len(phy_protocol_files) != 1:
        print(f"Expected one protocol in {phy_protocol_dir}, found {len(phy_protocol_files)}")
        continue

    phy_input_file = phy_protocol_files[0]
    # Skip if protocol doesn't exist
    if not phy_input_file.exists():
        print(f"Skipping {subj_folder}: protocol not found.")
        continue
    
    phy_runs = 6 # 6 physical runs
    phy_df = pd.read_csv(phy_input_file)
    phy_events = (phy_df[["Onset", "Durations", "Stim"]].rename(
            columns={
                "Onset": "onset",
                "Durations": "duration",
                "Stim": "trial_type",
            }))

    # Convert onset and duration from TRs to seconds
    phy_events["onset"] = (pd.to_numeric(phy_events["onset"], errors="raise") * TR)
    phy_events["duration"] = (pd.to_numeric(phy_events["duration"], errors="raise") * TR)

    phy_events["trial_type"] = (
        phy_events["trial_type"]
        .astype(str)
        .str.strip()
        .replace(condition_mapping)
    )

    # SAVE Physical event files
    for run_num in range(1, phy_runs+1):
        phy_output_file = (output_dir/f"{bids_sub}_task-physical_run-{run_num:02d}_events.tsv")
        # Save as tab-separated file
        phy_events.to_csv(
            phy_output_file,
            sep="\t",
            index=False,
            float_format="%.3f",
        )
        print(f"Saved {phy_output_file}")
              
    # --------------------------------------------------------
    # AMB
    # --------------------------------------------------------
    amb_run_TR = 248
    amb_fix_TR = 8
    amb_run_duration = amb_run_TR * TR
    amb_fixation_duration = amb_fix_TR * TR
    amb_final_fixation_onset = amb_run_duration - amb_fixation_duration
    amb_input_dir = (root/ subj_folder/ "Amb_MotQuart"/ "Output")
    amb_protocol_files = list(amb_input_dir.glob("*_amb_run*_key_presses.csv"))
    amb_files_by_run = {}

    for amb_input_file in amb_protocol_files:
        match = re.search(r"_amb_run(\d+)_key_presses\.csv$",amb_input_file.name,)
        if match is None:
            continue
        run_num = int(match.group(1))
        amb_files_by_run[run_num] = amb_input_file
    # --------------------------------------------------------
    # Process AMB runs
    # --------------------------------------------------------
    for run_num, amb_input_file in sorted(amb_files_by_run.items()):
        print(
            f"Reading {bids_sub} ambiguous run "
            f"{run_num:02d}: {amb_input_file.name}"
        )
        amb_df = pd.read_csv(amb_input_file)
        amb_output_file = (output_dir / f"{bids_sub}_task-ambiguous_run-{run_num:02d}_events.tsv")

        # Build BIDS events
        amb_events = (
            amb_df[["KeyPressedt", "Label"]].copy().rename(
                columns={
                    "KeyPressedt": "onset",
                    "Label": "trial_type",
                }
            )
        )
        amb_events["onset"] = pd.to_numeric(amb_events["onset"],errors="raise",)
        amb_events = (amb_events.sort_values("onset").reset_index(drop=True))
        amb_events["trial_type"] = (amb_events["trial_type"].astype(str).str.strip().replace(condition_mapping))

        # Duration = next onset - current onset
        amb_events["duration"] = (amb_events["onset"].shift(-1) - amb_events["onset"])
        # Last logged event ends when final fixation begins
        amb_events.loc[amb_events.index[-1],"duration",] = (amb_final_fixation_onset - amb_events.loc[amb_events.index[-1],"onset",])

        # Add initial fixation
        initial_fixation = pd.DataFrame(
            {
                "onset": [0.0],
                "duration": [amb_fixation_duration],
                "trial_type": ["fixation"],
            }
        )
        # Add final fixation
        final_fixation = pd.DataFrame(
            {
                "onset": [amb_final_fixation_onset],
                "duration": [amb_fixation_duration],
                "trial_type": ["fixation"],
            }
        )
        amb_events = pd.concat([initial_fixation, amb_events,final_fixation,],ignore_index=True,)
        # Sort once more
        amb_events = (amb_events.sort_values("onset").reset_index(drop=True))

        # Validate
        if (amb_events["duration"] <= 0).any():
            bad_rows = amb_events.loc[amb_events["duration"] <= 0]
            raise ValueError(f"Negative durations found:\n{bad_rows}")
        # Save
        amb_events.to_csv(
            amb_output_file,
            sep="\t",
            index=False,
            float_format="%.3f",
        )
        print(f"Saved {amb_output_file}")
        
    # --------------------------------------------------------
    # VOL
    # --------------------------------------------------------
    vol_run_TR = 210
    vol_run_duration = amb_run_TR * TR
    vol_input_dir = (root/ subj_folder/ "Volitional_MotQuart"/ "Protocols")
    vol_protocol_files = list(vol_input_dir.glob("*.csv"))
    vol_files_by_run = {}
    
    for vol_input_file in vol_protocol_files:
        match = re.search(r"_Volitional_Run(\d+)_protocol.csv$",vol_input_file.name,)
        if match is None:
            continue
        run_num = int(match.group(1))
        vol_files_by_run[run_num] = vol_input_file
        
    # -----------------------------------------------------------
    # Process VOL
    # -----------------------------------------------------------
    for run_num, vol_input_file in sorted(vol_files_by_run.items()):

        print(f"Reading {bids_sub} volitional run {run_num:02d}: {vol_input_file.name}")
        vol_df = pd.read_csv(vol_input_file)

        # Clean accidental whitespace from column names
        vol_df.columns = vol_df.columns.str.strip()
        print("Columns:", vol_df.columns.tolist())
        print(vol_df.head())
        vol_output_file = (output_dir/ f"{bids_sub}_task-volitional_run-{run_num:02d}_events.tsv")

        # Validate required columns
        required_columns = [
            "Trial",
            "Stim",
            "Duration",
            "Onset",
            "Instruct_V_H",
            "QuartetOrder",
            "illusory_physical",
            "ResponseKey",
            "ResponseTime",
            "expected_key",
            "success",
        ]

        missing_columns = [
            column
            for column in required_columns
            if column not in vol_df.columns
        ]

        if missing_columns:
            raise ValueError(
                f"{vol_input_file.name} is missing columns: "
                f"{missing_columns}\n"
                f"Available columns: {vol_df.columns.tolist()}"
            )

        # Keep only PrecueTime and DelayTime
        vol_events = vol_df.loc[vol_df["Stim"].isin(["PrecueTime", "DelayTime",])].copy()
        if vol_events.empty:
            raise ValueError(
                f"No PrecueTime or DelayTime rows found in "
                f"{vol_input_file.name}"
            )

        # Convert timing from TRs to seconds
        vol_events["Onset"] = pd.to_numeric(vol_events["Onset"], errors="raise",)
        vol_events["Duration"] = pd.to_numeric(vol_events["Duration"],errors="raise",)
        vol_events["onset"] = vol_events["Onset"] * TR
        vol_events["duration"] = vol_events["Duration"] * TR

        # Clean phase and instruction labels
        stim_mapping = {"PrecueTime": "precue", "DelayTime": "delay",}
        vol_events["phase"] = (vol_events["Stim"].astype(str).str.strip().replace(stim_mapping))
        vol_events["instructed_axis"] = (vol_events["Instruct_V_H"].astype(str).str.strip().str.lower())

        # trial_type examples:
        # precue_horizontal
        # delay_horizontal
        # precue_vertical
        # delay_vertical

        vol_events["trial_type"] = (vol_events["phase"]+ "_" + vol_events["instructed_axis"])

        # Preserve trial metadata
        vol_events["trial"] = pd.to_numeric(vol_events["Trial"], errors="raise",).astype(int)
        vol_events["quartet_order"] = (vol_events["QuartetOrder"].astype(str).str.strip())
        vol_events["stimulus_type"] = (vol_events["illusory_physical"].astype(str).str.strip().str.lower())

        # Convert success to consistent lowercase text
        # BIDS permits extra columns containing strings.
        vol_events["success"] = (
            vol_events["success"]
            .astype(str)
            .str.strip()
            .str.lower()
            .replace(
                {
                    "true": "1",
                    "false": "0",
                }
            )
        )
        # Preserve response information
        vol_events["response_key"] = (vol_events["ResponseKey"].astype(str).str.strip())
        vol_events["expected_key"] = (vol_events["expected_key"].astype(str).str.strip())
        vol_events["response_time"] = pd.to_numeric(vol_events["ResponseTime"],errors="coerce",)

        # The logged ResponseTime appears to be measured from the
        # beginning of the run. Keep it as response_onset rather
        # than calling it reaction time.
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
                vol_events[bids_column] = vol_events[
                    original_column
                ]

        # Select and order output columns
        output_columns = [
            "onset",
            "duration",
            "trial_type",
            "success",
        ]
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
            "expected_key",
        ]

        output_columns.extend(
            column
            for column in optional_output_columns
            if column in vol_events.columns
        )
        '''
        
        vol_events = vol_events[output_columns]
        
        # Sort and validate
        vol_events = (
            vol_events
            .sort_values(
                [
                    "onset",
                    "trial_type",
                ]
            )
            .reset_index(drop=True)
        )

        if vol_events["onset"].isna().any():
            raise ValueError(
                f"Missing onset values in {vol_input_file.name}"
            )

        if vol_events["duration"].isna().any():
            raise ValueError(
                f"Missing duration values in {vol_input_file.name}"
            )

        if (vol_events["onset"] < 0).any():
            raise ValueError(
                f"Negative onset found in {vol_input_file.name}"
            )

        if (vol_events["duration"] <= 0).any():
            bad_rows = vol_events.loc[
                vol_events["duration"] <= 0
            ]

            raise ValueError(
                f"Non-positive durations found in "
                f"{vol_input_file.name}:\n"
                f"{bad_rows}"
            )

        event_ends = (vol_events["onset"] + vol_events["duration"])

        if (event_ends > vol_run_duration + 1e-6).any():
            bad_rows = vol_events.loc[event_ends > vol_run_duration + 1e-6]

            raise ValueError(
                f"Events extend beyond the expected run duration "
                f"of {vol_run_duration:.3f} seconds in "
                f"{vol_input_file.name}:\n"
                f"{bad_rows}"
            )

        # Save BIDS events.tsv
        vol_events.to_csv(
            vol_output_file,
            sep="\t",
            index=False,
            na_rep="n/a",
            float_format="%.3f",
        )
        print(f"Saved {vol_output_file}")
        print(vol_events.head())
        
        
        
        