#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import os
import sys
import shutil 

# --- Configuration ---
# Define the paths and parameters (makes it easier to modify)

# Path to the executable script relative to this Python script's location
# Adjust if your script structure is different.
executable_path = "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/dl_binder_design/mpnn_fr/dl_interface_design.py"

# Input/Output directories relative to this Python script's location
pdb_input_dir = "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Backbone_Gen/outputs"
temp_dir = "/mnt/c/pdb_temp_work"
output_dir = "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Sequence_Gen_output"

# Parameters for the script
relax_cycles = 1
seqs_per_struct = 1

# --- End Configuration ---


def sequence_design():
    """Runs ProteinMPNN using the PDB files (input) corresponding to the backbones designed by RFdiffusion."""

    # Remove the checkpoint file (from a previous run)
    checkpoint_file = "check.point" # File to remove before starting
    try:
        os.remove(checkpoint_file)
    except OSError as e:
        print(f"Error removing '{checkpoint_file}': {e}", file=sys.stderr)
 

    # Remove previous temporary directory (if it exists)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)

    # Copy PBD files to a temporary directory (I need to do this because OneDrive gives errors)
    if not os.path.isdir(pdb_input_dir):
            raise FileNotFoundError(f"Source directory not found: {pdb_input_dir}")
    else:
        shutil.copytree(pdb_input_dir, temp_dir)


    # Construct the command arguments as a list
    command_to_run_in_conda = [
        "python", # Use the python from the 'binder_design' env
        executable_path,
        "-pdbdir", temp_dir,
        "-relax_cycles", str(relax_cycles),
        "-seqs_per_struct", str(seqs_per_struct),
        "-outpdbdir", output_dir
    ]

    # Full command for subprocess
    conda_run_command = ["conda",
        "run",
        "-n", "binder_design", #Name of the conda environment
        "--cwd", temp_dir, # Run the command with temp_work_dir as CWD
        *command_to_run_in_conda # Unpack the command list
    ]

    # Run the command using subprocess
    try:
        result = subprocess.run(
            conda_run_command,
            check=True, #Will raise an exception if something fails
            capture_output=True, #Captures the standard output and standard error.
            text=True #The captured output and errors are converted to strings.
        )

        print("\ProteinMPNN executed successfully.")
        print("--- Standard Output ---")
        print(result.stdout)
        print("-----------------------")
        if result.stderr:
            print("--- Standard Error (Warnings/Info) ---")
            print(result.stderr)
            print("------------------------------------")

    except Exception as error:
        print(f"\nAn unexpected error occurred: {error}", file=sys.stderr)
        sys.exit()

    finally:
        # Always attempt to remove the temporary directory
        if os.path.exists(temp_dir):
            try:
                shutil.rmtree(temp_dir)
            except OSError as error:
                print(f"Error removing temporary directory {temp_dir}: {error}")

# Standard Python practice: only run the main logic if the script is executed directly
if __name__ == "__main__":
    sequence_design()