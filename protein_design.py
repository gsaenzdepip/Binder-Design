import logging
import docker
import requests
import sys
import subprocess
# from pymol import cmd #Needs to activate 'pymol' environment in conda.
import os
import sys
import shutil 
import re
import __main__


log_filename  = "logging_info.log"
log_directory = "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python"
log_filepath = os.path.join(log_directory, log_filename)

logging.basicConfig(
    level=logging.INFO,  # Log everything from DEBUG upwards
    format='%(asctime)s - %(levelname)s - %(message)s',
    filename=log_filepath,
    filemode='w' # overwrites the file each time
)

logger = logging.getLogger()

class Backbone_Gen():
  def __init__(self, input_folder = "", output_folder = "", models_folder = "", pdb_id = "", target_id = ""):
    self.inputs = input_folder
    self.outputs = output_folder
    self.models = models_folder
    self.pdb_id = pdb_id
    self.target_id = target_id

  def get_pdb_file(self):
    '''
    Task:
      - Download the PBD file that will be used for the conditional binder generation.
    '''
    url = f"https://files.rcsb.org/view/{self.pdb_id}.pdb"
    try:
      response = requests.get(url)
      response.raise_for_status()
      with open(f'{self.inputs}/{self.pdb_id}.pdb', 'wb') as file:
        file.write(response.content)
      print("\nPDB file downloaded.\n")
    except requests.exceptions.RequestException as e:
      print(f"\nFailed to download PDB file: {e}\n")
      sys.exit()

  def remove_non_protein(self):
      """
      Removes everything except protein (amino acids) from the loaded structure
      and saves the cleaned structure as '<pdb_id>_cleaned.pdb'.

      IMPORTANT: Needs to activate 'pymol' environment in conda.

      :param pdb_id: The base PDB identifier (string) used in the saved file name.

      List of methods for Pymol cmd: https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha
      """
      # First, load your PDB file (if not already loaded)
      cmd.load(f'{self.inputs}/{self.pdb_id}.pdb', f'{self.pdb_id}')
      # Remove all non-protein molecules
      cmd.remove("not polymer.protein")
      cmd.rebuild()  # Update the view if necessary

      # Construct the filename and save the cleaned PDB
      filename = f"{self.pdb_id}_cleaned.pdb"
      cmd.save(f'{self.inputs}/{filename}')
      print(f"\nSaved cleaned structure as: {filename}\n")

  def chain_residue_ranges(self):
      """
      Shows the different chains and residue ranges in the protein of the loaded and cleaned PDB file.

      IMPORTANT: Needs to activate 'pymol' environment in conda.

      """
      cmd.load(f'{self.inputs}/{self.pdb_id}_cleaned.pdb', f'{self.pdb_id}')

      __main__._residues_by_chain = {}

      # Iterate over all atoms.
      # If 'resi' is numeric, convert to int; otherwise, leave it as a string.
      cmd.iterate("all",
                  "(__main__._residues_by_chain.setdefault(chain, set()).add(int(resi)) if resi.isdigit() "
                  "else __main__._residues_by_chain.setdefault(chain, set()).add(resi))")

      # For each chain, sort the residues and print the range.
      for chain, residues in __main__._residues_by_chain.items():
          try:
              sorted_residues = sorted(residues, key=lambda x: int(x) if isinstance(x, str) and x.isdigit() else x)
          except Exception:
              sorted_residues = sorted(residues)
          print(f"\nChain {chain}: residues {sorted_residues[0]}-{sorted_residues[-1]}")

  def backbone_gen(self, design_num = 20, hotspots = "", noise = 0, contigs = ""):
    '''
    Task:
      Designs backbones of binders against a target protein using RFdiffusion.
    '''
    try:
      logger.info("\nStarting the Backbone Design with RFdiffusion...\n")
      logger.info(f"\nConfiguration:\n"
                  f"Number of backbone designs: {design_num}\n"
                  f"Noise level: {noise}\n"
                  f"Contigs: {contigs}\n"
                  f"Hotspots: {hotspots}\n"
                  f"Target PDB ID: {self.pdb_id}\n"
                  f"Target name: {self.target_id}\n"
                  )

      client = docker.from_env()

      #Docker Volume Mounts and File Accessibility. Docker container does not have direct access to Windows file paths unless you mount the host directories into the container.
      volumes = {
      self.models: {'bind': '/models', 'mode': 'rw'},
      self.inputs: {'bind': '/inputs', 'mode': 'rw'},
      self.outputs: {'bind': '/outputs', 'mode': 'rw'},
      }

      # Define the command you want to run inside the container.
      # Note: Use Linux-style paths (the ones mounted in the container) for parameters.
      command = [
          f"inference.output_prefix=/outputs/{self.pdb_id}_design",
          f"inference.model_directory_path=/models",
          f"inference.input_pdb=/inputs/{self.pdb_id}_cleaned.pdb",
          f"inference.num_designs={design_num}",
          f"contigmap.contigs=[{contigs}]",
          f"ppi.hotspot_res=[{hotspots}]",
          f"denoiser.noise_scale_ca={noise}",
          f"denoiser.noise_scale_frame={noise}"
      ]

      # For GPU support, we can specify device_requests (Docker SDK version 4.4+)
      device_requests = [docker.types.DeviceRequest(count=-1, capabilities=[['gpu']])]

      # Run the container
      container = client.containers.run(
          image="rfdiffusion",    # Tag of the docker image
          command=command,
          volumes=volumes,
          detach=True,            # Run container in the background
          device_requests=device_requests,  # Enable GPU support
          # environment={"HYDRA_FULL_ERROR": "1"} # Allows to print errors
      )

      # print("Container started with ID:", container.id)

      # Optionally wait for the container to finish and then print the logs
      container.wait()
      logs = container.logs()
      print("Container logs:")
      print(logs.decode("utf-8"))

      container.remove()
      logger.info("\nRFdiffusion executed successfully.\n")

    except Exception as error:
       logger.info(f'An error happened during Backbone Generation: {error}')


class Sequence_Gen():

  def __init__(self, input_folder = "", output_folder = "", temp_dir = "", executable=""):
    self.inputs = input_folder
    self.outputs = output_folder
    self.temp_dir = temp_dir
    self.executable = executable

  def sequence_design(self, relax_cycles, seqs_per_struct, temperature):
      """
      Runs ProteinMPNN using the PDB files (input) corresponding to the backbones designed by RFdiffusion.
      """
      logger.info("\nStarting Sequence Design with ProteinMPNN...\n")
      logger.info(f"\nConfiguration:\n"
            f"Number of relax cycles: {relax_cycles}\n"
            f"Sequences per backbone: {seqs_per_struct}\n"
            f"Temperature: {temperature}\n"
            )
      
      # Remove the checkpoint file (from a previous run)
      checkpoint_file = "check.point" # File to remove before starting
      try:
          os.remove(checkpoint_file)
      except Exception as e:
          logger.info(f"'{checkpoint_file}': {e}")

      # Remove previous temporary directory (if it exists)
      if os.path.exists(self.temp_dir):
          shutil.rmtree(self.temp_dir)

      # Copy PBD files to a temporary directory (I need to do this because OneDrive gives errors)
      try:
        if not os.path.isdir(self.inputs):
                raise FileNotFoundError()
      except Exception as error:
         logger.info(f"Source directory with PDB files from RFdiffusion not found: {self.inputs}. {error}")
         sys.exit()
      
      shutil.copytree(self.inputs, self.temp_dir)

      # Construct the command arguments as a list
      command_to_run_in_conda = [
          "python", # Use the python from the 'binder_design' env
          self.executable, #Path to the ProteinMPNN executing .py file
          "-pdbdir", self.temp_dir,
          "-relax_cycles", str(relax_cycles),
          "-seqs_per_struct", str(seqs_per_struct),
          "-temperature", str(temperature),
          "-outpdbdir", self.outputs
      ]

      # Full command for subprocess
      conda_run_command = ["conda",
          "run",
          "-n", "binder_design", #Name of the conda environment
          "--cwd", self.temp_dir, # Run the command with temp_dir as CWD
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

          print("--- Standard Output ---")
          print(result.stdout)
          print("-----------------------")
          if result.stderr:
              print("--- Standard Error (Warnings/Info) ---")
              print(result.stderr)
              print("------------------------------------")

      except Exception as error:
          logger.info(f"\nAn unexpected error occurred: {error}", file=sys.stderr)
          sys.exit()

      finally:
          # Always attempt to remove the temporary directory
          if os.path.exists(self.temp_dir):
              try:
                  shutil.rmtree(self.temp_dir)
              except OSError as error:
                  logger.info(f"Error removing temporary directory {self.temp_dir}: {error}")

      logger.info("\nProteinMPNN executed successfully.\n")


class Structure_Prediction():

  def __init__(self, input_folder = "", output_folder = "", temp_dir = "", executable="", design_num=0):
    self.inputs = input_folder
    self.outputs = output_folder
    self.temp_dir = temp_dir
    self.executable = executable
    self.design_num = design_num
    self.results = {}
    self.results_dict()

  def results_dict(self):
      for index in range(self.design_num):
          self.results[f"design_{index}"] = {}

  def af2_initial_guess(self):
      """
      Runs AF2 'initial guess' using the PDB files (input) corresponding to the sequences designed by ProteinMPNN.
      """
      logger.info("\nStarting Structure Prediction with AF2 'initial guess' using default configuration...\n")
      
      # Remove the checkpoint file (from a previous run)
      checkpoint_file = "check.point" # File to remove before starting
      try:
          os.remove(checkpoint_file)
      except Exception as e:
          logger.info(f"'{checkpoint_file}': {e}")

      # Remove previous temporary directory (if it exists)
      if os.path.exists(self.temp_dir):
          shutil.rmtree(self.temp_dir)

      # Copy PBD files to a temporary directory (I need to do this because OneDrive gives errors)
      try:
        if not os.path.isdir(self.inputs):
                raise FileNotFoundError()
      except Exception as error:
         logger.info(f"Source directory with PDB files from ProteinMPNN not found: {self.inputs}. {error}")
         sys.exit()
      
      shutil.copytree(self.inputs, self.temp_dir)

      # Construct the command arguments as a list
      command_to_run_in_conda = [
          "python", # Use the python from the 'binder_design' env
          self.executable, #Path to the AF2 executing .py file
          "-pdbdir", self.temp_dir,
          "-outpdbdir", self.outputs
      ]

      # Full command for subprocess
      conda_run_command = ["conda",
          "run",
          "-n", "binder_design", #Name of the conda environment
          "--cwd", self.temp_dir, # Run the command with temp_dir as CWD
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

          for line in result.stdout.splitlines():                
          # Remove leading/trailing whitespace
            stripped_line = line.strip()

            pattern = r"design_\d+" #Looks for a pattern like 'design_1'
            match = re.search(pattern, stripped_line)
            if match:
              # If a match is found, extract the matched string
              design_number = match.group(0)

            # Check if the line looks like a dictionary string
            if stripped_line.startswith('{') and stripped_line.endswith('}'):
              dict_string = stripped_line
              parsed_dict = eval(dict_string)
              self.results[f"{design_number}"]["plddt_binder"] = parsed_dict["plddt_binder"]
              self.results[f"{design_number}"]["pae_interaction"] = parsed_dict["pae_interaction"]
              self.results[f"{design_number}"]["binder_aligned_rmsd"] = parsed_dict["binder_aligned_rmsd"]

          # print("--- Standard Output ---")
          # print(result.stdout)
          # print("-----------------------")
          # if result.stderr:
          #     print("--- Standard Error (Warnings/Info) ---")
          #     print(result.stderr)
          #     print("------------------------------------")

      except Exception as error:
          logger.info(f"\nAn unexpected error occurred: {error}", file=sys.stderr)
          sys.exit()

      finally:
          # Always attempt to remove the temporary directory
          if os.path.exists(self.temp_dir):
              try:
                  shutil.rmtree(self.temp_dir)
              except OSError as error:
                  logger.info(f"Error removing temporary directory {self.temp_dir}: {error}")


      logger.info("\nAF2 'initial guess' executed successfully.\n")

      logger.info(f'\nAF2 results:\n{self.results}\n')

  def ddg(self, index):
      '''
      Uses RosettaScripts to measure binding energy (ddg).
      '''

      # --- Configuration ---
      rosetta_executable = "/home/gsaenzdepip/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease"
      input_pdb = f"AF2_initial_guess_output/5O45_design_{index}_dldesign_0_cycle1_af2pred.pdb"
      input_xml = "ddg.xml"
      output_scorefile = "score.sc" # Default score file name

      # --- Determine expected output PDB name ---
      # Assumes Rosetta appends _0001 before the extension
      pdb_basename = os.path.basename(input_pdb)
      pdb_stem, pdb_ext = os.path.splitext(pdb_basename)
      expected_output_pdb = f"{pdb_stem}_0001{pdb_ext}"

      # --- Build the command ---
      command = [
          rosetta_executable,
          "-s", input_pdb,
          "-parser:protocol", input_xml,
          # Optional: Explicitly name the output score file
          # "-out:file:scorefile", output_scorefile
          # Optional: Prevent PDB output if you ONLY want the score
          # "-out:nooutput" # If you add this, the PDB removal below is unnecessary
      ]

      # --- Variables to store results ---
      ddg_value = None
      rosetta_success = False

      # --- Main Execution Block ---
      try:
        logger.info(f"\nRunning Rosetta command: {' '.join(command)}\n")
        subprocess.run(command, capture_output=True, text=True, check=True, timeout=600)
        logger.info(f"\nRosetta command completed successfully.\n")
        rosetta_success = True

      except Exception as e:
        print(f"An unexpected error occurred while running Rosetta: {e}")

      # --- Parse the Score File (only if Rosetta ran successfully) ---
      if rosetta_success:
          try:
              if not os.path.exists(output_scorefile):
                  raise FileNotFoundError(f"Output file '{output_scorefile}' not found. Did Rosetta create it?")

              with open(output_scorefile, 'r') as f:
                  lines = f.readlines()

              header_line_content = None
              header_line_index = -1
              last_data_line_content = None
              ddg_index = -1

              # Find header line index and content
              for i, line in enumerate(lines):
                  stripped_line = line.strip()
                  if stripped_line.startswith("SCORE:"):
                      header_line_content = stripped_line
                      header_line_index = i
                      break # Found the header

              if header_line_content is None:
                  print(f"Error: Could not find header line starting with 'SCORE:' in '{output_scorefile}'")
              else:
                  # Find the *last* non-empty line after the header that looks like data
                  for i in range(len(lines) - 1, header_line_index, -1):
                      stripped_line = lines[i].strip()
                      # Add checks to skip known non-data lines if necessary
                      if stripped_line and not stripped_line.startswith("SEQUENCE:") and not stripped_line.startswith("REMARK:"):
                          last_data_line_content = stripped_line
                          break # Found the last potential data line

                  if last_data_line_content is None:
                      print(f"Error: Could not find any potential data lines after the header in '{output_scorefile}'")
                  else:
                      headers = header_line_content.split()
                      try:
                          ddg_index = headers.index("ddg")
                      except ValueError:
                          print(f"Error: Column 'ddg' not found in header line of '{output_scorefile}'")
                          print(f"Available columns: {headers}")

                      if ddg_index != -1:
                          values = last_data_line_content.split()
                          if len(values) > ddg_index:
                              try:
                                  ddg_value = float(values[ddg_index])
                                  print(f"Successfully extracted ddg value: {ddg_value}")
                              except (ValueError, IndexError):
                                  print(f"Error: Could not convert value in column 'ddg' to a number on the data line.")
                                  print(f"Data line content: {last_data_line_content}")

          except Exception as e:
              print(f"An unexpected error occurred during score file parsing: {e}")

          finally:
              # Remove score file
              try:
                  if os.path.exists(output_scorefile):
                      os.remove(output_scorefile)
              except Exception as e:
                  print(f"Error removing file '{output_scorefile}': {e}")

              # Remove generated PDB file
              try:
                  # Check if the command included -out:nooutput before attempting removal
                  if "-out:nooutput" not in command:
                      if os.path.exists(expected_output_pdb):
                          os.remove(expected_output_pdb)
                          print(f"Removed '{expected_output_pdb}'")
                      else:
                          print(f"File '{expected_output_pdb}' not found, skipping removal.")
                  else:
                      print(f"PDB output was suppressed (-out:nooutput), skipping removal of '{expected_output_pdb}'.")
              except OSError as e:
                  print(f"Error removing file '{expected_output_pdb}': {e}")

          # --- Final Result ---
          if ddg_value is not None:
              self.results[f"Design_{index}"]["ddg"] = ddg_value

  def filter_designs(self):
    '''
    Filters designs based on pae_interaction, pLDDT, RMSD and, optionally, ddg.
    '''
    logger.info(f'\nStarting Filtering of Designs...\n\n')
    num_designs = 0
    design_id = []
    for design in self.results.keys():
        try:
          if self.results[design]['plddt_binder'] > 80 and self.results[design]['pae_interaction'] < 10 and self.results[design]['binder_aligned_rmsd'] < 1:
              if 'ddg' in self.results[design].keys() and self.results[design]['ddg'] < -40:
                  num_designs += 1
                  design_id.append(design)
              elif 'ddg' not in self.results[design].keys():
                  num_designs += 1
                  design_id.append(design) 

        except Exception:
            continue
        
    logger.info(f'\nNumber of successful designs: {num_designs}\n'
                f'\nRate of successful designs: {num_designs}/{self.design_num}\n'
                f'\nSuccessful designs IDs: {design_id}\n'
                )                 


if __name__ == "__main__":

  RFdif_settings = {
  "input_folder": "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Backbone_Gen/inputs", # Add directory where the PDB file will be downloaded
  "output_folder": "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Backbone_Gen/outputs", # Add directory where the designed proteins will be saved
  "models_folder": "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Backbone_Gen/models", # Add directory where the RFdiffusion models are saved
  "design_num": 100, # Number of backbone designs to generate
  "hotspots": "A56, A115, A123", #Hotspot residues that condition the diffusion process
  "noise": 0, #Noise to add to the inference (the lower the better)
  "contigs": "A17-139/0 75-100", # Configuration of the binder design
  "pdb_id": "5O45", # The PDB ID of the target protein
  "target_name": "PD-L1"
  }

  ProteinMPNN_settings = {
  "input_folder": "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Backbone_Gen/outputs",
  "temp_dir": "/mnt/c/pdb_temp_work",
  "executable": "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/dl_binder_design/mpnn_fr/dl_interface_design.py",
  "output_folder": "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Sequence_Gen_output",
  "relax_cycles": 1, #The number of relax cycles to perform on each structure (default: 1)
  "seqs_per_struct": 1, #The number of sequences to generate for each structure (default: 1)
  "temperature": 0.0001 #The sampling temperature to use when running ProteinMPNN (default: 0.0001)
  }
  AF2_settings = {
  "input_folder": "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Sequence_Gen_output",
  "temp_dir": "/mnt/c/pdb_temp_work",
  "executable": "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/dl_binder_design/af2_initial_guess/predict.py",
  "output_folder": "/mnt/c/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/AF2_initial_guess_output",
  }

############ NOTES #################
# This code should be run in WSL (or Linux) environment. In VS Code go to the >< sign (left bottom corner) and Connect to WSL.
# Initiate 'binder_design' environment in Conda, where ProteinMPNN and AF2 'initial guess' are installed. Choose the python interpreter in that same environment (3.10).
# Initiate 'pymol' environment in Conda when 'pymol' module is needed.


####################################

#1. Backbone design
  backbone = Backbone_Gen(input_folder = RFdif_settings["input_folder"], output_folder = RFdif_settings["output_folder"], models_folder = RFdif_settings["models_folder"], pdb_id = RFdif_settings["pdb_id"], target_id = RFdif_settings["target_name"])
#   backbone.get_pdb_file()
#   backbone.remove_non_protein() #Pymol required
#   backbone.chain_residue_ranges() #Pymol required
  backbone.backbone_gen(design_num = RFdif_settings["design_num"], hotspots = RFdif_settings["hotspots"], noise = RFdif_settings["noise"], contigs = RFdif_settings["contigs"])

#2. Sequence design
  sequence = Sequence_Gen(input_folder = ProteinMPNN_settings["input_folder"], output_folder = ProteinMPNN_settings["output_folder"], temp_dir = ProteinMPNN_settings["temp_dir"], executable = ProteinMPNN_settings["executable"])
  sequence.sequence_design(relax_cycles=ProteinMPNN_settings["relax_cycles"], seqs_per_struct=ProteinMPNN_settings["seqs_per_struct"], temperature=ProteinMPNN_settings["temperature"])

#3. Structure prediction
  folded_protein = Structure_Prediction(input_folder = AF2_settings["input_folder"], output_folder = AF2_settings["output_folder"], temp_dir = AF2_settings["temp_dir"], executable = AF2_settings["executable"], design_num = RFdif_settings["design_num"])
  folded_protein.af2_initial_guess()

#4. OPTIONAL: calculate binding energy (ddg)
  # for index in range(RFdif_settings["design_num"]):
  #   folded_protein.ddg(index)
  #   print(folded_protein.results)

#5. Filter structures based on pae_interaction, pLDDT, RMSD and, optionally, ddg
  folded_protein.filter_designs()
  







