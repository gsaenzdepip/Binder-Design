#This code should be run in WSL (or Linux) environment.

import docker
import requests
import sys
import subprocess
from pymol import cmd
import os
import sys
import shutil 
import __main__

class Backbone_Gen():
  def __init__(self, input_folder = "", output_folder = "", models_folder = "", model_num = 100, hotspots = "", noise = 0.5, contigs = "", pdb_id = ""):
    self.inputs = input_folder
    self.outputs = output_folder
    self.models = models_folder
    self.model_num = model_num
    self.hotspots = hotspots
    self.noise = noise
    self.contigs = contigs
    self.pdb_id = pdb_id

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

  def backbone_gen(self):
    '''
    Task:
      Designs backbones of binders against a target protein using RFdiffusion.
    '''
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
        f"inference.num_designs={self.model_num}",
        f"contigmap.contigs=[{self.contigs}]",
        f"ppi.hotspot_res=[{self.hotspots}]",
        f"denoiser.noise_scale_ca={self.noise}",
        f"denoiser.noise_scale_frame={self.noise}"
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

class Sequence_Gen():

  def __init__(self, input_folder = "", output_folder = "", pdb_id = ""):
    self.inputs = input_folder
    self.outputs = output_folder
    self.pdb_id = pdb_id

  def sequence_design():
      """Runs ProteinMPNN using the PDB files (input) corresponding to the backbones designed by RFdiffusion."""

      # Remove the checkpoint file (from a previous run)
      checkpoint_file = "check.point" # File to remove before starting
      try:
          os.remove(checkpoint_file)
      except OSError as e:
          print(f"'{checkpoint_file}': {e}", file=sys.stderr)

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
          executable_path, #Path to the ProteinMPNN executing .py file
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

          print("\nProteinMPNN executed successfully.\n")
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

  def fasta_gen(self):
    '''
    Task:
      - Creates two FASTA files from the ProteinMPNN output.
    '''
    # Open the input file and read all lines
    with open(f"{self.outputs}/seqs/{self.pdb_id}_model_0.fa", "r") as infile:
        lines = infile.readlines()

    # Note: Python lists are zero-indexed. So line 4 is index 3 and line 6 is index 5.
    sample1_seq = lines[3].strip()  # Sequence for sample 1 (line 4)
    # sample2_seq = lines[5].strip()  # Sequence for sample 2 (line 6)

    # Write sample 1 to its own FASTA file with a custom header
    with open(f"{self.outputs}/seqs/{self.pdb_id}_model_0_seq_1.fasta", "w") as out1:
        out1.write(f">{self.pdb_id}_seq_1\n")
        out1.write(sample1_seq + "\n")

    # Write sample 2 to its own FASTA file with a custom header
    # with open(f"{self.outputs}/seqs/{self.pdb_id}_model_0_seq_2.fasta", "w") as out2:
    #     out2.write(f">{self.pdb_id}_seq_2\n")
    #     out2.write(sample2_seq + "\n")


class Folding():

  def __init__(self):
     pass


  def colab_fold(self):
    # Initialize the Docker client from environment variables
    client = docker.from_env()

    # Define your paths
    weights_path = "C:/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Sequence_Evaluation/cache/colabfold"
    input_path = "C:/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Sequence_Gen/output/seqs"
    output_path = "C:/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Sequence_Evaluation/output"

    # Create a dictionary for volume bindings
    volumes = {
        weights_path: {'bind': '/colabfold', 'mode': 'rw'},
        input_path: {'bind': '/seqs', 'mode': 'rw'},
        output_path: {'bind': '/output', 'mode': 'rw'}
    }

    # Commands to be executed inside the container:

    command = "colabfold_batch --model-type alphafold2_ptm --num-recycle 0 --num-models 1 --data /colabfold /seqs/1WHQ_model_0_seq_1.fasta /output"
    # command = "colabfold_batch --msa-only --data /colabfold /seqs/1WHQ_model_0_seq_1.fasta /output" # Get only the MSA

    # command = "colabfold_batch --help"

    # For GPU support, we can specify device_requests (Docker SDK version 4.4+)
    device_requests = [docker.types.DeviceRequest(count=-1, capabilities=[['gpu']])]

    try:
        # Run the container
        container = client.containers.run(
            image="ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2",
            command=command,
            volumes=volumes,
            runtime="nvidia",  # Make sure NVIDIA Container Toolkit is installed
            device_requests=device_requests,
            detach=True,
            # tty=False,
            # stdin_open=False,
            # environment={"HYDRA_FULL_ERROR": "1"} # Allows to print errors
        )

        # Attach to logs
        log_stream = container.logs(stream=True)
        for log in log_stream:
            sys.stdout.write(log.decode("utf-8"))  # Writes log output without adding extra newlines
            sys.stdout.flush()  # Ensures immediate output

        # Wait for the container to finish
        container.wait()


        
    except docker.errors.DockerException as e:
        print(f"An error occurred: {e}")
    
    container.remove()

if __name__ == "__main__":

  RFdif_settings = {
  "input_folder": "C:/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Backbone_Gen/inputs", # Add directory where the PDB file will be downloaded
  "output_folder": "C:/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Backbone_Gen/outputs", # Add directory where the designed proteins will be saved
  "models_folder": "C:/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Backbone_Gen/models", # Add directory where the RFdiffusion models are saved
  "model_num": 20, # Number of designs to generate
  "hotspots": "A56, A115, A123",
  "noise": 0,
  "contigs": "A17-145/0 50-100", # Configuration of the protein designC
  "pdb_id": "5O45" # Add the PDB ID
  }

  ProteinMPNN_settings = {
  "input_folder": "C:/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Backbone_Gen/outputs",
  "output_folder": "C:/Users/gsaen/OneDrive - UW/Python/Structural_Bioinformatics_Course/Protein Design/Python/Sequence_Gen/output"
  }

# 1. Backbone design
  backbone = Backbone_Gen(input_folder = RFdif_settings["input_folder"], output_folder = RFdif_settings["output_folder"], models_folder = RFdif_settings["models_folder"], model_num = RFdif_settings["model_num"], hotspots = RFdif_settings["hotspots"], noise = RFdif_settings["noise"], contigs = RFdif_settings["contigs"], pdb_id = RFdif_settings["pdb_id"])
  # backbone.get_pdb_file()
  # backbone.remove_non_protein()
  # backbone.chain_residue_ranges()
  backbone.backbone_gen()

  sequence = Sequence_Gen(input_folder = ProteinMPNN_settings["input_folder"], output_folder = ProteinMPNN_settings["output_folder"], pdb_id=RFdif_settings["pdb_id"])
  # sequence.protein_mpnn()
  # sequence.fasta_gen()

  folded_protein = Folding()
  # folded_protein.colab_fold()







