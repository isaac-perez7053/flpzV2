import os
import shutil
import subprocess
from pathlib import Path

def process_mode(mode_name):
    """
    Converts and executes the logic of the given tcsh script in Python.

    Args:
        mode_name (str): The name of the mode.
    """
    # Paths and input preparation
    header_file = f"headerFile_{mode_name}.dat"
    smodes_dir = Path(f"SMODES_{mode_name}")
    boilerplate_dir = Path("boilerplate")
    smodes_dir.mkdir(exist_ok=True)

    # Move header file to SMODES directory
    shutil.move(header_file, smodes_dir / header_file)

    # Extract number of jobs (NumSAM)
    with open(smodes_dir / header_file, 'r') as f:
        lines = f.readlines()

    num_jobs_ind = None
    for line in lines:
        if "NumSAM" in line:
            num_jobs_ind = int(line.split()[1])
            break

    if num_jobs_ind is None:
        raise ValueError("NumSAM not found in the header file.")

    # Iterate and create directories for jobs
    joblist_file = Path("joblist")
    with joblist_file.open("w") as joblist:
        for ii in range(num_jobs_ind + 1):
            dist_dir = smodes_dir / f"dist_{ii}"
            dist_file = f"dist_{mode_name}_{ii}"

            # Copy boilerplate directory
            shutil.copytree(boilerplate_dir, dist_dir)

            # Insert the contents of `dist_file` into `template.abi` and update it
            template_file = dist_dir / "template.abi"
            with open(template_file, "r") as template:
                content = template.readlines()

            # Read the dist file and insert its content
            with open(dist_file, "r") as dist:
                dist_content = dist.read()

            updated_content = []
            for line in content:
                if "CELLDEF" in line:
                    updated_content.append(dist_content)
                else:
                    updated_content.append(line)

            # Remove "CELLDEF" placeholder
            updated_content = [line.replace("CELLDEF", "") for line in updated_content]

            # Write updated content to template.abi
            with open(template_file, "w") as template:
                template.writelines(updated_content)

            # Rename `template.abi` to `dist_II.abi`
            dist_abi_file = dist_dir / f"dist_{ii}.abi"
            template_file.rename(dist_abi_file)

            # Update `jobscript.sh` with the appropriate DISTNAME
            jobscript_file = dist_dir / "jobscript.sh"
            with open(jobscript_file, "r") as jobscript:
                jobscript_content = jobscript.read()

            updated_jobscript_content = jobscript_content.replace("DISTNAME", f"dist_{ii}")

            with open(jobscript_file, "w") as jobscript:
                jobscript.write(updated_jobscript_content)

            # Add job submission command to joblist
            joblist.write(f"cd {dist_dir}; sbatch jobscript.sh; cd -\n")

            # Clean up the dist file
            os.remove(dist_file)

    print(f"Boilerplate setup completed for mode: {mode_name}")

