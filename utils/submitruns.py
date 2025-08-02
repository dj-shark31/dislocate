from itertools import permutations
import subprocess
import os

def find_compositions(elements, delta_x, remove_binaries="false", remove_unaries="false"):
    """Find all possible compositions for the given elements with increment delta_x"""
    nelems = len(elements)
    bins = int(1.0 / delta_x)
    
    # Generate all possible compositions that sum to bins
    comps = []
    for i in range(bins + 1):
        for j in range(bins + 1 - i):
            if nelems == 2 and i + j == bins:
                comps.append([i, j])
            elif nelems == 3:
                k = bins - i - j
                if k >= 0:
                    comps.append([i, j, k])
                    
    # Convert to fractions
    fracs = [[x * delta_x for x in comp] for comp in comps]
    
    # Generate all permutations
    comps_to_run = set()
    for fr in fracs:
        for p in permutations(fr):
            comps_to_run.add(p)

    # Filter out unaries and binaries if requested
    filtered_comps = set()
    for comp in comps_to_run:
        # Skip compositions that are pure elements (unaries)
        if remove_unaries == "true" and 1.0 in comp:
            continue
            
        # For ternary systems, skip binary compositions (where one element is 0)
        # but keep unary compositions
        if remove_binaries == "true" and nelems == 3 and 0.0 in comp and not any(x == 1.0 for x in comp):
            continue
            
        filtered_comps.add(comp)
            
    return filtered_comps

def submit_savio_job(script_name, script_args_dict, job_name="job", partition="savio4_htc", 
                     cpus=10, time="5:00:00", qos="savio_lowprio", account="co_chrzangroup", nodes=1, ntasks_per_node=1):
    """
    Submit a generic job to Savio cluster
    
    Parameters:
        script_name (str): Name of the Python script to run
        script_args_dict (dict): Dictionary of arguments to pass to the script
        job_name (str): Name for the SLURM job
        partition (str): Savio partition to use
        cpus (int): Number of CPUs per task
        time (str): Wall time limit in format "HH:MM:SS"
        qos (str): Quality of service level
        account (str): Account to charge
    """
    
    # Convert dictionary arguments to command line format
    script_args = " ".join([f"--{k} {v}" for k,v in script_args_dict.items()])
    
    # Write batch script
    with open("runjob.sh", "w") as runjob:
        runjob.write("#!/bin/bash\n")
        runjob.write(f"#SBATCH --job-name={job_name}\n")
        runjob.write(f"#SBATCH --account={account}\n")
        runjob.write(f"#SBATCH --nodes={nodes}\n")
        runjob.write(f"#SBATCH --ntasks-per-node={ntasks_per_node}\n")
        runjob.write(f"#SBATCH --partition={partition}\n")
        runjob.write(f"#SBATCH --cpus-per-task={cpus}\n")
        runjob.write(f"#SBATCH --qos={qos}\n")
        runjob.write(f"#SBATCH --time={time}\n")
        runjob.write("#SBATCH --requeue\n\n")
        
        # Add command to run the script
        runjob.write(f"python {script_name} {script_args}\n")
        
    # Submit the job
    subprocess.run(["sbatch", "runjob.sh"])

def composition_dir(elements, composition):
    """Generate a directory name for a given composition"""
    output_dir = ''.join(elements) + "/"
    os.makedirs(output_dir, exist_ok=True)
    for i, elem in enumerate(elements):
        if len(elements) == 2:
            output_dir += f"{elem}{composition[i]:.2f}"
        elif len(elements) == 3:
            output_dir += f"{elem}{composition[i]:.3f}"
    os.makedirs(output_dir, exist_ok=True)
    return output_dir