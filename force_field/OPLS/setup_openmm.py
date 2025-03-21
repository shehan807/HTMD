"""
Module Name: 
Authors: Shehan M. Parmar
Date: 2024-11-25
Description:
Usage:
"""
import os
import argparse
import glob
import pandas as pd
import yaml
import subprocess
import shutil 
import mdtraj 
from tabulate import tabulate
from pathlib import Path 

from openff.toolkit import ForceField, Molecule, Topology
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange.components._packmol import UNIT_CUBE, pack_box
from openff.interchange.drivers import get_summary_data
from openff.toolkit.utils import get_data_file_path
from openmoltools import packmol

def parser():
    parser = argparse.ArgumentParser(description="Workflow script for molecular dynamics simulations")

    # Input arguments for configuration files
    parser.add_argument(
        "-m",
        "--molecule_file",
        type=str,
        required=True,
        help="Path to the CSV file containing molecule properties (e.g., molecule_config.csv)"
    )
    parser.add_argument(
        "-c",
        "--conditions_file",
        type=str,
        required=True,
        help="Path to the YAML file containing experimental conditions (e.g., config.yaml)"
    )
    parser.add_argument(
        "-s",
        "--sif_file",
        type=str,
        default="/storage/home/hhive1/sparmar32/projects/HTMD/force_field/OPLS/ligpargen-image/LPG.sif",
        help="Path to the .sif file for LigParGen (default: /storage/home/hhive1/sparmar32/projects/HTMD/force_field/OPLS/ligpargen-image/LPG.sif)"
    )
    parser.add_argument(
        "-w",
        "--water_model",
        type=str,
        default="tip3p",
        help="Name of water model (must be lower case): spce, swm4-dp, tip3p, TIP3P-FB, tip4p, tip4p-ew, tip4p-fb, tip5p"
    )

    args = parser.parse_args()

    # Load molecule data from CSV
    molecule_data = pd.read_csv(args.molecule_file)
    molecule_map = molecule_data.set_index("NAME")[["MOL", "CHARGE", "SMILES"]].to_dict("index")

    # Load conditions data from YAML
    with open(args.conditions_file, "r") as f:
        conditions = yaml.safe_load(f)

    sif_file = args.sif_file
    water_model = args.water_model 

    return molecule_map, conditions, sif_file, water_model

def run_ligpargen(molecule, charge, smiles, sif_file, water_model="tip3p",water_models_path="/storage/home/hhive1/sparmar32/projects/HTMD/force_field/OPLS/water_models"):
    mol_dir = f"{molecule}"
    library_dir = "/storage/home/hhive1/sparmar32/projects/HTMD/force_field/OPLS/il_xml_library/IL/2009IL/XML"
    library_xml_file = os.path.join(library_dir, f"{molecule}.xml")
    library_pdb_file = os.path.join(library_dir, f"{molecule}.pdb")
    if molecule != "HOH":
        if os.path.isfile(library_xml_file) and os.path.isfile(library_pdb_file):
            print(f"Found {molecule} in library; using existing XML and PDB.")
            xml_file = shutil.copy(library_xml_file, f"{molecule}.xml")
            pdb_file = shutil.copy(library_pdb_file, f"{molecule}.pdb") 
        else:
            cmd = f"apptainer exec --bind $(pwd):/opt/output {sif_file} bash -c 'ligpargen -n {molecule} -p {molecule} -r {molecule} -c {charge} -o 3 -cgen CM1A -s \"{smiles}\"'"
            subprocess.run(cmd, shell=True)
            xml_file = os.path.join(mol_dir, f"{molecule}.openmm.xml")
            pdb_file = os.path.join(mol_dir, f"{molecule}.openmm.pdb")

            xml_file = shutil.move(xml_file, f"{molecule}.xml")
            pdb_file = shutil.move(pdb_file, f"{molecule}.pdb")
            
            shutil.rmtree(mol_dir)
    elif molecule == "HOH":
        xml_file = os.path.join(water_models_path, water_model, f"{water_model}.xml")
        pdb_file = os.path.join(water_models_path, water_model, f"{water_model}.pdb")

        xml_file = shutil.copy(xml_file, f"{molecule}.xml")
        pdb_file = shutil.copy(pdb_file, f"{molecule}.pdb")
    return xml_file, pdb_file

def gather_xml_files(xml_files, output_dir="ffdir", script_path="/storage/home/hhive1/sparmar32/projects/HTMD/force_field/OPLS/XML_REFORMAT"):
    cmd = f"python {os.path.join(script_path,'LPG_reformat.py')} --xml_files {' '.join(xml_files)} --output {output_dir}; ls"
    subprocess.run(cmd, shell=True)
    output = glob.glob(str(Path(output_dir) / "*_sorted.xml"))
    return output

def molar_mass(mol):
    mass = 0.0
    for atom in mol.atoms:
        mass += atom.mass
    return mass

def get_num_molecules(comp1, comp2, conc, map):
    
    molecules1 = [Molecule.from_smiles(map[comp]["SMILES"], allow_undefined_stereo=True) for comp in comp1]
    molecules2 = [Molecule.from_smiles(map[comp]["SMILES"], allow_undefined_stereo=True) for comp in comp2]

    mass1 = [molar_mass(Molecule) for Molecule in molecules1]
    mass2 = [molar_mass(Molecule) for Molecule in molecules2]
    
    # This is not generalizable to 3-component systems or different materials, but works for now
    
    # assign heavier component 200 molecules (conc is for molecule1!)
    if sum(mass1) > sum(mass2):
        if conc == 0.0:
            n1 = 0
            n2 = 200
        else:
            n1 = 200
            n2 = int((1 / sum(mass2)) * (n1*sum(mass1) / conc - n1*sum(mass1)))
    elif sum(mass2) > sum(mass1):
        if (1-conc) == 0.0:
            n2 = 0
            n1 = 200
        else:
            n2 = 200
            n1 = int((1 / sum(mass1)) * ((n2*sum(mass2)) / (1-conc) - n2*sum(mass2)) )
    print("%%%%% MASS AND NUMBER %%%%%")
    print(mass1, mass2, n1, n2)
    
    num_molecules1 = len(molecules1) * [n1]
    num_molecules2 = len(molecules2) * [n2]

    comps         = comp1 + comp2
    molecules     = molecules1 + molecules2
    num_molecules = num_molecules1 + num_molecules2
    
    # avoid packmol error  ERROR: Number of molecules of type 1  set to less than 1
    for i, num_molecule in enumerate(num_molecules):
        if num_molecule == 0.0:
            comps[i] = None
            molecules[i] = None
            num_molecules[i] = None
    comps_final = [comp for comp in comps if comp is not None]
    molecules_final = [molecule for molecule in molecules if molecule is not None]
    num_molecules_final = [num_molecule for num_molecule in num_molecules if num_molecule is not None]
    return comps_final, molecules_final, num_molecules_final


def create_topology(molecules, num_molecules, target_density=800, system_dir = ".", system_file="system.pdb"):
    """
    TODO! Doesn't work for unique set of RES and HETATM names (!!!)
    Create a packed topology file (system.pdb) based on molecule list, 
    number of copies, and target density.
    """
    topology = pack_box(
        molecules=molecules,
        number_of_copies=num_molecules,
        target_density=target_density * unit.kilogram / unit.meter**3,
        box_shape="cubic",
        verbose=True
    )
    system_file_path = os.path.join(system_dir, system_file)
    topology.to_file(system_file_path)
    return system_file_path

def create_job(system_pdb_path, temp, all_xml, job_dir, slurm_job_name, template_dir="/storage/home/hhive1/sparmar32/projects/HTMD/jobs/templates"):
    """
    Creates job-specific SLURM scripts and updates parameters.
    """
    # Paths to the templates
    slurm_template = os.path.join(template_dir, "run.slurm")
    wrapper_template = os.path.join(template_dir, "slurm_wrapper.sh")

    # Paths to the job-specific scripts
    job_slurm = os.path.join(job_dir, "run.slurm")
    job_wrapper = os.path.join(job_dir, "slurm_wrapper.sh")

    # Create the job directory if it doesn't exist
    os.makedirs(job_dir, exist_ok=True)

    # Copy the slurm/wrapper script unchanged
    shutil.copy(slurm_template, job_slurm)
    shutil.copy(wrapper_template, job_wrapper)
    

    # Customize the SLURM/wrapper script
    with open(slurm_template, "r") as f:
        slurm_content = f.read()
    # Replace placeholders in the SLURM script
    slurm_content = slurm_content.replace("###JOB_NAME###", slurm_job_name)
    # Write the updated SLURM script to the job directory
    with open(job_slurm, "w") as f:
        f.write(slurm_content)
    
    with open(wrapper_template, "r") as f:
        wrapper_content = f.read()
    # Replace placeholders in the SLURM script
    wrapper_content = wrapper_content.replace("###TEMP###", str(temp))

    residue_xmls = [xml.replace(".xml", "_residues.xml") for xml in all_xml]
    wrapper_content = wrapper_content.replace("###RES_FILE###", ":".join(residue_xmls))
    wrapper_content = wrapper_content.replace("###PDB_FILE###", system_pdb_path)
    wrapper_content = wrapper_content.replace("###FF_FILE###", ":".join(all_xml))
    # Write the updated SLURM script to the job directory
    with open(job_wrapper, "w") as f:
        f.write(wrapper_content)

    print(f"Created SLURM script: {job_slurm}")
    print(f"Created wrapper script: {job_wrapper}")

    return job_slurm, job_wrapper

def print_summary(jobdir_list, jobdir_status):
    """
    Prints a formatted summary of job directories with path, status, and metadata.
    """
    # Prepare data for the table
    summary_data = []
    for i, path in enumerate(jobdir_list):
        # Extract metadata from the path
        name, conc, temp = path.split(os.sep)[-3:]
        summary_data.append({
            "Job#": i,
            "Name": name,
            "Concentration": conc,
            "Temperature": temp,
            "Status (ns)": jobdir_status[i]
        })

    # Print the summary table
    print(tabulate(summary_data, headers="keys", tablefmt="grid"))


def main():
    """
    """
    # Parse command-line arguments
    molecule_map, conditions, sif_file, water_model = parser()
    jobdir_list = []
    jobdir_status = []
    for mixture in conditions["mixtures"]:
        name = mixture["name"]
        comp1 = mixture["component1"]
        comp2 = mixture["component2"]
        mols = comp1 + comp2
        xmls = []
        pdbs = []
        for mol in mols:
            MOL = molecule_map[mol]["MOL"]
            CHARGE = molecule_map[mol]["CHARGE"]
            SMILES = molecule_map[mol]["SMILES"]

            mol_xml, mol_pdb = run_ligpargen(MOL, CHARGE, SMILES, sif_file, water_model=water_model)
            print(f"generated {mol_xml} and {mol_pdb}")
            xmls.append(mol_xml)
            pdbs.append(mol_pdb)
        print(f"xmls: {xmls}")
        xmls = list(set(xmls))
        all_xml = gather_xml_files(xmls)

        # check if ffdir exists in mixture "name" directory (i.e., save only one copy of ff)
        name_dir  = os.path.join(os.getcwd(), str(name))
        name_ffdir = os.path.join(name_dir, 'ffdir')
        if os.path.exists(name_ffdir):
            diff = subprocess.run(['diff', '-rq', 'ffdir', name_ffdir], capture_output=True, text=True)
            if diff.returncode == 0: 
                print(f"'ffdir' exists in {name_dir} with NO DIFFERENCE.")
            else:
                print(f"WARNING: the existing 'ffdir' is different!!")
                print(diff.stdout)
        else: 
            os.makedirs(name_dir)
            root_ffdir = shutil.copytree('ffdir', name_ffdir)
            print(f"created {root_ffdir}")

        for conc in mixture["concentrations"]:
            for temp in mixture["temperatures"]:
                dir = os.path.join(os.getcwd(), str(name), str(conc), str(temp))
                jobdir_list.append(dir)
                if os.path.exists(dir):
                    try:
                        sim_out = os.path.join(dir, "simulation_output")
                        num_dcds = [file for file in os.listdir(sim_out) if file.endswith(".dcd")]
                        jobdir_status.append(f"{len(num_dcds)*5} ns")
                        print(f"{dir} already exists with {len(num_dcds)*5} ns complete! Skipping to next directory.\n")
                    except FileNotFoundError:
                        jobdir_status.append(f"{0} ns")
                        print(f"{dir} exists without 'simulation_output'. Still skipping to next directory.")
                    continue
                else:
                    jobdir_status.append(f"{0} ns")
                    os.makedirs(dir)
                    print(f"Directory {dir} created")

                # run openff system.pdb generator
                comps, molecules, num_molecules = get_num_molecules(comp1, comp2, conc, molecule_map)
                pdbs = [molecule_map[comp]["MOL"]+".pdb" for comp in comps]
                print(pdbs, comps, molecules, num_molecules)
                system_pdb = packmol.pack_box(pdbs, num_molecules)
                system_pdb_path = os.path.join(dir, 'system.pdb')
                system_pdb.save_pdb(system_pdb_path)

                print(f"Created system.pdb in {system_pdb_path}.")
                # copy system.pdb to dir
                
                slurm_job_name = f"{name}_{conc}_{temp}"
                all_xml = [os.path.join(name_dir, xml) for xml in all_xml] 
                
                create_job(system_pdb_path, temp, all_xml, dir, slurm_job_name)
        _ = [os.remove(x) for x in xmls if os.path.isfile(x)]
        _ = [os.remove(p) for p in pdbs if os.path.isfile(p)]
        shutil.rmtree("ffdir")

    print_summary(jobdir_list, jobdir_status)
if __name__ == "__main__":
    main()
