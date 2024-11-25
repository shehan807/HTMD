"""
Module Name: 
Authors: Shehan M. Parmar
Date: 2024-11-25
Description:
Usage:
"""
import os
import argparse
import pandas as pd
import yaml
import subprocess
import shutil 
import mdtraj 

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
        default="/storage/home/hcoda1/4/sparmar32/p-jmcdaniel43-0/scripts/MD_Workflow/force_field/OPLS/ligpargen-image/LPG.sif",
        help="Path to the .sif file for LigParGen (default: /storage/home/hcoda1/4/sparmar32/p-jmcdaniel43-0/scripts/MD_Workflow/force_field/OPLS/ligpargen-image/LPG.sif)"
    )

    args = parser.parse_args()

    # Load molecule data from CSV
    molecule_data = pd.read_csv(args.molecule_file)
    molecule_map = molecule_data.set_index("NAME")[["MOL", "CHARGE", "SMILES"]].to_dict("index")

    # Load conditions data from YAML
    with open(args.conditions_file, "r") as f:
        conditions = yaml.safe_load(f)

    sif_file = args.sif_file

    return molecule_map, conditions, sif_file

def run_ligpargen(molecule, charge, smiles, sif_file):
    cmd = f"apptainer exec --bind $(pwd):/opt/output {sif_file} bash -c 'ligpargen -n {molecule} -p {molecule} -r {molecule} -c {charge} -o 3 -cgen CM1A -s \"{smiles}\"'"
    subprocess.run(cmd, shell=True)
    mol_dir = f"{molecule}"
    xml_file = os.path.join(mol_dir, f"{molecule}.openmm.xml")
    pdb_file = os.path.join(mol_dir, f"{molecule}.openmm.pdb")

    xml_file = shutil.move(xml_file, f"{molecule}.xml")
    pdb_file = shutil.move(pdb_file, f"{molecule}.pdb")

    shutil.rmtree(mol_dir)
    return xml_file, pdb_file

def merge_xml_files(xml_files, output_dir="ffdir", script_path="/storage/home/hcoda1/4/sparmar32/p-jmcdaniel43-0/scripts/MD_Workflow/force_field/OPLS/XML_REFORMAT"):
    cmd = f"python {os.path.join(script_path,'LPG_reformat.py')} --xml_files {' '.join(xml_files)} --output {output_dir}; ls"
    subprocess.run(cmd, shell=True)
    base_names = [
        os.path.splitext(os.path.basename(file))[0]
        for file in xml_files
    ]
    output = os.path.join(output_dir, "_".join(base_names) + ".xml")
    return output

def molar_mass(mol):
    mass = 0.0
    for atom in mol.atoms:
        mass += atom.mass
    return mass

def get_num_molecules(comp1, comp2, conc, map):
    
    molecules1 = [Molecule.from_smiles(map[comp]["SMILES"]) for comp in comp1]
    molecules2 = [Molecule.from_smiles(map[comp]["SMILES"]) for comp in comp2]

    mass1 = [molar_mass(Molecule) for Molecule in molecules1]
    mass2 = [molar_mass(Molecule) for Molecule in molecules2]
    
    # This is not generalizable to 3-component systems or different materials, but works for now
    
    # assign heavier component 200 molecules (conc is for molecule1!)
    if sum(mass1) > sum(mass2):
        if conc == 0.0:
            n1 = 0
        else:
            n1 = 200
        n2 = int((1 / sum(mass2)) * (n1*sum(mass1) / conc - n1*sum(mass1)))
    elif sum(mass2) > sum(mass1):
        if (1-conc) == 0.0:
            n2 = 0
        else:
            n2 = 200
        print(mass1, mass2, n2)
        n1 = int((1 / sum(mass1)) * ((n2*sum(mass2)) / (1-conc) - n2*sum(mass2)) )
    
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

def create_job():
    """
    Copy necessary job scripts from /template directory to 
    current directory.
    TO-DO: Eventually, this should all be replaced by jobflow libraries. 
    """
    job_name =  

def main():
    """
    """
    # Parse command-line arguments
    molecule_map, conditions, sif_file = parser()
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

            mol_xml, mol_pdb = run_ligpargen(MOL, CHARGE, SMILES, sif_file)
            print(f"generated {mol_xml} and {mol_pdb}")
            xmls.append(mol_xml)
            pdbs.append(mol_pdb)
        print(f"xmls: {xmls}")
        merged_xml = merge_xml_files(xmls)
        print(f"generated {merged_xml}\n")

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
                if os.path.exists(dir):
                    try:
                        sim_out = os.path.join(dir, "simulation_output")
                        num_dcds = [file for file in os.listdir(sim_out) if file.endswith(".dcd")]
                        print(f"{dir} already exists with {len(num_dcds)*5} ns complete! Skipping to next directory.\n")
                    except FileNotFoundError:
                        print(f"{dir} exists without 'simulation_output'. Still skipping to next directory.")
                    continue
                else:
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

        _ = [os.remove(x) for x in xmls if os.path.isfile(x)]
        _ = [os.remove(p) for p in pdbs if os.path.isfile(p)]
        shutil.rmtree("ffdir")

if __name__ == "__main__":
    main()
