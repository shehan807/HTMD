import sys, os
from jobflow import job, Flow, run_locally
import argparse
import matplotlib.pyplot as plt

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, current_dir)
sys.path.insert(0, parent_dir)
from analyze import get_rho, get_rdf
from inputs import input_df, create_OSG_csv
from visualize import _figure_settings
def get_cation_anion(system):
    """
    Extracts cation and anion from the system name.
    Example: "p4444_dmbs_water" --> ("p4444", "dmbs")
    """
    parts = system.split("_")
    return parts[0], parts[1]

def main():
    # Shades of red for P4444 systems
    p4444_colors = {
        "bzso3": "#e41a1c",  # red
        "tso":   "#f16975",  # light red
        "dmbs":  "#99000d",  # dark red
        "tmbs":  "#fb6a4a",  # salmon
    }
    
    # Shades of orange for P4448 systems
    p4448_colors = {
        "bzso3": "#ff7f00",  # orange
        "tso":   "#fdae6b",  # light orange
        "dmbs":  "#e6550d",  # burnt orange
        "tmbs":  "#fd8d3c",  # apricot
    }
    
    # Shades of green for P5555 systems
    p5555_colors = {
        "bzso3": "#4daf4a",  # green
        "tso":   "#a1d99b",  # light green
        "dmbs":  "#006d2c",  # dark green
        "tmbs":  "#74c476",  # mint green
    }
 
    base_dir =  "/storage/home/hhive1/sparmar32/projects/HTMD/force_field/OPLS/example_lcst"
    systems = [
                "p4444_bzso3_water",
                "p4444_tso_water",
                "p4444_dmbs_water",
                "p4444_tmbs_water",
                "p4448_bzso3_water",
                "p4448_tso_water",
                "p4448_dmbs_water",
                "p4448_tmbs_water",
                "p5555_bzso3_water",
                "p5555_tso_water",
                "p5555_dmbs_water",
                "p5555_tmbs_water",
    ]
    system_labels = {
        "p4444_bzso3_water": r"$\rm P^+_{P_{4444}}\cdots O^-_{BzSO_3}$",
        "p4444_tso_water":   r"$\rm P^+_{P_{4444}}\cdots O^-_{TSO}$",
        "p4444_dmbs_water":  r"$\rm P^+_{P_{4444}}\cdots O^-_{DMBS}$",
        "p4444_tmbs_water":  r"$\rm P^+_{P_{4444}}\cdots O^-_{TMBS}$",
        
        "p4448_bzso3_water": r"$\rm P^+_{P_{4448}}\cdots O^-_{BzSO_3}$",
        "p4448_tso_water":   r"$\rm P^+_{P_{4448}}\cdots O^-_{TSO}$",
        "p4448_dmbs_water":  r"$\rm P^+_{P_{4448}}\cdots O^-_{DMBS}$",
        "p4448_tmbs_water":  r"$\rm P^+_{P_{4448}}\cdots O^-_{TMBS}$",
        
        "p5555_bzso3_water": r"$\rm P^+_{P_{5555}}\cdots O^-_{BzSO_3}$",
        "p5555_tso_water":   r"$\rm P^+_{P_{5555}}\cdots O^-_{TSO}$",
        "p5555_dmbs_water":  r"$\rm P^+_{P_{5555}}\cdots O^-_{DMBS}$",
        "p5555_tmbs_water":  r"$\rm P^+_{P_{5555}}\cdots O^-_{TMBS}$",
    }
    # LCST systems to highlight
    lcst_systems = {
        "p4444_dmbs_water",
        "p4444_tmbs_water",
        "p4444_tso_water",
    }
    
    max_gr = 0  # Keep track of the max g(r) value across all systems

    CONC = "0.15"
    TEMP = "300"
    fig, ax = _figure_settings() 
    for system in systems:
        print(f"RDF for {system} @ {float(CONC)*100}%% H2O, {TEMP} K")
        input_job = input_df(
                    system=system, 
                    conc=CONC, 
                    temp=TEMP, 
                    dir=base_dir, 
                    pdbfile="npt_final-4.pdb", 
                    dcdfile="md_npt-4.dcd", 
                    resfiles="",
                    ff_files="",
        )

        rho_job = get_rho(input_job)
        selection = {
            'O-H': [
                '(name O* and not resname HOH)',  # O atoms not in HOH or TBP
                '(name P*)'                            # P atoms in TBP
            ],
        }
        
        r, gr = get_rdf(input_job, selections=selection)
        max_gr = max(max_gr, max(gr))
        
        cation, anion = get_cation_anion(system)
        if cation == "p4444":
            color = p4444_colors.get(anion, "#d73027")
        elif cation == "p4448":
            color = p4448_colors.get(anion, "#f46d43")
        elif cation == "p5555":
            color = p5555_colors.get(anion, "#1a9850")
        else:
            color = "gray"  # fallback
        linestyle = "-"
        
        if system in lcst_systems:
            ax.plot(r, gr, label=system_labels[system], color=color,
                    linestyle='-', linewidth=8.0, alpha=0.9)
        else:
            ax.plot(r, gr, label=system_labels[system], color=color,
                    linestyle='-', linewidth=3.0, alpha=0.3)
         


        del input_job
    ax.set_xlim(2, 8.0)
    ax.set_ylim(0, max_gr * 1.1)
    ax.axhline(1.0, color='gray', linestyle='--', linewidth=1)
    label_fontsize = 26
    tick_fontsize = 24
    legend_fontsize = 22
    
    ax.set_xlabel("r (Å)", fontsize=label_fontsize)
    ax.set_ylabel("g(r)", fontsize=label_fontsize)
    ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)
    ax.legend(frameon=False, fontsize=legend_fontsize)
    
    plt.tight_layout()
    plt.savefig("rdfs.png", dpi=300, bbox_inches='tight')
        
        

if __name__ == "__main__":
    main()
