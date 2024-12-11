# Steps towards generating SAPT FF Paramaters

## 1. Compute GRAC shift of both monomers  

In Psi4, SAPT(DFT) dimer calculations require a user-defined GRAC shift, 

```bash
    sapt_dft_grac_shift_a GRAC_SHIFT_A
    sapt_dft_grac_shift_b GRAC_SHIFT_B
```

which are computed for each \textit{monomer}. The GRAC shift is defined as:

$$ GRAC = IP + HOMO $$

where HOMO is the Kohn-Sham HOMO energy and IP is the ionization potential,

$$ IP = E_{2} - E_{1} $$

for which the energy difference can either be determined via NIST Chemistry Webbook or by the difference between the molecule as given, and the molecule after one electron has bee removed (e.g., the energy difference between a neutral molecule and its cation).

## 2. Create dimer geometries

There are two ways of creating ~1000s of dimer configurations:

### a. Random generation 

### b. Sampled from low-fidelity MD simulation

## 3. Compute dimer energies in Psi4

### a. Fit charges 

### b. Create Psi4 input files 

## 4. Compare with existing SAPT FF terms (optional, but recommended)

## 5. Obtain polarizabilities via CAMCasp

## 6. Fit to SAPT functional form 

