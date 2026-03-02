
**must be replaced** with the full content of a geometry file.

You may either:

- Use the original `molecular_geometry.in`, or  
- Use the automatically generated `geometry.in` produced by `rototraslation.py`.

This substitution is required for the calculations to run correctly.

------------------------------------------------------------------------

## Modifying the Molecular Configuration

To obtain energy profiles corresponding to different intermolecular arrangements (e.g., parallel-displaced, T-shaped, rotated configurations, or distance scans), it is possible to:

- Apply rigid roto-translations to the two benzene molecules using `rototraslation.py`  
- Generate a new `geometry.in` file  
- Replace the `<molecular_geometry>` block in the input files accordingly  

No further modifications to the input files are required, provided the molecular geometry block is correctly replaced.

------------------------------------------------------------------------

## Scientific Reference

The input files were **kindly provided by the authors** of the following article:

*Substituent and Heteroatom Effects on π--π Interactions: Evidence That Parallel-Displaced π-Stacking is Not Driven by Quadrupolar Electrostatics*,  
Journal of the American Chemical Society.  
https://pubs.acs.org/doi/10.1021/jacs.4c13291

If you use these inputs in your research, please cite the original article accordingly.

------------------------------------------------------------------------

## Software Requirements

- Q-Chem (version 6.3 or compatible)  
- Python 3.x (for running `rototraslation.py`)  
- Sufficient computational resources for SAPT/XSAPT calculations  
