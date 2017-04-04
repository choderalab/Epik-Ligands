# Epik-Ligands
This repository contains ligands parametrized for MCCE using Epik to predict protonation states and OpenEye for charge assignment. 

See [this page](summary.md) for a summary of the results.

Each run has the following files


* Output files containing multiple protonation states and their charge parameters
	* MOLECULE_NAME_pH_X.Y-charged_output.mol2 
	* MOLECULE_NAME_pH_X.Y-charged_output.pdb

* Data sources from Ligand Expo
	* MOLECULE_NAME_pH_X.Y-rcsb_download.pdb
	* MOLECULE_NAME_pH_X.Y-rcsb_download.sdf

* Processed input for Epik
	* MOLECULE_NAME_pH_X.Y-before_epik.mol2
	
* Raw output from Epik
  * MOLECULE_NAME_pH_X.Y-epik.log
  * MOLECULE_NAME_pH_X.Y-epik.mae
  * MOLECULE_NAME_pH_X.Y-after_epik.mol2
  * MOLECULE_NAME_pH_X.Y-after_epik.sdf
  * MOLECULE_NAME_pH_X.Y-state-penalties.out

* Procedure log
	* log.txt 
