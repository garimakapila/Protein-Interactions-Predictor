# Protein Interactions Predictor

## Description
This repository contains the code I wrote for my independent research project when I worked with Yu Lab @ Cornell in Spring 2018.

This program transfers protein interaction annotations across organisms based on sequence similarity at various sites including interface sites, domain sites, or DNA binding sites. Sequences are compared by calculating features including # matching residues, # similar residues, or differences in residue features such as hydrophobicity, polarity, or solvent accessible surface area. You can find the rankings of these features in the [/P-Value](https://github.com/garimakapila/Protein-Interactions-Predictor/tree/master/P-Values) folder.

The steps are automated with messages and progress bars to indicate the status in the command line: data collection from the web, data filtering, calculations, clustering interactions, and generating a powerpoint file with visualizations.

Databases used:
* [`UniProt`](http://www.uniprot.org)
* [`HINT`](http://hint.yulab.org)
* [`Interactome INSIDER`](http://interactomeinsider.yulab.org)
* [`The European Bioinformatics Institute`](https://www.ebi.ac.uk)

Organisms used:
* Escherichia Coli, Strain K12
* Drosophila Melanogaster
* Mus Musculus
* Homo Sapiens
* Schizosaccharomyces Pombe, Strain 972H
* Caenorhabditis Elegans
* Arabidopsis Thaliana
* Saccharomyces Cerevisiae, Strain S288C

## Instructions
1. Have [`Python 2.7`](https://www.python.org/download/releases/2.7/)

2. Setup BLAST: follow installation and configuration steps from [`NCBI BLAST Manual`](https://www.ncbi.nlm.nih.gov/books/NBK279671/)

3. Have [`PyPI`](https://pypi.org/) and run the following to install requirements:
```
pip install -r requirements.txt
```

4. Running 'python run.py' prints instructions for how to use this program:
```
Enter two species to transfer interactions from, or one species to transfer
all other organism protein interactions from, or "arguments.txt" as input
arguments and an option for interfaces, examples:
> HS SC ALL
> HS HQ
> arguments.txt ALL

Organism options: 
	EC	Escherichia Coli, Strain K12
	DM	Drosophila Melanogaster
	MM	Mus Musculus
	HS	Homo Sapiens
	SP	Schizosaccharomyces Pombe, Strain 972H
	CE	Caenorhabditis Elegans
	AT	Arabidopsis Thaliana
	SC	Saccharomyces Cerevisiae, Strain S288C
Interface options: 
	ALL
	HQ

Enter "quit" to exit
Next time, you can directly run, for example, "python run.py HS SC HQ"

>
```

Or you can directly run any of the following examples:
```
python run.py HS SC ALL
python run.py HS HQ
python run.py arguments.txt ALL
```

## Results

There are a total of 56 possible organism-organism mappings. You can find a few examples in the [/Graphs](https://github.com/garimakapila/Protein-Interactions-Predictor/tree/master/Graphs) folder.

An example of a plot comparing the BLOSUM scores between non-interactacting (0) vs. interacting (1) pairs. 

![alt text](https://raw.githubusercontent.com/garimakapila/Protein-Interactions-Predictor/master/plot.png)

After picking the top features and doing dimensionality reduction, these were the PCA results:

![alt text](https://raw.githubusercontent.com/garimakapila/Protein-Interactions-Predictor/master/pca.png)


