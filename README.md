
1. Uses [`Python 2.7`](https://www.python.org/download/releases/2.7/)

2. Setup BLAST: follow installation and configuration steps from [`NCBI BLAST Manual`](https://www.ncbi.nlm.nih.gov/books/NBK279671/)

4. Have [`PyPI`](https://pypi.org/) and run the following to install requirements:
```
pip install -r requirements.txt
```

5.
Doing 'python run.py' prints instructions for how to use this program:
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


