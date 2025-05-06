---
Final: "Bamboo Circovirus Sequence Analysis Project"
Author: "Jonathan Fink"
Date: "05/05/2025"
Accession: "MF497827"

---

# Introduction

<img src="BambooRatImage.jpg" alt="Bamboo Rat" width="300"/>

- **Viral classification:**  
  - *ICTV classification*: [Insert ICTV classification here with citation].  
  - *Baltimore classification*: Class 2, single stranded DNA virus.

- **Physical size:**  
  - There was no physical size in nanometers that I was able to find online, however it was 2112 base pairs and most other circoviruses are 15 to 25 nanometers long which compared to the human cell which is much smaller and covid is roughly 50 to 140 nm which makes it smaller as well..

- **Shape and envelope:**  
  - The virus is circular in shape and has no envelope.

- **Discovery and outbreaks:**  
  - The virus was discovered in 2014 in Fujian, China and there is no information on when the last outbreak occurred.

- **Host range:**  
  - Bamboo rat circovirus was found in a bamboo rat and there is no evidence that it has infected any other species, however circoviruses are known to infect many kinds of vertebrates.

- **Cell entry:**  
  - Circoviruses penetrate the cell by binding to the hosts receptors and then perform receptor mediated endocytosis to get inside the cell. Once inside it replicates in the nucleus.

- **Replication strategy:**  
  - [Your virus] [has its own replication machinery/relies on host machinery] and replicates by [describe process] [insert citation].

- **Release mechanism:**  
  - Viral progeny are released by [cell lysis, budding, excretion, etc.] [insert citation].

- **Latency:**  
  - [State if the virus shows latency and describe, or state if not] [insert citation].

- **Equilibrium and antigenic shift:**  
  - [Discuss if the virus is in equilibrium with humans and if it shows antigenic shift] [insert citation].

- **Vaccines:**  
  - [State if vaccines are available, type/mechanism of vaccine] [insert citation].

- **Antiviral drugs:**  
  - [List available antivirals and their mechanisms] [insert citation].

# Methods

1. **First, I downloaded the viral sequence by accession number, and selected XXX close relatives to identify a most recent common ancesstor**  Use the provided accession number (Column X) to download the viral genome sequence. Show the python code for each of the functions you performed.

```python
#Insert python code after each step
```
...Follow the additional methods instructions in the text. Make sure to include your code where relevant! For example, if you used a slurm script, make sure to paste it in here.

e.g. 
Align your sequences using the MAFFT slurm script
```bash
#!/bin/bash
#SBATCH --job-name=mafft_align
mafft --auto input.fasta > aligned.fast
```
     

# Results and Discussion

Present your results, including tables, figures, and summary statistics. This should be written more like a text document. Make sure you follow the instructions about which plots to paste within here and reference them as figures (e.g., Fig. 1) when you are describing them.

(10 points) Hydrophobicity plot agains the E.coli proteome

<img src="hydrophobicityNEWVIRUScomparison.png" alt="Hydrophobicity comparison" width="600"/>
(10 points) Identify any of your outlier hydrophobic proteins and BLAST them. What did these sequences annotate as? Might they have any important function to viral entry? 
(10 points) Plot the genome size of your virus relative to other viruses (see code from lab 9.12b). 
(20 points) Phylogeny and model selection. Use figtree to root the tree to your outgroup, make it look nice by ordering the nodes, increasing tip label size. State and interpret the best fit model used to infer this phylogeny. Display the bootstrap values. Discuss about your results. What are the three closest relatives of your virus, does it suggest a host switch? Are the branches well-supported by bootstrap values?
<img src="FINALTREE.pdf" alt="FigTree plot, relatives" width="600"/>


# References Cited

Use at least five references, formatted using a numbering system (1) or (Doe, et al. 2024). List the references in the appropriate bibliography format here. Use markdown to have all of your text references link to the bottom part of the page to where the reference is listed.
