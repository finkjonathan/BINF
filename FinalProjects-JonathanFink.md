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
  - The virus was discovered in 2014 in Fujian, China and has been found all over the Guangdong Province in China. There is no news on recent outbreaks.

- **Host range:**  
  - Bamboo rat circovirus was found in a bamboo rat and there is no evidence that it has infected any other species, however circoviruses are known to infect many kinds of vertebrates. Many new circoviruses have been identified in recent years and they have become a large threat with farmed animals.

- **Cell entry:**  
  - Circoviruses penetrate the cell by binding to the hosts receptors and then perform receptor mediated endocytosis to get inside the cell. Once inside it replicates in the nucleus.

- **Replication strategy:**  
  - The virus replicates using a rolling-circle mechanism

- **Release mechanism:**  
  - If this virus follows other circoviruses, it likely uses the cell lysis release mechanism.

- **Latency:**  
  - There is no latency in the host cells and like most circoviruses, it replicates actively. 

- **Equilibrium and antigenic shift:**  
  - Bamboo rat Circoviruses is not at an equilibrium with humans and only affects bamboo rats and other rodents. An antigenic shift is not likely to happen as the virus has a low mutation rate. 

- **Vaccines:**  
  - As of now there are no vaccines created to counteract bamboo rat circovirus. 

- **Antiviral drugs:**  
  - There is no information on any antiviral drugs that have been created.

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

**Hydrophobicity Plot**

<img src="hydrophobicityNEWVIRUScomparison.png" alt="Hydrophobicity comparison" width="600"/>

In this image we plot the hydrophobicity of the bamboo circovirus against an E.coli proteome. We can see that the hydrophobicity sits roughly around 0 and is fairly simular to the E.coli, however we only have 5 points plotted of the circovirus versus the 7,767 plotted for the E.coli meaning ours is a less accruate measure.

-----------------------------------------------------------------------------------------------------------------------------------
(10 points) Identify any of your outlier hydrophobic proteins and BLAST them. What did these sequences annotate as? Might they have any important function to viral entry? 

-----------------------------------------------------------------------------------------------------------------------------------
**Bamboo Circovirus Compared to other Viruses:**


<img src="viral_genome_histogramMF.png" alt="Genome size compared to other viruses" width="600"/>

The Bamboo Circovirus falls roughly in the middle of the pack when compared to other viruses based on its genome size of 2112bp. 

-----------------------------------------------------------------------------------------------------------------------------------

**Phylogeny Tree**

<img src="FigTree.png" alt="FigTree plot, relatives" width="600"/>

For my tree I used the L-INS-i model to best fit my data. It stands for local pairwise alignment, iterative refinement, incorporating structural information. In MAFFT, it is one of the more accurate models, but the downside is how slow it is and computationally intensive. We can see our virus with its accession code: MF497827. Our outgroup is the MH782429.1 (Bovdisa virus) rooted at the top. Its three closest relatives are the dipodfec virus, rodent circovirus, and sonfela circovirus 1. There is no evidence of a host switch and the bootstrap values are not good at all.

-----------------------------------------------------------------------------------------------------------------------------------

# References Cited

Use at least five references, formatted using a numbering system (1) or (Doe, et al. 2024). List the references in the appropriate bibliography format here. Use markdown to have all of your text references link to the bottom part of the page to where the reference is listed.
