# the COG-or

## Introduction
Package for improving the functional annotation of bacterial genomes, classification of protein-coding sequences into clusters of orthologous groups, and visualization 
of the final annotated genome. The package uses the outputs of the tools that assign COG to protein
coding sequences, namely eggNOG-mapper, Operon-mapper,and Batch CD-Search. The COG-or includes functions to 
process these outputs and improve the annotation. It outputs a new processed file in a suitable format 
that is ready to be visualized in DNAPlotter program. 

The tools for a genome annotation annotation are available at:
eggNOG-mapper:
Operon-mapper
Batch CD-Search:

## Installation 
Installers for the latest released version are available at the Python Package Index (PyPI):

```
pip install COGor
```

## Usage
The COG-or package includes functions which can be called as follows:

### PROGRAM PROCESSING

```
batch_splitter(organism_name,CDS_file)
batch_merger(organism_name,file1,file2)
batch_processor(organism_name,batch_file)
em_processor(organism_name, eggNOGmapper_file, CDS_file)
om_processor(organism_name, Operon_ORF_file, Operon_COG_file)
```

### ANNOTATION IMPROVEMENT
```
consensus(om_file,em_file,batch_file,fasta_file, get_pseudo=True, get_ncrna=True, gff_file)
```

Functions can be used in order that is shown in the diagram.

<img src="diagram.png" width="450" height="400">

### VISUALIZATION
```
track_manager()
get_legend()
```

After uploading the required file to DNAPlotter, user can obtain similar image as the one shown as an example below.

<img src="genome_map.png" width="600" height="350">
