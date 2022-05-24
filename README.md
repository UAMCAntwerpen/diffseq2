# Diffseq2: a Pymol plugin to visualise differeces between two proteins


## Description
Highlight sequence differences between two proteins. There can be no missing loops in the proteins, except at the beginning or end of the proteins.
Makes use of the python 'align' algorithm, this works best if the two proteins are highly similar.


maxdiff/lessdiff option controls what to do with parts of the sequence that are not considered aligned by the pymol alignment algorithm
maxdiff: resiudes that are not considered aligned by pymol, but have the same residue name are considered as different. (default)
lessdiff: resiudes that are not considered aligned by pymol, but have the same residue name can be considered as sufficiently different.

In the case of lessdif we identify three different secenarios:
        1) the non aligned block has the same length in A and  in B and has the exact same residue names
        2) the non aligned block has the same length in A and  in B and has some different residue names 
        3) the non aligned block has  not the same length in A and  in B

We have chosen to implement 'lessdiff' in such a way that only case 1 will not be included in the difference objects. Case 2 and 3 are seen as suffieciently different


## Usage

Pymol command line:

    diffseq2 proteinA, proteinB
    diffseq2 proteinA, proteinB, maxdiff 
    diffseq2 proteinA, proteinB, lessdiff
    
## Installation

One time use: place python script in same folder as where you launch pymol and enter 'load diffseq2.py' in the pymol command line.

Installation as plugin: Pymol --> Plugin --> Plugin Manager --> Install New Plugin --> Choose file and then select diffseq2.py file location.





diffseq2 is written by Olivier Beyens 
This standalone script builds further on Joao Rodrigues 'diffseq' script.
