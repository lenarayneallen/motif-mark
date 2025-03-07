# motif-mark-oop.py:

This python script generates a png figure depicting the location of motifs (provided in a text file) within sequences (provided in a FASTA format). This script will work with up to 10 FASTA sequences (of no greater than 1000bp), and 5 motifs. The produced figure depicts introns as black lines and and exons as yellow boxes. The motifs are depicted as color-coded vertical boxes. All renderings of sequences, introns, exons, and motifs are to scale. However, overlapping motifs are not staggered in the output figures. 

[Script](./motif-mark-oop.py)

[FASTA input File](./Figure_1.fasta)

[Motif input file](./Figure_1_motifs.txt)

[Output figure](./Figure_1.png)


## Input File Specifications:
### FASTA file:
Headers should be formatted as: 
`>GENENAME chr#:start_position-end_position`
Where #, start_position, and end_position are all integers. 

Headers of sequences which are reverse complements should be formatted as: 
`>GENENAME chr#:start_position-end_position (reverse complement)`

In sequences, exons should be represented by capital bases (A,C,G,T), and introns should be represented by lowercase bases (a,c,g,t).


### Motifs file:

The Motifs file should be a txt file with one motif per line. Motifs are assumed to be case-insensitive.

For motifs containing uracil, al "U"s will be converted to "T"s. 

Motifs containing the [ambiguous nucleotides](https://en.wikipedia.org/wiki/Nucleic_acid_notation) W, S, M, K, R, Y, B, D, H, V, and N are acceptable and will be properly interpreted.  

## Usage:
```bash
#general usage

./motif-mark-oop.py -m <motif_file.txt> -f <fasta_file.txt>

#to run with test files in this repository
./motif-mark-oop.py -m Fig_1_motifs.txt -f Figure_1.fasta
```
## About the script:

### Classes: 
- InputFile
- Sequence
- Motif
- Exon
- DrawingMaterials

### Pycario

This python script uses the [`pycairo`](https://pycairo.readthedocs.io/en/latest/) package to generate the output figure. This output file will have the same prefix as the FASTA file (i.e. `figA.fasta` -> `figA.png`)