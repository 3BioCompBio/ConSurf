# ConSurf Python implementation

This `consurf.py` file provides a Python implementation of the famous ConSurf pipeline, which uses the [rate4site algorithm](https://pubmed.ncbi.nlm.nih.gov/15201400/) to predict the conservation/evolutionary rate at each position of a protein sequence. Briefly, the pipeline consits in first looking for homologous sequences, clustering them, generting a Multiple Sequence Alignement (MSA), and finally deriving the conservation score at each position in the MSA using rate4site.

When comparing the evolutionary scores predicted by our implementation with the ones predicted by the ConSurf website on a set of 30 randomly selected proteins, we obtained a Pearson correlation of 0.95. The differences are the result of the database version used as these will result in slightly different lists of homologous sequences.      

## Dependencies

The ConSurf pipeline has 4 major dependencies: `phmmer`, `cdhit`, `mafft` and `rate4site`. They are all easily installable using `apt-get install`. If these are not in your PATH, you should modify the paths in `soft_exec_paths` at the begining of the file. In addition, you will need a sequence database such as NR or UniRef90 in fasta format in which the homologous sequences will be searched.

## Usage

It only requires the input sequence as a FASTA file and the path to the sequence database, so the command looks something like

```
python consurf.py /home/user/Downloads/sequence.fasta -database=/home/user/Downloads/NR.fasta
```

To see all available parameters, run `python consurf.py -h`