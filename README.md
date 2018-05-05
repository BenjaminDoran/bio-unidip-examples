# BioDip Testing Examples

## Cloning Repo

```
git clone <repo url>
git submodule init
git submodule update
```

## Setup

### PYTHON 

We recommend using the Anaconda 3.6 version of Python with these packages:

BioPython
Numpy
Pandas
Matplotlib


### MUSCLE

Go to https://www.drive5.com/muscle/
to download the MUSCLE alignment tool for your platform. unzip the folder and place the executable in the base folder of this repo. All notebooks that use MUSCLE have a variable that holds the path to the executable.

### MEME

`runmeme.sh`  starts and runs a docker conatiner with meme installed. The file also contains an example command on how to run MEME on the raw sequences. Adjust the input file and output folder to run on trimmed sequences.

## Examples

- **numeric.ipynb** provides examples of generated random numeric samples.
- **letters.ipynb** provides examples looking at generated symbolic data.
- **aligned.ipynb** provides examples of isolating motifs after global alignments.
- **real_data.ipynb** retrives sequences from JASPER and NCBI
- **real_data_align.ipynb** provides examples of UniDip applied to real genomic sequences.
- **util.py** provides helper functions for data transforms and string generation.

