\# Sussex eDNA Metabarcoding Pipeline



An automated environmental DNA (eDNA) metabarcoding workflow that processes raw NCBI Sequence Read Archive (SRA) files into taxonomically annotated species tables.



The pipeline performs:



\- SRA → FASTQ conversion

\- Primer trimming (Cutadapt)

\- Sequence dereplication (USEARCH)

\- UNOISE3 denoising

\- NCBI BLAST identification

\- NCBI Taxonomy annotation



\---



\## Requirements



\### Software



\- Python 3.11+

\- SRA Toolkit

\- Cutadapt

\- USEARCH v12+

\- BLAST+

\- Internet connection (for NCBI queries)



\### Python packages



```bash

pip install -r requirements.txt

```



\---



\## Installation



Clone the repository



```bash

git clone https://github.com/YOUR\_USERNAME/YOUR\_REPOSITORY.git

cd YOUR\_REPOSITORY

```



Make the scripts executable



```bash

chmod +x 00\_sra\_to\_fastq.sh

chmod +x run\_pipeline.sh

```



\---



\## Input



Create a folder named



```

SR normalised Clark et al 2024

```



Place all downloaded `.sra` files inside.



\---



\## Run the pipeline



```bash

./run\_pipeline.sh

```



The pipeline automatically executes:



```

Stage 0  SRA → FASTQ

Stage 1  Primer trimming

Stage 2  Dereplication

Stage 3  Denoising

Stage 4  BLAST annotation

Stage 5  Taxonomy annotation

```



If a stage fails, execution stops and reports which stage encountered the error.



\---



\## Output



Each analysis creates a timestamped run directory inside



```

runs/

```



The final annotated species table is



```

master\_species\_composition.csv

```



To open the newest results directly from the terminal:



```bash

xdg-open "$(find runs -name master\_species\_composition.csv | sort | tail -n1)"

```



Or simply display its location:



```bash

find runs -name master\_species\_composition.csv | sort | tail -n1

```



\---



\## Repository Structure



```

00\_sra\_to\_fastq.sh

01\_trim.py

02\_dereplicate.py

03\_denoise.py

04\_blast.py

05\_taxonomy.py

run\_pipeline.sh

requirements.txt

README.md

```



\---



\## License



MIT License



\---



\## Author



Developed as part of an undergraduate research project at the University of Sussex.

