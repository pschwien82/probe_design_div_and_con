### Dependencies:
The code has three dependencies that are easily installed with `pip`.
#### AnyTree for using the binary tree and rending it into ASCII
```bash
$ pip install anytree
```
see https://anytree.readthedocs.io/en/2.8.0/installation.html

#### Biopython for reading the fasta file
```bash
$ pip install biopython
``` 
see https://biopython.org/wiki/Download

#### MatplotLib for plotting probe histogram data
```bash
$ pip install matplotlib
```
see https://matplotlib.org/stable/users/installing.html

### Example invocations:
```bash
$ python3 design_probes_div_and_con.py --help
usage: design_probes_div_and_con.py [-h] [--k K] [--verbose] proteome

positional arguments:
  proteome    path to fasta file with all protein sequences to design probes
              for

optional arguments:
  -h, --help  show this help message and exit
  --k K       set the probe's target size in amino acids
  --verbose   print progress information during divide & conquer
```
For the proteome input file, I used [UniProt's latest reference proteome](https://www.uniprot.org/proteomes/UP000005640) including isoforms available at https://www.uniprot.org/uniprot/?include=true&format=fasta&force=true&query=proteome:UP000005640
```bash
$ wget "https://www.uniprot.org/uniprot/?include=true&format=fasta&force=true&query=proteome:UP000005640" > UP000005640_plus_iso.fasta

$ python3 design_probes_div_and_con.py UP000005640_plus_iso.fasta --k 3
Reading protein file...
Read 5000 protein sequences.
Read 10000 protein sequences.
Read 15000 protein sequences.
Read 20000 protein sequences.
Read 25000 protein sequences.
Read 30000 protein sequences.
Read 35000 protein sequences.
Read 40000 protein sequences.
Read 45000 protein sequences.
Read 50000 protein sequences.
Read 55000 protein sequences.
Read 60000 protein sequences.
Read 65000 protein sequences.
Read 70000 protein sequences.
Read 75000 protein sequences.
Read 80000 protein sequences.
Read 85000 protein sequences.
Read 90000 protein sequences.
Read 95000 protein sequences.
Read 97766 proteins.

Checking AA alphabet (non-canonical AAs will be ignored during probe design)...
M Methionine (Met) - 2.1799% (852634)
A Alanine (Ala) - 6.9282% (2709818)
W Tryptophan (Trp) - 1.2221% (477995)
S Serine (Ser) - 8.4291% (3296870)
N Asparagine (Asn) - 3.5950% (1406101)
Q Glutamine (Gln) - 4.8023% (1878300)
V Valine (Val) - 6.0189% (2354160)
T Threonine (Thr) - 5.4754% (2141596)
E Glutamic Acid (Glu) - 7.1940% (2813778)
F Phenylalanine (Phe) - 3.5458% (1386856)
I Isoleucine (Ile) - 4.3363% (1696032)
L Leucine (Leu) - 9.8363% (3847276)
R Arginine (Arg) - 5.6472% (2208768)
G Glycine (Gly) - 6.5096% (2546084)
Y Tyrosine (Tyr) - 2.6020% (1017715)
P Proline (Pro) - 6.3269% (2474642)
H Histidine (His) - 2.5809% (1009473)
D Aspartic Acid (Asp) - 4.8265% (1887802)
K Lysine (Lys) - 5.7527% (2250047)
C Cysteine (Cys) - 2.1696% (848596)
Non-canonical: U - 0.0001% (58)
Non-canonical: X - 0.0212% (8273)

Thu Apr  1 20:51:45 2021  Counting kmers...
done in 35.7 seconds, 8523 unique 3-mers found (80.04% of 22^3 possible)

Thu Apr  1 20:52:21 2021  Removing non-canonical kmers...
Eliminated the following 523 kmers due to containing non-canonical amino acids:
['SLU', 'LUG', 'UGT', 'CGU', 'GUK', 'UKL', 'GAU', 'AUG', 'UGY', 'LNU', 'NUS', 'USS', 'SCU', 'CUG', 'UGS', 'SPX', 'PXG', 'XGK', 'GLU', 'LUS', 'USY', 'GGU', 'GUG', 'UGR', 'VTU', 'TUS', 'TAX', 'AXT', 'XTC', 'GUQ', 'UQL', 'ASU', 'SUY', 'UYL', 'RSU', 'SUC', 'UCC', 'ITU', 'TUQ', 'UQC', 'CSU', 'SUQ', 'UQG', 'CQU', 'QUR', 'URL', 'AAU', 'AUQ', 'UQI', 'SUR', 'URU', 'RUK', 'UKN', 'KKU', 'KUE', 'UEU', 'EUP', 'UPS', 'VSU', 'SUG', 'KGU', 'UGC', 'TSU', 'URV', 'GCU', 'CTU', 'TUP', 'UPP', 'SYU', 'YUG', 'UGL', 'SDU', 'DUL', 'ULG', 'FPX', 'PXN', 'XNM', 'MPX', 'PXL', 'XLA', 'AYU', 'SQU', 'QUG', 'UGK', 'WCU', 'CUI', 'UIL', 'SSU', 'URS', 'HDU', 'DUA', 'UAR', 'RFU', 'FUI', 'UIF', 'ATU', 'KRU', 'UKK', 'ILU', 'LUL', 'ULV', 'XSY', 'XNT', 'XDL', 'XAR', 'XTK', 'XIQ', 'XVT', 'XGL', 'XEG', 'XQM', 'XLH', 'XLK', 'XER', 'XVY', 'XKK', 'XPA', 'XGG', 'XDP', 'XLR', 'XSK', 'XHL', 'XAE', 'XNY', 'XKS', 'XKE', 'XTL', 'XLL', 'XKL', 'XSQ', 'XAT', 'XRR', 'XFI', 'XAD', 'XRT', 'XCG', 'XAP', 'XSE', 'XKR', 'XPQ', 'XDS', 'XVC', 'XQR', 'XQL', 'XWS', 'XDG', 'XLT', 'XLE', 'XFG', 'XSL', 'XKG', 'XST', 'XNI', 'XGH', 'XEM', 'XEL', 'XTA', 'XWQ', 'XSC', 'XKC', 'XGY', 'XAQ', 'XED', 'XSN', 'XSS', 'XKI', 'XQE', 'XFE', 'XLY', 'XSD', 'XAA', 'XTW', 'XII', 'XMF', 'XRI', 'XPD', 'XNR', 'XID', 'XPS', 'XKT', 'XTP', 'XFS', 'XDK', 'XPP', 'XLI', 'XGW', 'XSP', 'XRF', 'XTT', 'XSR', 'XNL', 'XNN', 'XEV', 'XKH', 'XHH', 'XEH', 'XNA', 'XLV', 'XKN', 'XVS', 'XRG', 'XGE', 'XRA', 'XLQ', 'XVH', 'XIS', 'XLN', 'XKA', 'XGV', 'XGI', 'XNF', 'XNG', 'XFV', 'XLG', 'XFF', 'XLS', 'XCK', 'XFA', 'XVK', 'XHT', 'XDR', 'XHR', 'XGR', 'XDV', 'XRS', 'XYR', 'XGA', 'SAX', 'AXL', 'XLD', 'XPG', 'XAL', 'XVE', 'XYI', 'XGS', 'XEI', 'XIY', 'XYF', 'XVP', 'XCL', 'XPI', 'XTQ', 'XLC', 'XRD', 'XSV', 'XMH', 'XQD', 'XRV', 'XMS', 'XAS', 'XWD', 'XEN', 'XDY', 'XWL', 'XHG', 'XMP', 'XME', 'XPR', 'XML', 'XLP', 'XGC', 'XYY', 'XHA', 'XQH', 'XNS', 'XVN', 'XNH', 'XQS', 'XES', 'XIF', 'XTI', 'XVA', 'XIN', 'XVG', 'XEE', 'XAI', 'XAF', 'XET', 'XRP', 'XYA', 'XCE', 'XYE', 'XTG', 'XFL', 'XHI', 'XPK', 'XLF', 'XEK', 'XAG', 'XQY', 'XTV', 'XSH', 'XFP', 'XMN', 'XVV', 'XCV', 'XYK', 'XAM', 'XYS', 'XMM', 'XCA', 'XMR', 'XEY', 'XGT', 'XHE', 'XTD', 'XFY', 'XRK', 'XQG', 'XRQ', 'XTH', 'XAW', 'XYV', 'XRH', 'XPL', 'XIL', 'XRE', 'XTR', 'XAH', 'XPE', 'XIE', 'XKV', 'XGQ', 'XMA', 'XCR', 'XQQ', 'XAK', 'XFR', 'XYT', 'XQN', 'XEA', 'XQI', 'XCQ', 'XMT', 'XQP', 'XDD', 'XQT', 'XDW', 'XEC', 'XGD', 'XHV', 'XRY', 'XEQ', 'XTY', 'XVR', 'XNK', 'XDF', 'XPY', 'XIK', 'XYP', 'XDA', 'XSG', 'XWY', 'XGP', 'XYD', 'XYG', 'XVL', 'XPV', 'XAV', 'XQV', 'XVQ', 'XIP', 'XVF', 'XWT', 'XRL', 'XQF', 'XVD', 'XTM', 'XYL', 'XWG', 'XQK', 'XGF', 'XLM', 'XGM', 'XCC', 'XIG', 'XHF', 'XVW', 'XMI', 'XSI', 'XPF', 'XSA', 'XIC', 'XYQ', 'XDI', 'XPM', 'XTS', 'XDC', 'XNP', 'XIR', 'XPW', 'XKY', 'ISX', 'SXL', 'XMK', 'XIH', 'XFT', 'XFH', 'XNC', 'XCS', 'XNQ', 'XCT', 'XIA', 'XCM', 'XDN', 'XMD', 'XCF', 'XWM', 'XCP', 'XFQ', 'XTN', 'XRM', 'XQA', 'XHS', 'XKF', 'XVI', 'XYM', 'XKQ', 'XAY', 'XIT', 'XTF', 'XNV', 'XNE', 'XGN', 'XAN', 'XHQ', 'XHY', 'XDE', 'XSM', 'XWK', 'XPN', 'XEF', 'XND', 'XKM', 'XKD', 'XFK', 'XTE', 'XFN', 'XMC', 'XWN', 'XCW', 'XMG', 'XFM', 'XYW', 'CRX', 'RXL', 'XHK', 'XHP', 'XSF', 'XCD', 'XPH', 'XCI', 'XPT', 'XIM', 'XHN', 'XDT', 'XCN', 'XRN', 'XDQ', 'XWR', 'XIV', 'XWP', 'XRW', 'XEW', 'XYN', 'XNW', 'XWA', 'XDM', 'XAC', 'XFD', 'XKP', 'XEP', 'XSW', 'XHW', 'XHM', 'XWV', 'XIW', 'XMQ', 'PVX', 'VXL', 'XMV', 'XCY', 'XPC', 'XQC', 'XHC', 'XFW', 'XDH', 'XWC', 'XLW', 'IDX', 'DXI', 'RRX', 'RXA', 'XCH', 'XRC', 'XVM', 'SVX', 'XWI', 'XFC', 'XYC', 'LQX', 'QXP', 'AXP', 'XHD', 'XQW', 'XKW', 'QSX', 'DLX', 'LXL', 'KKX', 'KXR', 'PAX', 'XWF', 'FAX', 'AXR', 'XMW', 'XYH', 'XWW', 'DKX', 'HAX', 'AXQ', 'XMY', 'XWE']
done in 0.0 seconds

Thu Apr  1 20:52:21 2021  Building protein lookup table...
done in 39.6 seconds

Thu Apr  1 20:53:03 2021  Calculating probes...
done in 6.3 minutes

Thu Apr  1 20:59:23 2021  Validating results...
done in 3.3 seconds

Writing results...

A total of 4761 unique (2375157 total) kmers have been used to uniquely identify 96925 of 97766 proteins (841 could not be uniquely resolved).
Average number of probes to uniquely identify a protein is 24.51.

Writing ASCII tree...

Writing Newick tree...
```
### Output files:

#### Primary result file
    {input-filename}_{date}_{k}-mer_run_results.txt

The output format looks like this:

    [...]
    29577        +LLL,+SLL,+ALL,+ASL,-VSL,-LLT,-SSS,+LLE,-SSL,+LLG,+LVP,-FLR,+TVL,-SMS,+LHS
    29578        +LLL,+SLL,+ALL,+ASL,-VSL,-LLT,-SSS,+LLE,-SSL,+LLG,+LVP,-FLR,+TVL,+SMS,-QVA
    29576        +LLL,+SLL,+ALL,+ASL,-VSL,-LLT,-SSS,+LLE,-SSL,+LLG,+LVP,-FLR,+TVL,+SMS,+QVA,-RGE
    16371;16364  +LLL,+SLL,+ALL,+ASL,+VSL,+LST,+FLR,+LPK,+AEQ,+PQK,+IKL,+YLR,-MLG,-LAF,-FLC
    [...]
The first column is the sequence ID of the protein (ID is the sequence number `[0-n[` as read in order from fasta input file) and the comma separated unique probe set for this sequence.
The +/- sign in front of the AA triplet (single AA code) symbolizes whether or not this probe will or will not cause a signal when hybridized with the sequence.
In rare cases, multiple semi-colon separated sequence IDs are present, this indicates that these sequences could not be further resolved (see example row #4).


#### Binary tree in newick format.
    {input-filename}_{date}_{k}-mer.newick

Internal nodes are named with AA triplet and +/- sign. Leafs are one (or rarely multiple "-"-separated) sequence IDs (ID is the sequence number `[0-n[` as read in order from fasta input file). There are multiple online tools out there to visualize newick trees, however I have yet to find one that can handle a tree of the full human genome with k>=3 smoothly. [ITOL](https://itol.embl.de/) and [phylo.io](https://phylo.io/) work somewhat and so does [NCBI's Tree Viewer](https://www.ncbi.nlm.nih.gov/projects/treeview/) but their features are limited (e.g. can't label internal nodes). [Iroki](https://www.iroki.net/viewer) has the best features but can only visualize results from running a small subset of the whole proteome for debugging purposes.

#### ASCII rendered representation of the binary tree
    {input-filename}_{date}_{k}-mer_tree.txt

This file is not produced if the system's recursion limit is reached (k>=4 in my case)

it looks like this:
```
root
|-- +LLL
|   |-- -SSS
|   |   |-- +EEE
|   |   |   |-- -LAL
|   |   |   |   |-- -LGL
|   |   |   |   |   |-- -SLL
|   |   |   |   |   |   |-- -AAA
|   |   |   |   |   |   |   |-- +PPP
|   |   |   |   |   |   |   |   |-- -LLA
|   |   |   |   |   |   |   |   |   |-- -LSL
|   |   |   |   |   |   |   |   |   |   |-- -SAL
|   |   |   |   |   |   |   |   |   |   |   |-- -KAK
|   |   |   |   |   |   |   |   |   |   |   |   |-- -RRL
|   |   |   |   |   |   |   |   |   |   |   |   |   |-- -EKL
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -KKL
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -LSS
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -LRE
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -RAL
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -KVA
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- +EVE
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -LLV
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -EGK
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -STV
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -LIL
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -TRL
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -KEL
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -TLP
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -CCR
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -RRG
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -LLS
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -FFG
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -RRS
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -EER
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -EKE
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -NRR
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -GKG
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -FDR
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |-- -IKL
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   +-- 750
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   +-- +IKL
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |       +-- 795
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   +-- +FDR
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |       +-- 112
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   +-- +GKG
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |       |-- -LRG
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |       |   +-- 3002
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |       +-- +LRG
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |           +-- 1586
[...]
```
