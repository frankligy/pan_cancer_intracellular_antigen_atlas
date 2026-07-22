# Validating vaccine candidates by Immunopeptidome (using Tesorai)

Leveraging large-scale publicly available immunopeptidomic data, one can rapidly search their antigen of interest for various applications. One of them will be validating the existence of cancer vaccine candidates.

<p align="center">
  <img src="../images/vaccine.png" alt="vaccine" width="400" height="400">
</p>


### Step 1: Fasta File

You shall have a list of antigen peptide sequence as fasta files, and let's say you are interested in melanoma (SKCM), you shall append this list with [our search space](https://genome.med.nyu.edu/public/yarmarkovichlab/ImmunoVerse/search_space_tesorai/SKCM/db_fasta_tesorai/). You don't need to use all aberrations, you can only use the canonical one (ensembl_protein.fasta) if desirable. Lastly, add [contaminants fasta](https://genome.med.nyu.edu/public/yarmarkovichlab/ImmunoVerse/database/contaminants_polish.fasta). 

Essentially, you just combine all fasta together, this can be easily done by manual copy-and-paste, or linux command `cat`, or your favorate coding language, let's assume after this step, you shall have a fasta named `combined_SKCM_pan_cancer.fasta`.


### Step 2: Upload Fasta file to Tesorai platform

We host all public raw files on [Tesorai platform](https://console.tesorai.com/) under `tesorai_immunoverse` project, you shall follow the instructions on creating an account. If you don't have access to the workspace, email me (li2g2@mail.uc.edu) or the tesorai support. 

Now, upload your fasta file to the project and SKCM workspace.

<p align="center">
  <img src="../images/tesorai_vaccine.png" alt="vaccine" width="500" height="500">
</p>


### Step 3: Run Tesorai

Now, just easily submit the job, and the search results will be sent to your email once finished, it takes about 1 hour!

<p align="center">
  <img src="../images/run_step1.png" alt="vaccine" width="500" height="500">
</p>

<p align="center">
  <img src="../images/run_step2.png" alt="vaccine" width="500" height="500">
</p>









