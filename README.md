# Bioelectrodynamics of Yeast Cells

This repository provides code implementation to evaluate bioelectrodynamic properties of yeast cells.

## Description
This project evaluates the bioelectric dynamic properties of a yeast cell model.
The yeast cell model iMM904, accesible in [BiGG Models repository](http://bigg.ucsd.edu/models/iMM904),  is included within the Matlab file `iMM904.mat`.
The model of the genes, metabolites and metabolic reactions within the cells.

The model is described with the following variables:

Related to the model:
- `description`: Used to store a free‐text summary of the model. 
This summary typically includes details such as: Organism or System (which organism or cell type the model represents), Version and Scope (The version of the model and what it covers), Purpose and Context (The intended use of the model, e.g., for metabolic engineering, gene essentiality analysis, etc.), References (Citations to publications or databases that describe or support the model).
- `modelVersion`: Specify the version of the model reconstruction.
Typically, it is a string that indicates a version number, revision date, or other identifier that distinguishes one version of the model from another.
This metadata is important for: Tracking Changes (It helps document updates, corrections, or improvements made to the model over time),Reproducibility (Users can reference the exact model version used in their analyses, ensuring that results can be reproduced or compared with other studies), Curation and Communication (It provides clarity on which iteration of the model is being used, which is especially useful when models are updated or when multiple versions exist).

Related to the metabolites:
- `mets`: cell array of strings that lists all the metabolite identifiers included in the model.
Each entry corresponds to one metabolite, and the identifiers often include information about the metabolite as well as its cellular compartment.
This field is fundamental because it defines the set of chemical species that participate in the metabolic reactions of the model.
For example, "AUXOSC10fthf_u[u]" which is the first entry stand for "10_Formyltetrahydrofolate".
The notation within brakets denotes the compartment, such as '[c]' for cytosol, '[e]' for extracellular, '[x]' denotes an exchange compartment as a pseudo‐compartment used for metabolites that are taken up from or secreted to the environment, '[u]' usually stands for unassigned or undefined (or sometimes universal), meaning that the metabolite is not assigned to one of the standard compartments (such as cytosol, mitochondrion, etc.).
In some models it may be used as a catch‐all compartment when the precise localization is not specified.

- `metNames`: Cell array of strings that provides the full, human‐readable names for each metabolite in the model.
Each entry in `metNames` corresponds to the metabolite identifier listed in the `mets` field, giving a clearer description of the metabolite (for example, "10_Formyltetrahydrofolate" instead of an abbreviated identifier like "AUXOSC10fthf").
This field is useful for model visualization, reporting, and ensuring that users can easily interpret the metabolic components during analysis or when presenting results.
The metabolites' names are also labeled with specific naming conventions such as
  - '_c': for metabolites in the cytosol,
  - '_e': for metabolites in the extracellular space, and
  -  '_m' for metabolites in the mitochondria


- `metFormulas`: Cell array of strings that contains the chemical formulas for each metabolite listed in the model's mets field.
These formulas—expressed in standard chemical notation (e.g., "C6H12O6" for glucose)—provide the elemental composition of each metabolite.
They are used for tasks such as checking mass and charge balance in reactions, performing thermodynamic calculations, and aiding in the annotation and curation of the model.
- `metCharges`: Numerical vector in a COBRA model where each element represents the net electrical charge of the corresponding metabolite (as listed in the model’s "mets" field).
This information is essential for ensuring that reactions are charge balanced and for supporting further thermodynamic or mass balance analyses within the model.
- `metHMDBID`: Cell array that contains the Human Metabolome Database (HMDB) identifiers for each metabolite in the model.
These unique identifiers allow users to cross-reference the metabolites in the COBRA model with the HMDB, facilitating the retrieval of detailed metabolite information—such as chemical structures, properties, and related biochemical data—from the external database.
- `metInChIString`:  Cell array of strings where each string is the IUPAC International Chemical Identifier (InChI) for a corresponding metabolite (as listed in the model's 'mets' field).
An InChI string is a standardized, computer-readable representation of a metabolite's chemical structure, encoding details like the molecular formula, connectivity, stereochemistry, and charge.
This information is useful for cross-referencing metabolites with external chemical databases, ensuring consistency in metabolic reconstructions, and supporting further chemical analyses.
- `metKEGGID`: Cell array of strings that contains the KEGG compound identifiers for each metabolite listed in the model s 'mets' field. 
These identifiers (such as "C00031" for glucose) link model metabolites to entries in the Kyoto Encyclopedia of Genes and Genomes (KEGG) database.
This cross-referencing facilitates model curation, helps ensure consistency with published metabolic data, and supports further pathway and functional analyses by connecting the model to a well-established resource for metabolic information.
- `metChEBIID`: Cell array of strings where each string is the ChEBI identifier corresponding to a metabolite in the model's 'mets' field.
These identifiers link the metabolites to entries in the Chemical Entities of Biological Interest (ChEBI) database, allowing users to retrieve detailed chemical information, such as molecular structure, synonyms, and chemical properties.
This cross-referencing is valuable for model curation, ensuring consistency with established chemical databases, and facilitating further analyses such as metabolic pathway mapping and network integration.
- `metPubChemID`: cell array of strings where each entry represents the PubChem Compound Identifier (CID) for the corresponding metabolite in the model's 'mets' field.
These identifiers link the metabolites to their records in the PubChem database, which provides extensive information about each compound—including molecular structure, chemical properties, and biological activity.
This cross-referencing aids in model curation and enables further integration with other biochemical databases and computational tools.

Related to metabolic reactions:
- `rxns`: List of reactions identifiers. 
It is a cell array of char streams.
- `rxnNames`: List of reactions names. 
It is a cell array of char streams.
- `lb`: Array of doubles for the lower bounds of the reactions.
- `ub`: Array of doubles for the upper bounds of the reactions.
- `S` matrix: the stoichiometric matrix
  - If a component of `S`, denoted as $S_{ij}$, is negative, metabolite $i$ is a reactant in reaction $j$, and the absolute value represents the amount consumed.​
  - $S_{ij}>0$ denotes metabolite $i$ is a product in reaction $j$, and the absolute value represents the amount produced.​
  - $S_{ij}=0$ denotes metabolite $i$ does not participate in reaction $j$.

- `b`: Array of doubles for the accumulation (positive) or depletion (negative) of the corresponding metabolites.
The value of '0' Indicates no concentration change.
- `csense`: List of indicators whether the `b` vector is a lower bound ('G'), upper bound ('L'), or hard constraint 'E' for the metabolites. It is a cell array of char.
- `c`: Array of double to define the linear objective.
- `osenseStr`: The objective sense either 'max' for maximisation or 'min' for minimisation.


Related to gene description:
- `genes`: List of genes within the model. 
- `rxnGeneMat`: Matrix form of gene–reaction links, as an associations that link genes (which encode enzymes) to metabolic reaction.
Ìt is a binary matrix with one row per reaction and one column per gene in the model​.
An entry of '1' in this matrix indicates that a given gene is associated with (i.e. participates in or is required for) the corresponding reaction, while '0' means no association.
- `rules`: The Gene-protein-reaction rules in a computer readable format present in your model.
- `grRules`: Describe how the genes combine to enable the reaction.
It is given as a logic rule.
For instance, if a reaction can be catalyzed by geneX OR geneY (isoenzymes), the row will have '1's for both genes, and the rule would be “geneX or geneY”


    
<figure>
    <p align="center">
        <img src="https://github.com/tkn-tub/semeco_explained/blob/main/figures/prototrophic_cell.png" alt="nn" width="200">
    </p>
</figure>
<p align="center">
Fig. 1: Illustration of a Protrophic cell used in the SeMeCo experiment.
</p> 

## Usage

This project directly runs from the file `A_Master_File.mlx` located in the folder code.
As this code uses some COBRA toolbox functions, it is needed first to install it in matlab.
Description for [installing COBRA in Matlab is in this link](https://opencobra.github.io/cobratoolbox/stable/installation.html).

## License
![Licence](https://img.shields.io/github/license/larymak/Python-project-Scripts)

## References
<a name="fn1">[1]</a>:Yu, J.S.L., Correia-Melo, C., Zorrilla, F. et al. Microbial communities form rich extracellular metabolomes that foster metabolic interactions and promote drug tolerance. Nat Microbiol 7, 542–555 (2022). 
[Link to paper in Nature Microbiology](https://doi.org/10.1038/s41564-022-01072-5)

## Contact Information

- **Name:** Jorge Torres Gómez

    [![GitHub](https://img.shields.io/badge/GitHub-181717?logo=github)](https://github.com/jorge-torresgomez)

    [![Email](https://img.shields.io/badge/Email-jorge.torresgomez@ieee.org-D14836?logo=gmail&logoColor=white)](mailto:jorge.torresgomez@ieee.org)

    [![LinkedIn](https://img.shields.io/badge/LinkedIn-torresgomez-blue?logo=linkedin&style=flat-square)](https://www.linkedin.com/in/torresgomez/)

    [![Website Badge](https://img.shields.io/badge/Website-Homepage-blue?logo=web)](https://www.tkn.tu-berlin.de/team/torres-gomez/)

- **Name:** Mohammad Tauqeer Alam

    [![Email](https://img.shields.io/badge/Email-mtalam@uaeu.ac.ae-D14836?logo=gmail&logoColor=white)](mailto:mtalam@uaeu.ac.ae)

    [![LinkedIn](https://img.shields.io/badge/LinkedIn-blue?logo=linkedin&style=flat-square)](linkedin.com/in/mohammad-tauqeer-alam-49a2116)

    [![Website Badge](https://img.shields.io/badge/Website-Homepage-blue?logo=web)](https://www.uaeu.ac.ae/en/cos/profile.shtml?email=mtalam@uaeu.ac.ae)
