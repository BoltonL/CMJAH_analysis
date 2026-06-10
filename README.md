Below is a GitHub README drafted in a format appropriate for an epidemiology/infectious diseases analysis repository accompanying a peer-reviewed publication.

# Discordant Empiric Antimicrobial Therapy, Mortality and Weighted-Incidence Syndromic Combination Antibiograms (WISCA)-Estimated Coverage of Neonatal Healthcare-Associated Bloodstream Infections

## Overview

This repository contains the analysis scripts used for the manuscript:

**Discordant empiric antimicrobial therapy, mortality and weighted-incidence syndromic combination antibiograms (WISCA)-estimated coverage of neonatal healthcare-associated bloodstream infections**

### Study Aim

Neonatal healthcare-associated bloodstream infections (HA-BSIs) are a major cause of morbidity and mortality, particularly in settings with high levels of antimicrobial resistance (AMR). Empiric antibiotic regimens are commonly prescribed before microbiological results become available; however, discordance between empiric therapy and pathogen susceptibility may adversely affect clinical outcomes.

## Repository Contents

This repository contains scripts used for:

* Data cleaning and preparation
* Bloodstream infection episode classification
* Antimicrobial susceptibility processing
* WISCA generation and estimation of regimen coverage
* Manuscript figures/tables

## WISCA Methodology

The WISCA implementation scripts included in this repository represent precursor versions of functionality that subsequently informed development of the WISCA tools available within the AMR for R package.

Current WISCA functionality is available through:

**AMR for R package**
[https://amr-for-r.org/articles/WISCA.html](https://amr-for-r.org/articles/WISCA.html)

Researchers interested in applying WISCA methods are encouraged to explore the actively maintained implementation available through the AMR package.

## Software Requirements

Analyses were conducted using R.

Required packages may include (depending on script version):

* tidyverse
* AMR
* ggplot2
* lubridate

Please refer to individual scripts for package-specific requirements.

## Reproducibility

The analytical workflow is provided to promote transparency and reproducibility of the analysis methods used in the manuscript.

Because the original study data are not publicly available, execution of the complete pipeline will require access to appropriately structured datasets and associated data permissions.

## Ethical Considerations and Data Availability

The analysis scripts contained within this repository are publicly available under the MIT License.

However, the clinical and microbiological datasets used in the analyses are subject to:

* Institutional ethics approvals
* Data governance requirements
* Data sharing agreements with the relevant data custodians

Consequently, the underlying data cannot be publicly released through this repository.

Researchers interested in accessing the data should contact the relevant data custodians and obtain all necessary approvals before any data sharing can be considered.

## Related Resources

### AMR for R WISCA Documentation

[https://amr-for-r.org/articles/WISCA.html](https://amr-for-r.org/articles/WISCA.html)

## Citation

If you use these scripts or adapt the methodology, please cite:

> [Authors]. Discordant empiric antimicrobial therapy, mortality and weighted-incidence syndromic combination antibiograms (WISCA)-estimated coverage of neonatal healthcare-associated bloodstream infections. [Journal, Year].

## License

This repository is distributed under the MIT License.

See the LICENSE file for details.

## Contact

For questions regarding the analysis methods, WISCA implementation, or repository contents, please open a GitHub issue or contact the corresponding author (vindanac@nicd.ac.za).

