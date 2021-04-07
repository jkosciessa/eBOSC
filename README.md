# eBOSC: extended Better OSCillation Detection

[![DOI](https://zenodo.org/badge/138600886.svg)](https://zenodo.org/badge/latestdoi/138600886)

## Overview
--------

**eBOSC** (extended Better OSCillation detection) is a toolbox (or a set of scripts) that can be used to detect the occurrence of rhythms in continuous signals (i.e., at the single trial level). It uses a static aperiodic ‘background’ spectrum as the basis to define a ‘power threshold’ that continuous signals have to exceed in order to qualify as ‘rhythmic’. As such, it leverages the observation that stochastic components of the frequency spectrum of neural data are aharacterized by a '1/f'-like power spectrum. An additional ‘duration threshold’ can be set up in advance, or rhythmic episodes can be filtered by duration following detection to ensure that detected rhythmic episodes have a rather sustained vs. transient appearance.

## Documentation
-------------

A project wiki for eBOSC is available [here](https://github.com/jkosciessa/eBOSC/wiki).

* [Motivation](https://github.com/jkosciessa/eBOSC/wiki/Pitfalls)
* [Tutorial](https://github.com/jkosciessa/eBOSC/wiki/Tutorial)
* [Version update/Legacy information](https://github.com/jkosciessa/eBOSC/wiki/Legacy)

Simulation scripts and data files regarding the 2020 NeuroImage paper can be found at https://github.com/jkosciessa/eBOSC_resources_NI2020.

## Installation
-------------

After downloading, simply add the toolbox path. The toolbox requires a version of MATLAB and has been tested with R2017b.

Note: If you use the 'Download ZIP' button to retrieve the repository, the .mat example file may not be downloaded due to a problem in GitHub's lfs implementation. Using 'git clone' to copy the repository should work to retrieve the file.

## Credits
-------------

If you find the method useful, please cite the following papers:

Kosciessa, J. Q., Grandy, T. H., Garrett, D. D., & Werkle-Bergner, M. (2020). Single-trial characterization of neural rhythms: Potential and challenges. NeuroImage, 206, 116331. http://doi.org/10.1016/j.neuroimage.2019.116331

Whitten, T. A., Hughes, A. M., Dickson, C. T., & Caplan, J. B. (2011). A better oscillation detection method robustly extracts EEG rhythms across brain state changes: The human alpha rhythm as a test case. NeuroImage, 54(2), 860–874. http://doi.org/10.1016/j.neuroimage.2010.08.064

## License
-------------

eBOSC is an extension of the BOSC library and partially uses scripts from the original toolbox. These functions are included in the 'external' folder of the current package.

The eBOSC library (and the original BOSC library) are free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The eBOSC library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
