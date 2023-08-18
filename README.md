## ReadMe
This repository contains the code and analysis necessary to reproduce the project I'm currently referring to as `MetaCon`. The goal of the project is to simulate mutualistic network assembly processes across a series of metapatches, and then see how either landscape structural properties or dispersal patterns affect which species persist through time. I'm especially interested in parameter spaces that allow for the coexistance of generalist and specialist species at the same time. 

The original model was largely inspired by the one-patch system presented by Becker et al. in their 2022 American Naturalist paper (https://doi.org/10.1086/720421), but has been changed in a number of ways during its development to better suit our goals. 

The general repository structure is as follows: (Note, this is not all files, but at least points you to some of the most important ones to start getting oriented)

```{bash}

├── README.md. \
├── .gitignore \
├── Analysis \
│   ├── DispersalSimulationAnComp.R: The meat of our analysis, this is a large function set up to run a simulation of our model. See documentation for details; right now this defaults to returning some summary outputs, but variations (such as `BigDispersalSimulationAnComp.R`) return different objects. The general workflow is to source this function into other rmd files, where we specify paramaters/landscapes to simulate across \
│   ├── Simulate_Across_Modularity.R: This file constructs our landscapes across a gradient of modularity (method details contained within), and then sources the above simulation function to simulate dynamics across them. This file utilizes the `parallel` package to paralalelize the simulations aross all but two cores of the machine it is running on, but be mindful of your machine\'s memory usage, especially when saving larger output objects. \
├── Data \
│   ├── This folder is ignored by default by the .gitignore, as the data object created from our simulations are almost always to big for simulation. However, with sufficient computing power you can recreate our analyses by setting the seed to be the same. \
├── Figures \
├── Manuscript \
│   ├── Also empty right now, as manuscript development is occuring on Overleaf. However, once the manuscript is developed some accessible version will live here. \

```

## License

![CC-BYNCSA-4](https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png)

Unless otherwise noted, the content in this repository is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-nc-sa/4.0/).
