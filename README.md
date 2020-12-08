This repository provides data and pipeline to create Fig.2C in the paper.

To plot the figure, simply change directory to ```GND_16S_23S``` and use the script ```plot.R```. You will need the packages ggplot2, scales, and chgpt.

To rerun the pipeline, change directory to ```GND_16S_23S``` and run

```bash
  bash full_pipeline.sh
```

Note that you will need to install the following packages:
+ TreeN93 <https://github.com/niemasd/TreeN93>
+ tqDist 
+ TreeSwift
+ Newick utilities
