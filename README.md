This repository provides data and pipeline to create Fig.2C in the paper. Please note that this pipeline developed and run on macOS, so there is some assumption about the OS system. The pipeline would run on Linux systems (perhaps with little modification), but will NOT work on Windows.

To plot the figure, simply change directory to ```GND_16S_23S``` and use the script ```plot.R```. You will need the packages ggplot2, scales, and chgpt.

To rerun the pipeline, you first need to install the following packages and make sure they are in your path:
+ [TreeN93](https://github.com/niemasd/TreeN93)
+ [tqDist](https://users-cs.au.dk/cstorm/software/tqdist/)
+ [treeswift](https://github.com/niemasd/TreeSwift)
+ [Newick utilities](http://cegg.unige.ch/newick_utils)

Then change directory to ```GND_16S_23S``` and run

```bash
  bash full_pipeline.sh
```

You should see many files created that look similar to the files in the ```Data.tar.gz```. Please note that you will NOT see the exact same files as in ```Data`.tar.gz`` because there is randomness in our algorithm. Nevertheless, you will get a similar figure.


