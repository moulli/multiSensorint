# multiSensorint
Program that allows to compare multi-sensory brain scans, through the use of a grid. Comparison is made possible with a graphical interface. Loading a ZBG object allows then to chose subsets (e.g. visual or vestibular stimuli) and to compare their most active neurons. 


## The ZBG class
@ZBraingrid is a class that aims to facilitating manipulation of grids associated to datasets. 

### Defining a ZBG object
A ZBG object has three intrinsic values: 

* a *method* [character|string]: defining name of object,

* a *size* [array of integers, length 1 or 3]: vector of the number of voxels in the grid, along all 3 dimensions. If one wants to define grid based on voxel size, one can use grid dimensions to compute final grid size. Grid dimensions are ([0, 0.496], [0, 1.122], [0, 0.276]), reference brain from the [Z Brain Atlas](https://engertlab.fas.harvard.edu/Z-Brain/home/),

* an *orientation* [character]: indicating where axes (x, y, z) point to. Orientation should be a character (e.g. 'RAS' or 'LPI'), with R:right, L:left, A:anterior, P:posterior, S:superior, I:inferior.

### Adding a dataset to a ZBG object
The function addDataset allows to add a dataset to an already existing object. It takes two mandatory arguments and one optional argument. First argument is the ZBG object, second argument is the dataset structure (explained below), and third argument is to disable information printing while new data set is added (just type 'information_off' in this case). 

Dataset structure is a structure containing mandatory information on new dataset to add. As grid size is set when ZBG object is created, new data will already adjust this size. The fields required are:

* *name* [character|string]: defining name of new dataset,

* *path* [character|string]: defining path to actual data (HDF5 file for us),

* *comment* [character|string]: defining key words associated to dataset. Comment is important because it is used when discriminating ZBG subset from ZBG object, so chose wisely,

* *orientation* [character]: see orientation from section above. Orientation is automatically adjusted to be the same as the one from the ZBG object,

* *coordinates* [matrix of float|double, size (n x 3)]: matrix with as many rows as there are neurons to add (n), and 3 columns for the three coordinates along x, y, and z axes,

* *correlation* [array of float|double, length n]: vector with a length equal to the number of neurons from coordinates (n). Correlation does not have to be actual correlation, but rather a number linking neurons' signals to stimulus. We discuss this aspect in a following section. If user only wants to compute projection on grid, then they should enter random values.

When adding a new dataset to ZBG object, addDataset computes for each neuron the voxel to which it belongs. To each neuron is therefore attributed an index in the grid. Then, correlations are averaged for each voxel. It is this averaged correlation value that is later used in the graphical interface. **When we speak about correlation later, we refer to this averaged correlation**.

### ZBG methods 

We will not make an extensive list for all ZBG methods here, but rather describe the most important and useful ones. User can visualize all methods in the ZBG class Matlab file. We already went through addDataset, which allows to add new datasets to an already existing object. Basic operations are provided, such as addition (merge two ZBG objects), length (returns number of datasets in ZBG object) and indexing (get new ZBG object with datasets corresponding to indexes). Following are some of the other methods:

* *subset*: takes an additional input argument (key words, as a character). Returns a new ZBG object containing datasets whose comments include the key words. This method is later used to discriminate datasets in the graphical interface,

* *clean*: detects if a dataset occurs more than once, and delete potential doublons,

* *flatten*: returns a new ZBG object with only one dataset, which contains the averaged correlation of all datasets from former ZBG object,

* *downgrid*: takes an additional input argument (new number of voxels). Returns a new ZBG object with new grid size and datasets adapted to it,

* *getLabels*: returns labels vector for all voxels in ZBG object,

* *get3Dcoord*: returns coordinates matrix for all voxels in ZBG object,

* *scat3*: scatter3 of all voxels, 

* *corAnalysis*: boxplot for each dataset, allowing to visualize difference in correlations from one dataset to another.


## Graphical interface

### Design and controls

### Loading a ZBG object and selecting subsets

### Isovalues for correlation

### Common neurons between two datasets

### Plotting neurons' signals, and HDF5 architecture

### Visualizing brain regions

### Further analyses


## Discussion on methods

### Interfish reliability

This atlas does not have the pretention to compare zebra fish brains neuron by neuron. Depending on the fish, brain organisation can differ, even if brain functions are fairly the same. In order to certify that this way of comparing brains is relevant, we made this particular study.

For each stimulus (we had auditory, vestibular step, vestibular sine, hot and cold), we took the 1% of most active neurons, based on the F-statistic. Using the [Z Brain Atlas](https://engertlab.fas.harvard.edu/Z-Brain/home/), we deduced for each region (294 in total) the number of most active neurons part of it. By dividing this 294-long array by the total number of 1% most active neurons, we obtained a distribution, comparable across fish and stimuli.

We used the Kolmogorov-Smirnov test to compare these distributions, telling if two regions vectors follow the same law. Kolmogorov-Smirnov test checks maximum distance between the two cumulative distributions against a table. 

![Example of cumulative distributions for acoustic and thermotaxis stimuli](/README_img/cumul.png)

We used this test to compare all the experiments one to another for a given stimulus. In order to do that, we used the p-value for the Kolmogorov-Smirnov test between two regions vectors. The p-value gave the following information: it was the probability to find samples this extreme taken from the same distribution (the null hypothesis). The highest the p-value, the more possible null hypothesis was accepted. In the end, it was possible to make a boxplot of the distribution of p-values per stimulus, as plotted below.

![Distributions of p-values across stimuli](/README_img/boxplot.png) 

How to interpret these plots? Some stimuli have a great repeatability in terms of most active regions, some have a good repeatability with a wide range, and some have a poor repeatability with a wide range. When we take a closer look at the p-values inside a stimulus, we notice that there can be experiments with a low p-value, no matter what other experiment we compare it to. We could then delete these 'bad' experiments from the initial dataset, in order to only keep experiments with a robust regions repartition for most active neurons. 

We did that using a threshold of 0.5 for the p-value. This means that for each dataset related to a particular stimulus, we made a list of all the p-values from Kolmogorov-Smirnov tests with all the other datasets related to the same stimulus. The we computed the averaged p-value. If this average was lower than 0.5, we got rid of the dataset, otherwise we kept it. This resulted in keeping datasets only if they looked similar in terms of regions to datasets related to the same stimulus. Using this technique, we went from 55 datasets (all stimuli considered) to 46 datasets. This allowed some interesting cleaning of the data.

![Distribution of p-values across stimuli after removing datasets](/README_img/boxplot2.png)

### What correlation to use
