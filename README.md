# multiSensorint
Program that allows to compare multi-sensory brain scans, through the use of a grid. Comparison is made possible with a graphical interface. Loading a ZBG object allows then to chose subsets (e.g. visual or vestibular stimuli) and to compare their most active neurons. 


## The ZBG class
@ZBraingrid is a class that aims to facilitating manipulation of grids associated to datasets. 

### Defining a ZBG object
A ZBG object has three intrinsic values: 

* a *method* [character|string], defining name of object,

* a *size* [array of integers, length 1 or 3], which is a vector of the number of voxels in the grid along all 3 dimensions. If one wants to define grid based on voxel size, one can use grid dimensions to compute final grid size. Grid dimensions are ([0, 0.496], [0, 1.122], [0, 0.276]), reference brain from the [Z Brain Atlas](https://engertlab.fas.harvard.edu/Z-Brain/home/),

* an *orientation* [character], indicating where axes (x, y, z) point to. Orientation should be a character (e.g. 'RAS' or 'LPI'), with R:right, L:left, A:anterior, P:posterior, S:superior, I:inferior.

### Adding a dataset to a ZBG object
The function addDataset allows to add a dataset to an already existing object. It takes two mandatory arguments and one optional argument. First argument is the ZBG object, second argument is the dataset structure (explained below), and third argument is to disable information printing while new data set is added (just type 'information_off' in this case). 

Dataset structure is a structure containing mandatory information on new dataset to add. As grid size is set when ZBG object is created, new data will already adjust this size. Following are the fields required:

* *name* [character|string]: defining name of new dataset,

* *path* [character|string]: defining path to actual data (HDF5 file for us),

* *comment* [character|string]: defining key words associated to dataset. Comment is important because it is used when discriminating ZBG subset from ZBG object, so chose wisely,

* *orientation* [character]: see orientation from section above. Orientation is automatically adjusted to be the same as the one from the ZBG object,

* *coordinates* [matrix of float|double, size (n x 3)]: matrix with as many rows as there are neurons to add (n), and 3 columns for the three coordinates along x, y, and z axes,

* *correlation* [array of float|double, length n]: vector with a length equal to the number of neurons from coordinates (n). Correlation does not have to be actual correlation, but rather a number linking neurons' signals to stimulus. We discuss this aspect in a following section. If user only wants to compute projection on grid, then they should enter random values.

When adding a new dataset to ZBG object, addDataset computes for each neuron the voxel to which it belongs. To each neuron is therefore attributed an index in the grid. Then, correlations are averaged for each voxel. It is this averaged correlation value that is later used in the graphical interface. **When we speak about correlation later, we refer to this averaged correlation**.

### ZBG methods 

We will not make an extensive list for all ZBG methods here, but rather describe the most important and useful ones. User can visualize all methods in the ZBG class Matlab file. We already went through addDataset, which allows to add new datasets to an already existing object. Basic operations are provided, such as addition (merge two ZBG objects), length (returns number of datasets in ZBG object) and indexing (get new ZBG object with datasets corresponding to indexes). Following are some of the other methods:

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

### What correlation to use
