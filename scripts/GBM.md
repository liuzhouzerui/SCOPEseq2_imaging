# Example workflow 

#### import packages

```
import pickle
from scopeseq.Bead_Intensity import BeadIntensity
from scopeseq.Well_Bead_Cell import WellBeadCell
from scopeseq.Cell_Features import CellFeatures
```
### image analysis

#### parameters - demultiplexing
bead_folder - a bead folder that contains de-multiplexing fluorescent measurement. \
Each file-name should contain: \
 (1) 00000,00001,00002... as de-multiplexing round indicator; \
 (2) 0000, 0001, 0002... as image number indicator; \
 (3) c1,c2,c3 or BF,CY5,CY3 or user_defined other names as image channel indicator; \
 (4) .tif_Results.xls as file extension.
 
 ```
project_folder = "~/GBM/"
bead_folder = project_folder + "image/bead/"
```
 channels - a list of image channel indicator. \
 The first element is the name of brightfield; \
 The second element is the name of S-probe associated CY5; \
 The third element is the name of Q-probe associated CY3; 
 
 ```
channels = ["BF", "CY5", "CY3"]
```
 or
 ```
channels = ["c1", "c2", "c3"]
```

 lanes - the total number of lanes of the device. \
 patches - the number of patches of each lane used for sequencing. \
 first_border - the x coordinate (pixel) of the first border of the patches. \
 border_gap - the pixels between the two borders of the patch. \
 round_number - the number of rounds used for de-multiplexing, default 8. \
 d_th_bead - the distance threshold for registering bead, default 40. if using low resolution for 
 de-multiplexing image, then it's 20. \
 d_th_cell - the distance threshold for registering bead, default is d_th_bead.

```
lanes = 1
patches = 10
first_border = 600
border_gap = 7400
round_number = 8
d_th_bead = 40
```
obc_ref_fn - the obc barcode reference filename.

```
obc_ref_fn = project_folder + "image/192_8mer_seq_reference.csv"
```

mode - 'all' for iteration, 'max' for only max intensity change. 
```
mode = 'all'
```

#### parameters - register bead and cell information to wells
well_folder - a well folder that contains single well measurement, and the landmark coordinates file. \
Each single well measurement file-name should contain: \
(1) 0000, 0001, 0002... as image number (lane) indicator; \
(2) channel name for the well measurement; \
(3) .tif_Results.xls as file extension. \
note 1: only well-location from the cell scan is used in the analysis for registering. 
If the cell features are well-based measurement, it can also be used as a well measurement. \
note 2: only 1 measurement (1 channel and 1 round image) for each lane is needed. 

The landmark coordinates file: \
the file-name should end with "UL_BR_{lane#}.csv". \
the file should have measurements from imageJ, with following orders: 
cell.scan_UL_well, cell.scan_BR_well, bead.scan_UL_well, bead.scan_BR_well.

```
well_folder = project_folder + "image/well/"
well_channel = 'c2'
```

cell_folder - a cell folder that contains single cell measurement. \
Each file-name should contain: \
(1) 0000, 0001, 0002... as image number (lane) indicator; \
(2) a list as the channel name for the cell measurement; \
(3) .tif_Results.xls as file extension. \
note: if the cell measurement is well based, 
then you need to manually filter out those empty wells in the `WellBeadCell.obc_cell` before sequencing linkage.

```
cell_folder = project_folder + "image/cell/"
cell_image_channels = ['c2', 'c3']
```

d_th_cell - the distance threshold, used in linking bead and cell to wells.
```
d_th_cell = 40
```


#### objects - bead intensity

The following code would generate a BeadIntensity object for each single lane, \
and store as an object using pickle. \
The BeadIntensity object contains: \
(1) self.bead_intensity_folder - the bead folder \
(2) self.channel_name - the channels of the image \
(3) self.n_lane - the number of lane, which the object related to \
(4) self.total_lanes - the total lanes of the device \
(5) self.n_patch_per_lane - the number of patches of each lane used for sequencing. \
(6) self.round - the number of rounds used for de-multiplexing \
(7) self.d_th - the distance threshold for register \
(8) self.bead_area - the area of the bead \
(9) self.bead_position - the location of the bead (pixel) \
(10) self.probe - the probe intensity of the beads \
(11) self.bg - the background intensity of the beads \
(12) self.patch - the patch-id of the bead \
(13) self.obc - the optical barcode of the bead \
e.g. 0_19_23 means the 0th patch, 19th S-probe, 23th Q-probe \

```
for i in range(lanes):
    print('doing lane ' + str(i) + '...')
    # create bead-intensity object for each single lane
    bead_intensity = BeadIntensity(bead_folder, channels, n_lane=i, total_lanes=lanes, total_patches_per_lane=patches
                                   d_th=d_th_bead)
    bead_intensity.generate_bead_intensity()
    bead_intensity.assign_patch(first_border=first_border, border_gap=border_gap)
    bead_intensity.obc_calling(obc_ref_fn, no_signal_th=None, mode=mode)
    # write bead intensity
    f = open(image_folder + 'bead_intensity_' + str(i) + '.obj', 'wb')
    pickle.dump(bead_intensity, f)
    f.close()
```

#### objects - well_bead_cell link

The following code would generate a WellBeadCell object for each single lane, \
and store as an object using pickle. It uses well location as a reference, links the bead and the cell to their closest well. \
The WellBeadCell object contains: \
(1) self.n_lane - the number of lane, which the object related to \
(2) self.well_folder - the well folder \
(3) self.well_position - the well locations, as a reference for linking beads and cells \
(4) self.bead - the lane associated BeadIntensity object \
(5) self.rotation_matrix - the rotation matrix for aligning the bead and cell \
(6) self.bead_rotated_position - the rotated bead position \
(7) self.well_bead_link - bead id linked to the well \
(8) self.cell_folder - the cell folder \
(9) self.cell_position - the cell position \
(10) self.cell - the cell features, stored as a dict. key is the channel name, value is the measurement dataframe \
(11) self.well_cell_link - cell id linked to the well \
(12) self.obc_cell - the bead optical barcode linked to the cell id

```
for i in range(lanes):
    print('doing lane ' + str(i) + '...')
    well_bead_cell = WellBeadCell(well_folder, n_lane=i)
    well_bead_cell.initialize_well(well_channel)
    # load bead obc calling result
    f = open(image_folder + 'bead_intensity_' + str(i) + '.obj', 'rb')
    bead_intensity = pickle.load(f)
    f.close()
    well_bead_cell.link_bead(bead_intensity, d_th=d_th_cell)
    # load cell features
    well_bead_cell.link_cell(cell_folder, cell_image_channels)
    well_bead_cell.link_obc_cell()
    # write well-bead-cell link result to file
    f = open(image_folder + 'well_bead_cell_' + str(i) + '.obj', 'wb')
    pickle.dump(well_bead_cell, f)
    f.close()

```
Optional: The following code is an example of well-based cell measurement. (in mixed species experiment)\
 Manually filtering out those empty wells in the `WellBeadCell.obc_cell` before sequencing linkage. 
 
```
plt.hist(np.log2(well_bead_cell.cell['c2']['c2_Mean']))
plt.show()
plt.hist(np.log2(well_bead_cell.cell['c3']['c3_Mean']))
plt.show()
th_c2 = 2**9.5  # mannually chosen threshold for calling cells in the c2 (calcein green human) channel
th_c3 = 2**11  # mannually chosen threshold for calling cells in the c3 (calcein red mouse) channel

obc_cell = well_bead_cell.obc_cell
obc_cell['c2_mean'] = well_bead_cell.cell['c2'].loc[obc_cell['cell_num'], 'Mean'].values
obc_cell['c3_mean'] = well_bead_cell.cell['c3'].loc[obc_cell['cell_num'], 'Mean'].values
obc_cell = obc_cell[(obc_cell['c2_mean'] > th_c2) | (obc_cell['c3_mean'] > th_c3)]

well_bead_cell.obc_cell = pd.DataFrame(obc_cell[['obc', 'cell_num']].values, columns=['obc', 'cell_num'], 
                                       index=range(obc_cell.shape[0]))
```


### linking image optical barcode to sequencing cell ids.

#### parameters

cell_id - a list or an array of the cell ids. \
if it's None, then it will use cell_id_fn to generate cell ids. 

cell_id_fn - the filename of the cell id file. \
The cell id file should contain: one row of cell ids, white-delimiter separated. \
note: this file could be generated by `head -n 1 count-matrix.txt > cell_id.txt`. It's ok to contain `gid`, `gene` 
and `\n`, as they will be removed.

```
cell_id_fn = project_folder + "seq/PJ069.cell_ID.txt"
```
optional: diffex_folder - the output folder from the diffexcluster analysis. \
the diffex _folder should contain: \
(1) a cluster file, end with ".pg.txt" \
(2) a umap coordinate file, end with ".umap.txt"

```
diffex_folder = project_folder + "diffex/"
```

#### objects
The following code would generate a CellFeatures object for the whole device, \
and store as an object using pickle. \
The CellFeatures object contains: \
(1) self.cell_id - the cell ids \
(2) self.seqobc_cell - the sequencing cell id and the linked cell number id \
(3) self.cluster - the cluster of the cell \
(4) self.umap - the umap coordinate of the cell \
(5) self.count - the count matrix of the cell, cell id (row) by gene (col) 

```
cell_features = CellFeatures(cell_id=None, cell_id_fn=cell_id_fn)
for i in range(lanes):
    print('doing lane ' + str(i) + '...')
    f = open(image_folder + 'well_bead_cell_' + str(i) + '.obj', 'rb')
    well_bead_cell = pickle.load(f)
    f.close()
    cell_features.seq_image_link(well_bead_cell)
```
```
# link the cell feature relationship to cluster
cell_features.load_cluster(diffex_folder)
```

```
f = open(image_folder + 'cell_features.obj', 'wb')
pickle.dump(cell_features, f)
f.close()
```

### for down-stream analysis
this will create a dataframe. cell-id as the row.names, 
the cluster (if apply), umap coordiantes (if apply), expression (if apply) and the cell features as col.names.

```
data = cell_features.to_dataframe()
```
