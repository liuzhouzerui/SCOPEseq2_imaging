import pickle

from scopeseq.Bead_Intensity import BeadIntensity
from scopeseq.Well_Bead_Cell import WellBeadCell
from scopeseq.Cell_Features import CellFeatures

# bead intensity
project_folder = "/Users/liuzhouzerui/Desktop/sims/scope_aws/GBM/"
image_folder = project_folder + "image/"
bead_folder = project_folder + "image/bead/"
channels = ["c1", "c2", "c3"]
lanes = 1
patches = 10
first_border = 600
border_gap = 7400
obc_ref_fn = project_folder + "image/192_8mer_seq_reference.csv"
for i in range(lanes):
    print('doing lane ' + str(i) + '...')
    # create bead-intensity object for each single lane
    bead_intensity = BeadIntensity(bead_folder, channels, n_lane=i, total_lanes=lanes, total_patches_per_lane=patches)
    bead_intensity.generate_bead_intensity()
    bead_intensity.assign_patch(first_border=first_border, border_gap=border_gap)
    bead_intensity.obc_calling(obc_ref_fn)
    # write bead intensity
    f = open(image_folder + 'bead_intensity_' + str(i) + '.obj', 'wb')
    pickle.dump(bead_intensity, f)
    f.close()

# well-bead-cell link
well_folder = project_folder + "image/cell/"
cell_folder = project_folder + "image/cell/"
well_channel = 'c2'
cell_image_channels = ['c2', 'c3']
for i in range(lanes):
    print('doing lane ' + str(i) + '...')
    well_bead_cell = WellBeadCell(well_folder, n_lane=i)
    well_bead_cell.initialize_well(well_channel)
    # load bead obc calling result
    f = open(image_folder + 'bead_intensity_' + str(i) + '.obj', 'rb')
    bead_intensity = pickle.load(f)
    f.close()
    well_bead_cell.link_bead(bead_intensity)
    # load cell features
    well_bead_cell.link_cell(cell_folder, cell_image_channels)
    well_bead_cell.link_obc_cell()
    # write well-bead-cell link result to file
    f = open(image_folder + 'well_bead_cell_' + str(i) + '.obj', 'wb')
    pickle.dump(well_bead_cell, f)
    f.close()

# image-sequencing link
# link cell_ID, cell_feature
cell_id_fn = project_folder + "PJ069.cell_ID.txt"
cell_features = CellFeatures(cell_id_fn=cell_id_fn)
for i in range(lanes):
    print('doing lane ' + str(i) + '...')
    f = open(image_folder + 'well_bead_cell_' + str(i) + '.obj', 'rb')
    well_bead_cell = pickle.load(f)
    f.close()
    cell_features.seq_image_link(well_bead_cell)
# the cell feature relationship to cluster
diffex_folder = project_folder + "diffex/"
cell_features.load_cluster(diffex_folder)
f = open(image_folder + 'cell_features.obj', 'wb')
pickle.dump(cell_features, f)
f.close()

data = cell_features.to_dataframe()

