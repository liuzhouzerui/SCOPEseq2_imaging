import os
import fnmatch
import pandas as pd
import numpy as np

from scopeseq.method import find_unique_match_position


class CellFeatures:
    def __init__(self, cell_id=None, cell_id_fn=None):
        """

        :param cell_id: a list or an array of cell_id from the sequencing run
        :param cell_id_fn: the cell_id filename. the file should contain one row of cell_id, '\t' separated.
                            ['gid','gene','\n'] will be removed.
        """
        if cell_id is None:
            f = open(cell_id_fn, 'r')
            cell_id = f.read().split('\t')
            f.close()
            for x in ['gid', 'gene', '\n']:
                cell_id.remove(x)
        self.cell_id = cell_id
        self.seqobc_cell = pd.DataFrame(index=self.cell_id)
        self.cluster = None
        self.umap = None
        self.count = None

    def seq_image_link(self, well_bead_cell):
        """

        :param well_bead_cell: a WellBeadCell object.
        :return:
        """
        print("linking sequencing barcode to cell# ...")
        obc_position = np.array([find_unique_match_position(x, well_bead_cell.obc_cell['obc']) for x in self.cell_id])
        cell_num = well_bead_cell.obc_cell.iloc[obc_position[obc_position > 0], :]
        self.seqobc_cell.loc[cell_num['obc'], 'lane'] = well_bead_cell.bead.n_lane
        self.seqobc_cell.loc[cell_num['obc'], 'cell_num'] = cell_num['cell_num'].values
        for key in well_bead_cell.cell.keys():
            for column in well_bead_cell.cell[key].columns.values:
                self.seqobc_cell.loc[cell_num['obc'], key + '_' + column] = \
                    well_bead_cell.cell[key].loc[cell_num['cell_num'], column].values

    def load_cluster(self, diffex_folder):
        """

        :param diffex_folder:
        :return:
        """
        cluster_fn = fnmatch.filter(os.listdir(diffex_folder), "*.pg.txt")[0]
        self.cluster = pd.read_csv(diffex_folder + cluster_fn, header=None, names=['cluster'])
        umap_fn = fnmatch.filter(os.listdir(diffex_folder), "*.umap.txt")[0]
        self.umap = pd.read_csv(diffex_folder + umap_fn, header=None, names=['X', 'Y'], sep='\t')

    def load_counts(self, count):
        """
        :param count: count matrix, row.names are cell_id, col.names are gene_id
        :return:
        """
        self.count = count

    def to_dataframe(self):
        """

        :return:
        """
        # channels is a list
        data = pd.DataFrame(index=self.cell_id)
        if self.cluster is not None:
            data['cluster'] = self.cluster.values
        if self.umap is not None:
            data['umap_x'] = self.umap['X'].values
            data['umap_y'] = self.umap['Y'].values
        if self.count is not None:
            data[self.count.columns] = self.count
        if not self.seqobc_cell.empty:
            data[self.seqobc_cell.columns] = self.seqobc_cell
        return data
