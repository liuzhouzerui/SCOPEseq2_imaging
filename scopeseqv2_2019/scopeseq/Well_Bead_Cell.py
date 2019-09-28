import os
import fnmatch
import pandas as pd
import numpy as np

from scopeseq.method import register, find_rotation_matrix, find_unique_match_position


class WellBeadCell:
    def __init__(self, well_folder, n_lane=0):
        """

        :param well_folder: a well folder should contain the UL_BR_{n_lane}.csv
        :param n_lane:
        """
        self.n_lane = n_lane
        self.well_folder = well_folder
        self.well_position = None
        self.bead = None
        self.bead_rotated_position = None
        self.rotation_matrix = None
        self.well_bead_link = None
        self.cell = {}
        self.cell_folder = None
        self.cell_position = None
        self.well_cell_link = None
        self.obc_cell = None

    def initialize_well(self, channel):
        """

        :param channel: channel name used for the well image
        :return:
        """
        file_seq = str(10000 + self.n_lane)[1:5]
        file_ext = channel + ".tif_Results.xls"
        file_pattern = '*' + file_seq + '*' + file_ext
        file_name = fnmatch.filter(os.listdir(self.well_folder), file_pattern)[0]
        print("initializing... " + file_name)
        well = pd.read_csv(self.well_folder + file_name, sep="\t")
        self.well_position = well[['XM', 'YM']]
        self.well_bead_link = np.repeat(-1, well.shape[0])
        self.well_cell_link = np.repeat(-1, well.shape[0])

    def rotate_bead_position(self):
        """
        rotate de-multiplexing image to match the well-cell scan image
        :return:
        """
        # rotate bead position to match well position
        # the well position is the same with the cell position
        self.bead_rotated_position = pd.DataFrame(-1, columns=['XM', 'YM'], index=range(self.bead.bead_position.shape[0]))
        landmark_fn = fnmatch.filter(os.listdir(self.well_folder), '*UL_BR_' + str(self.n_lane) + '.csv')[0]
        landmark_all = pd.read_csv(self.well_folder + landmark_fn)
        # rotate for each small parts, parts number = landmark.shape[0]/4
        parts_number = int(landmark_all.shape[0]/4)
        if parts_number == 1:
            landmark = landmark_all
            target_vector = np.array(
                [landmark.iloc[1, 9] - landmark.iloc[0, 9], landmark.iloc[1, 10] - landmark.iloc[0, 10]])
            initial_vector = np.array(
                [landmark.iloc[3, 9] - landmark.iloc[2, 9], landmark.iloc[3, 10] - landmark.iloc[2, 10]])
            self.rotation_matrix = find_rotation_matrix(target_vector, initial_vector)
            self.bead_rotated_position = np.array(landmark.iloc[0, 9:11]) + \
                np.dot(np.array(self.bead.bead_position - landmark.iloc[2, 9:11]), self.rotation_matrix) * \
                target_vector / np.dot(initial_vector, self.rotation_matrix)
        else:
            borders = np.arange(self.bead.borders[0], self.bead.borders[-1]+(self.bead.borders[-1]-self.bead.borders[0])
                                / parts_number, (self.bead.borders[-1]-self.bead.borders[0])/parts_number)
            bead_parts = self.bead.bead_position['XM'].apply(lambda x: np.argmax((borders > x) * 1) - 1)
            for i in range(parts_number):
                landmark = landmark_all.iloc[(i*4):((i+1)*4), :]
                target_vector = np.array(
                    [landmark.iloc[1, 9] - landmark.iloc[0, 9], landmark.iloc[1, 10] - landmark.iloc[0, 10]])
                initial_vector = np.array(
                    [landmark.iloc[3, 9] - landmark.iloc[2, 9], landmark.iloc[3, 10] - landmark.iloc[2, 10]])
                self.rotation_matrix = find_rotation_matrix(target_vector, initial_vector)
                self.bead_rotated_position[bead_parts == i] = np.array(landmark.iloc[0, 9:11]) + \
                    np.dot(np.array(self.bead.bead_position[bead_parts == i] - landmark.iloc[2, 9:11]), self.rotation_matrix) * \
                    target_vector / np.dot(initial_vector, self.rotation_matrix)

    def link_bead(self, bead, d_th=None):
        """

        :param bead: BeadIntensity object
        :param d_th:
        :return:
        """
        print("linking bead information to wells...")
        if d_th is None:
            d_th = self.bead.d_th
        self.bead = bead
        self.rotate_bead_position()
        d_position, d_inrange = register(self.bead_rotated_position, self.well_position, d_th)
        self.well_bead_link[d_position[d_inrange]] = np.arange(d_position.size)[d_inrange]

    def link_cell(self, cell_folder, channels):
        """

        :param cell_folder:
        :param channels: channel names for the cell images, e.g. ['GFP', 'TRITC']
        :return:
        """
        print("linking cell information to wells...")
        # channel is a list
        self.cell_folder = cell_folder
        for channel in channels:
            file_seq = str(10000 + self.n_lane)[1:5]
            file_ext = channel + ".tif_Results.xls"
            file_pattern = '*' + file_seq + '*' + file_ext
            file_name = fnmatch.filter(os.listdir(self.cell_folder), file_pattern)[0]
            cell_property = pd.read_csv(self.cell_folder + file_name, sep="\t")
            self.cell[channel] = cell_property
        self.cell_position = cell_property[['XM', 'YM']]
        del cell_property
        # register cell to wells, w. doublet removal
        d_position, d_inrange = register(self.well_position, self.cell_position, self.bead.d_th, doublet_removal=True)
        self.well_cell_link[np.arange(d_position.size)[d_inrange]] = d_position[d_inrange]

    def link_obc_cell(self):
        """

        :return:
        """
        print("linking cell to bead...")
        position = np.array([find_unique_match_position(x, self.well_bead_link) for x in range(self.bead.obc.size)])
        self.obc_cell = pd.DataFrame(-1, columns=['obc', 'cell_num'], index=range(sum(position > 0)))
        self.obc_cell['obc'] = self.bead.obc[np.arange(position.size)[position > 0]].values
        self.obc_cell['cell_num'] = self.well_cell_link[position[position > 0]]
        self.obc_cell = self.obc_cell.iloc[np.where(self.obc_cell['cell_num'] > 0)]
        self.obc_cell.index = np.arange(self.obc_cell.shape[0])
