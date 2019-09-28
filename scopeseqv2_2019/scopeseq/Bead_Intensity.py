import os
import fnmatch
import pandas as pd
import numpy as np

from scopeseq.method import register, assign_obc


class BeadIntensity:
    def __init__(self, bead_intensity_folder, channel_name, n_lane=0, total_lanes=1, total_patches_per_lane=10, round_number=8, d_th=40):
        """

        :param bead_intensity_folder:
        :param channel_name: should list as [bf, CY5_Sprobe, CY3_Qprobe]
        :param n_lane:
        :param total_lanes:
        :param total_patches_per_lane:
        :param round_number: number of probe hybridization round
        :param d_th: the distance threshold of well-bead, well-cell
        """
        self.bead_intensity_folder = bead_intensity_folder
        self.channel_name = channel_name
        self.n_lane = n_lane
        self.total_lanes = total_lanes
        self.n_patch_per_lane = total_patches_per_lane
        self.bead_area = None
        self.bead_position = None
        self.probe = None
        self.bg = None
        self.borders = None
        self.patch = np.array([])
        self.obc = []
        self.obc_s = []
        self.obc_q = []
        self.round = round_number
        self.d_th = d_th

    def find_bead_intensity_file(self, iter_round, channel, bp):
        """

        :param iter_round:
        :param channel:
        :param bp:
        :return:
        """
        file_round = '0000' + str(iter_round - 1)
        if bp == 'bg':
            file_seq = str(10000 + (iter_round - 1) * self.total_lanes * 2 + self.n_lane)[1:5]
        elif bp == 'probe':
            file_seq = str(10000 + (iter_round - 1) * self.total_lanes * 2 + self.total_lanes + self.n_lane)[1:5]
        else:
            print("no such file!")
            file_seq = ''
        file_ext = self.channel_name[channel] + ".tif_Results.xls"
        file_pattern = '*' + file_round + '*' + file_seq + '*' + file_ext
        file_name = fnmatch.filter(os.listdir(self.bead_intensity_folder), file_pattern)[0]
        return file_name

    def initialize(self, iter_round, channel, bp):
        """

        :param iter_round:
        :param channel:
        :param bp:
        :return:
        """
        print("initializing...")
        file_name = self.find_bead_intensity_file(iter_round, channel, bp)
        print(file_name)
        intensity_matrix_oneround = pd.read_csv(self.bead_intensity_folder + file_name, sep='\t')
        total_bead_num = intensity_matrix_oneround.shape[0]
        self.bead_area = intensity_matrix_oneround['Area']
        self.bead_position = intensity_matrix_oneround[['XM', 'YM']]
        probe_header = []
        for c in self.channel_name[1:3]:
            for i in range(1, self.round+1, 1):
                probe_header.append('probe_' + str(i) + '_' + c)
        self.probe = pd.DataFrame(-1, columns=probe_header, index=range(total_bead_num))
        bg_header = []
        for c in self.channel_name[1:3]:
            for i in range(1, self.round+1, 1):
                bg_header.append('bg_' + str(i) + '_' + c)
        self.bg = pd.DataFrame(-1, columns=bg_header, index=range(total_bead_num))
        if bp == 'bg':
            self.bg[bp + '_' + str(iter_round) + '_' + self.channel_name[channel]] = intensity_matrix_oneround[
                'Mean'].values
        if bp == 'probe':
            self.probe[bp + '_' + str(iter_round) + '_' + self.channel_name[channel]] = intensity_matrix_oneround[
                'Mean'].values
        del intensity_matrix_oneround

    def register_bead(self, iter_round, channel, bp):
        """

        :param iter_round:
        :param channel:
        :param bp:
        :return:
        """
        print("registering..." + bp + '_' + str(iter_round) + '_' + self.channel_name[channel])
        file_name = self.find_bead_intensity_file(iter_round, channel, bp)
        print(file_name)
        intensity_matrix_oneround = pd.read_csv(self.bead_intensity_folder + file_name, sep='\t')
        d_position, d_inrange = register(intensity_matrix_oneround[['XM', 'YM']], self.bead_position, self.d_th)
        if bp == 'bg':
            self.bg.loc[d_position[d_inrange], bp + '_' + str(iter_round) + '_' + self.channel_name[channel]] \
                = intensity_matrix_oneround.iloc[np.arange(d_position.size)[d_inrange], 2].values
        if bp == 'probe':
            self.probe.loc[d_position[d_inrange], bp + '_' + str(iter_round) + '_' + self.channel_name[channel]] \
                = intensity_matrix_oneround.iloc[np.arange(d_position.size)[d_inrange], 2].values
        del intensity_matrix_oneround

    def generate_bead_intensity(self):
        """
        de-multiplexing image intensity.
        :return:
        """
        for iter_round in range(1, self.round+1, 1):
            for channel in [1, 2]:
                for bp in ['bg', 'probe']:
                    if self.bead_position is None:
                        self.initialize(iter_round, channel, bp)
                    else:
                        self.register_bead(iter_round, channel, bp)

    def assign_patch(self, first_border=600, border_gap=7400):
        """

        :param first_border:
        :param border_gap:
        :return:
        """
        print('assigning patches...')
        self.borders = np.arange(first_border, first_border + (self.n_patch_per_lane + 1) * border_gap, border_gap)
        self.patch = self.bead_position['XM'].apply(lambda x: np.argmax((self.borders > x) * 1)-1)
        self.patch = self.patch + self.n_patch_per_lane * self.n_lane

    def obc_calling(self, barcode_ref_fn, no_signal_th=None, mode='all'):
        """

        :param barcode_ref_fn: a file contains Barcode_S and Barcode_Q
        :param no_signal_th:
        :param mode:
        :return:
        """
        print('optical barcode calling...')
        barcode_ref_table = pd.read_csv(barcode_ref_fn, dtype=str)
        self.obc_s = self.probe.iloc[:, 0:self.round].apply(
            lambda x: assign_obc(x, barcode_ref_table['Barcode_S'], no_signal_th=no_signal_th, mode=mode), axis=1).values
        self.obc_q = self.probe.iloc[:, self.round:self.probe.shape[1]].apply(lambda x: assign_obc(x, barcode_ref_table[
            'Barcode_Q'], no_signal_th=no_signal_th, mode=mode), axis=1).values
        self.obc = self.patch.astype('str') + '_' + self.obc_s.astype('str') + '_' + self.obc_q.astype('str')
