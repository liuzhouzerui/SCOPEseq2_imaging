import os
import fnmatch
import pandas as pd
import numpy as np

from scopeseq.method import register, assign_obc, assign_obc_stage, quantile_linear_transformation


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
        self.probe_normalized = None
        self.bg = None
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

    def intensity_normalization(self, method='ql'):
        """

        :param method: 'ql', quantile linear transformation
        :return:
        """
        if method == 'ql':
            round = [1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8]
            self.probe_normalized = pd.DataFrame(0, columns=self.probe.columns, index=self.probe.index)
            for i in range(len(round)):
                self.probe_normalized.iloc[:, i], m1, m2 = quantile_linear_transformation(self.probe.iloc[:, i], round[i])
                print(self.probe.columns[i], m1, m2)

    def assign_patch(self, first_border=600, border_gap=7400):
        """

        :param first_border:
        :param border_gap:
        :return:
        """
        print('assigning patches...')
        borders = np.arange(first_border, first_border + (self.n_patch_per_lane + 1) * border_gap, border_gap)
        self.patch = self.bead_position['XM'].apply(lambda x: np.argmax((borders > x) * 1)-1)
        self.patch = self.patch + self.n_patch_per_lane * self.n_lane

    def obc_calling(self, barcode_ref_fn, no_signal_th=None, mode='all'):
        """
        :param barcode_ref_fn: a file contains Barcode_S and Barcode_Q
        :param no_signal_th:
        :param mode: 'max' is one round core-by-core method, 'all' is full round core-by-core method
        :return:
        """
        print('optical barcode calling using core-by-core method...')
        barcode_ref_table = pd.read_csv(barcode_ref_fn, dtype=str)
        if self.probe_normalized is None:
            self.obc_s = self.probe.iloc[:, 0:self.round].apply(
                lambda x: assign_obc(x, barcode_ref_table['Barcode_S'], no_signal_th=no_signal_th, mode=mode),
                axis=1, result_type="expand")
            self.obc_q = self.probe.iloc[:, self.round:self.probe.shape[1]].apply(
                lambda x: assign_obc(x, barcode_ref_table[
                    'Barcode_Q'], no_signal_th=no_signal_th, mode=mode), axis=1, result_type="expand")
        else:
            self.obc_s = self.probe_normalized.iloc[:, 0:self.round].apply(
                lambda x: assign_obc(x, barcode_ref_table['Barcode_S'], no_signal_th=no_signal_th, mode=mode),
                axis=1, result_type="expand")
            self.obc_q = self.probe_normalized.iloc[:, self.round:self.probe_normalized.shape[1]].apply(
                lambda x: assign_obc(x, barcode_ref_table[
                    'Barcode_Q'], no_signal_th=no_signal_th, mode=mode), axis=1, result_type="expand")
        self.obc = self.patch.astype('str') + '_' + self.obc_s[0].astype('str') + '_' + self.obc_q[0].astype('str')
        self.obc_round = self.obc_s[1].astype('str') + '_' + self.obc_q[1].astype('str')


    def obc_calling_stage(self, barcode_ref_fn):
        """
        barcode calling using stage-by-stage method
        :param barcode_ref_fn:
        :return:
        """
        print('optical barcode calling using stage-by-stage method...')
        barcode_ref_table = pd.read_csv(barcode_ref_fn, dtype=str)
        obc = pd.DataFrame(None, index=self.probe.index, columns=self.probe.columns)
        # determine cut-off of each stage and assign on-off
        for i in self.probe.columns:
            chist, cedge = np.histogram(np.log2(self.probe[self.probe[i] > 0][i]), bins=50)
            md = np.median(cedge)
            mn1 = cedge[np.argmax(chist[cedge[0:-1] < md])]
            mn2 = cedge[np.argmax(chist[cedge[0:-1] > md]) + len(cedge[cedge <= md])]
            # the shortest bin between two mean as thresh
            thresh = cedge[np.argmin(chist[(cedge[0:-1] > mn1) & (cedge[0:-1] < mn2)]) + len(cedge[cedge <= mn1])]
            obc[i][(np.log2(self.probe[i]) > 0) & (np.log2(self.probe[i]) <= thresh)] = 0
            obc[i][np.log2(self.probe[i]) > thresh] = 1
            obc[i][self.probe[i] == -1] = -1
        self.obc_s_stage = obc.iloc[:, 0:self.round].apply(
            lambda x: assign_obc_stage(x, barcode_ref_table['Barcode_S']), axis=1).values
        self.obc_q_stage = obc.iloc[:, self.round:self.probe.shape[1]].apply(
            lambda x: assign_obc_stage(x, barcode_ref_table['Barcode_Q']), axis=1).values
        self.obc_stage = self.patch.astype('str') + '_' + self.obc_s_stage.astype(
            'str') + '_' + self.obc_q_stage.astype('str')

