import numpy as np
import pytest
import shutil
import os

class TestFitGaussian:
    def setup_method(self):
        self.cwd = os.getcwd()
        # create tmp directory and get in
        if os.path.isdir('tmp'):
            shutil.rmtree('tmp')
        os.mkdir('tmp')

    def teardown_method(self):
        # clean all result files
        os.chdir(self.cwd)
        shutil.rmtree('tmp')

    def test_fit_gaussian_command(self):
        input_file = '006.wf-p0.dat'
        output_file = '006.wf-p0.gbs'
        shutil.copy(os.path.join('test_files',input_file), 'tmp/')
        os.chdir('tmp')
        os.system('fit-gaussians -p '+input_file)
        assert os.path.isfile(output_file)

    def test_fit_gaussian(self):
        # prepare an input file
        input_file = '006.wf-p0.dat'
        output_file = '006.wf-p0.gbs'
        shutil.copy(os.path.join('test_files',input_file), 'tmp/')
        os.chdir('tmp')
        # main
        from thunder_ase.utils import fit_gaussian
        fit_gaussian(input_name=input_file, plot=True)
        assert os.path.isfile(output_file)

    def test_fit_gaussian2(self):
        # prepare an input file
        input_file = '006.wf-s0.dat'
        output_file = '006.wf-s0.gbs'
        shutil.copy(os.path.join('test_files',input_file), 'tmp/')
        os.chdir('tmp')
        # main
        from thunder_ase.utils import fit_gaussian
        fit_gaussian(input_name=input_file, plot=True)
        assert os.path.isfile(output_file)

    def test_orthogonality(self):
        from thunder_ase.utils.basis_set import get_full_wf, read_wf
        input_file1 = 'test_files/006.wf-p0.dat'
        input_file2 = 'test_files/006.wf-s0.dat'
        wf1 = get_full_wf(input_file1)
        wf2 = get_full_wf(input_file2)
        #wf1 = read_wf(input_file1).T
        #wf2 = read_wf(input_file2).T
        # 对齐 wf1 和 wf2
        assert wf1[0][0] - wf1[0][1] == wf2[0][0] - wf2[0][1]  # x should have same interval
        length_max = max(wf1.shape[1], wf2.shape[1])
        new_wf1 = np.zeros(length_max)
        new_wf2 = np.zeros(length_max)
        new_wf1[:wf1.shape[1]] = wf1[1]
        new_wf2[:wf2.shape[1]] = wf2[1]

        orthogonality = np.dot(new_wf1, new_wf2)
        print(orthogonality)