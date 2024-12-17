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
        # prepare input file
        input_file = '006.wf-p0.dat'
        output_file = '006.wf-p0.gbs'
        shutil.copy(os.path.join('test_files',input_file), 'tmp/')
        os.chdir('tmp')
        # main
        from thunder_ase.utils import fit_gaussian
        fit_gaussian(input_name=input_file, plot=True)
        assert os.path.isfile(output_file)