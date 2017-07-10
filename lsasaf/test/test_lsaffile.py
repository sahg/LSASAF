import numpy as np

import lsasaf

def test_read_raw_dataset_compressed():
    data_dir = 'data/sample_data_files/'

    dssf_fname = data_dir + 'S-LSA_-HDF5_LSASAF_MSG_DSSF_SAfr_201706091400'

    dssf_file = lsasaf.DSSFFile(dssf_fname)
    decomp_dset = dssf_file.read_raw_dataset('DSSF')

    dssf_fname = data_dir + 'S-LSA_-HDF5_LSASAF_MSG_DSSF_SAfr_201706091400.bz2'

    dssf_file = lsasaf.DSSFFile(dssf_fname)
    comp_dset = dssf_file.read_raw_dataset('DSSF')

    assert np.allclose(decomp_dset, comp_dset)

def test_sample_raw_dataset_compressed():
    lat = [-22.0, -25.0, 20, 17]
    lon = [14.97, 30.0, -15, -5]

    data_dir = 'data/sample_data_files/'

    dssf_fname = data_dir + 'S-LSA_-HDF5_LSASAF_MSG_DSSF_SAfr_201706091400'

    dssf_file = lsasaf.DSSFFile(dssf_fname)
    decomp_dset = dssf_file.sample_raw_dataset('DSSF', lat, lon)

    dssf_fname = data_dir + 'S-LSA_-HDF5_LSASAF_MSG_DSSF_SAfr_201706091400.bz2'

    dssf_file = lsasaf.DSSFFile(dssf_fname)
    comp_dset = dssf_file.sample_raw_dataset('DSSF', lat, lon)

    assert np.allclose(decomp_dset, comp_dset)

if __name__ == '__main__':
    test_read_raw_dataset_compressed()
    test_sample_raw_dataset_compressed()
