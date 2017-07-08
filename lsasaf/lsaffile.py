"""Utilities for reading LSA-SAF product files.

This module contains some functions for reading the data products produced by
the Land Surface Analysis Satellite Applications Facility
(http://landsaf.meteo.pt). The data products are distributed in HDF5 format and
therefore the module requires the PyTables package. Numpy is also used.

"""
import os
from bz2 import BZ2File
from datetime import datetime

import h5py
import numpy as np
import tables as h5

from numpy import ma

__all__ = ['LSAFFile', 'DSSFFile']

class LSAFFile:
    """Base class for LSASAF file readers

    """
    def __init__(self, fname):
        self.fname = fname
        self.metadata = {}  # Only evaluate on demand

    def read_metadata(self, dset_name=None):
        if dset_name is None:
            dset_name = '/'
            self._read_metad(dset_name)
        else:
            self._read_metad(dset_name)

        return self.metadata[dset_name]

    def slot_time(self):
        """Parse an LSA-SAF file name to get the slot time.

        A datetime object containing the slot time (in UTC) is returned.

        """
        # HDF5_LSASAF_MSG_DSLF_SAfr_200702211000.h5 or
        # S-LSA_-HDF5_LSASAF_MSG_DSLF_SAfr_200707260830 etc.
        indx = self.fname.rfind('_') + 1
        year = int(self.fname[indx:indx+4])
        month = int(self.fname[indx+4:indx+6])
        day = int(self.fname[indx+6:indx+8])
        hour = int(self.fname[indx+8:indx+10])
        minute = int(self.fname[indx+10:indx+12])

        return datetime(year, month, day, hour, minute)

    def read_raw_dataset(self, dset_name):
        """Read a raw dataset as it appears on file

        Reads the requested dataset as it appears in the HDF5 file. No
        shifting, scaling, masking or datatype changes are applied. The data
        returned by this function must be interpreted in conjunction with the
        dataset metadata and LSASAF documentation.

        Parameters
        ----------
        dset_name : string
            The name of the dataset to be read.

        Returns
        -------
        data : numpy ndarray

            A numpy ndarray object containing the dataset as stored on disk.
            The datatype is as stored and no shifting, scaling or masking is
            applied.

        """
        with h5py.File(self.fname) as h5file:
            data = np.array(h5file[dset_name][...])

        return data

    def read_dataset(self, dset_name):
        with h5py.File(self.fname) as h5file:
            ds = h5file[dset_name]

            data = ds[...]

            offset = ds.attrs.get("OFFSET")
            scale = ds.attrs.get("SCALING_FACTOR")
            missing = ds.attrs.get("MISSING_VALUE")

            if (scale is not None) and (offset is not None):
                data = data/scale + offset

            if missing is not None:
                data[data == missing] = np.nan

        return data

    def _read_metad(self, dset_name):
        """Write file metadata into class dict"""
        with h5py.File(self.fname) as h5file:
            self.metadata[dset_name] = {k: v for k, v in
                                        zip(h5file[dset_name].attrs.keys(),
                                            h5file[dset_name].attrs.values())
                                        }

class DSSFFile(LSAFFile):
    """docstring for DSSFFile."""
    def __init__(self, fname):
        super().__init__(fname)

    def read_dataset(self):  # Override parent class version
        """Get a masked array containing the DSSF values.

        Sea, space and severly contaminated pixels are masked out. The mask is
        defined according to the bitfield specified in
        SAF_LAND_MF_PUM_DSSF_1.4.pdf. The masked array returned by this
        function contains the DSSF in W/m^2.

        """
        data = super().read_dataset('/DSSF')
        flags = self.read_raw_dataset('/DSSF_Q_Flag')

        # Mask based on the quality flags [THIS IS STILL IN-PROGRESS]
        data[flags == 0] = np.nan  # ocean pixel
        data[flags == 2] = np.nan  # space pixel

        return data

class DSLFFile(LSAFFile):
    """docstring for DSLFFile."""
    def __init__(self, fname):
        super().__init__(fname)

    def read_dataset(self):  # Override parent class version
        """Get a masked array containing the DSLF values.

        Sea, space and severly contaminated pixels are masked out. The masked
        array returned by this function contains the DSLF in W/m^2.

        """
        data = super().read_dataset('/DSLF')
        flags = self.read_raw_dataset('/Q_FLAGS')

        # Mask based on the quality flags [THIS IS STILL IN-PROGRESS]
        data[flags == 0] = np.nan  # sea or space pixel
        data[flags == 4] = np.nan  # T2m missing
        data[flags == 12] = np.nan  # CMa - pixel non processed
        data[flags == 92] = np.nan  # CMa - Undefined
        data[flags == 156] = np.nan  # TPW information missing
        data[flags == 44] = np.nan
        # ^^ CTTH_EFFECTIVE missing (CMa - pixel contaminated by clouds)
        data[flags == 60] = np.nan
        # ^^ CTTH_EFFECTIVE missing (CMa - Cloud filled)
        data[flags == 76] = np.nan
        # ^^ CTTH_EFFECTIVE missing (CMa - contaminated by snow/ice)
        data[flags == 812] = np.nan
        # ^^ Td2m missing (CMa - pixel contaminated by clouds)
        data[flags == 828] = np.nan  # Td2m missing (CMa - Cloud filled)
        data[flags == 844] = np.nan
        # ^^ Td2m missing (CMa - contaminated by snow/ice)
    #    data[flags == 11422] = np.nan  # Below Nominal (CMa - Cloud-free)
    #    data[flags == 19614] = np.nan  # Nominal (CMa - Cloud-free)
    #    data[flags == 27806] = np.nan  # Above Nominal (CMa - Cloud-free)
    #    data[flags == 13102] = np.nan
    #    # ^^ Below Nominal (CMa - pixel contaminated by clouds)
    #    data[flags == 21294] = np.nan
    #    # ^^ Nominal (CMa - pixel contaminated by clouds)
    #    data[flags == 29486] = np.nan
    #    # ^^ Above Nominal (CMa - pixel contaminated by clouds)
    #    data[flags == 13118] = np.nan  # Below Nominal (CMa - Cloud filled)
    #    data[flags == 21310] = np.nan  # Nominal (CMa - Cloud filled)
    #    data[flags == 29502] = np.nan  # Above Nominal (CMa - Cloud filled)
    #    data[flags == 13134] = np.nan
    #    # ^^ Below Nominal (CMa - contaminated by snow/ice)
    #    data[flags == 21326] = np.nan
    #    # ^^ Nominal(CMa - contaminated by snow/ice)
    #    data[flags == 29518] = np.nan
    #    # ^^ Above Nominal (CMa - contaminated by snow/ice)

        return data

def _read_raw(file_name, data_node_name, quality_node_name):
    """Return the raw data and quality control flags.

    This function returns the data as stored in the HDF5 data and q_flag
    arrays. The scaling factors are applied. Use this function if you need to
    do your own (non-standard) masking of the LSA-SAF data. Numpy arrays are
    returned with the same shape as the HDF5 data arrays. The returned data
    array has type float32 and the flags array has the same type as the data in
    the HDF5 file.

    """
    h5file = h5.openFile(file_name)

    node = h5file.getNode(data_node_name)
    data = node.read()
    data = np.asarray(data, np.float32)
    if (node._v_attrs.SCALING_FACTOR != 1):
        data /= node._v_attrs.SCALING_FACTOR

    node = h5file.getNode(quality_node_name)
    flags = node.read()

    h5file.close()

    return data, flags

def read_lst(file_name):
    """Get a masked array containing the LST values.

    Sea, space and severly contaminated pixels are masked out. The masked array
    returned by this function contains the LST in degrees Centigrade.

    """
    # _read_raw() requires an uncompressed HDF5 file
    if file_name[-3:] == 'bz2':
        # create a temp file
        temp_fname = 'temp.h5'

        bz2_file = BZ2File(file_name)
        fp = open(temp_fname, 'wb')
        fp.write(bz2_file.read())
        fp.close()
        bz2_file.close()

        data, flags = _read_raw(temp_fname, '/LST', '/Q_FLAGS')

        os.remove(temp_fname)
    else:
        data, flags = _read_raw(file_name, '/LST', '/Q_FLAGS')

    # mask based on the quality flags
    data = ma.masked_where(flags == 0, data)# sea pixel
    data = ma.masked_where(flags == 4, data)# corrupted pixel
    data = ma.masked_where(flags == 12, data)# CMa - pixel non processed
    data = ma.masked_where(flags == 44, data)# CMa - pixel contaminated by clouds
    data = ma.masked_where(flags == 60, data)# CMa - Cloud filled
    data = ma.masked_where(flags == 76, data)# CMa - contaminated by snow/ice
    data = ma.masked_where(flags == 92, data)# CMa - Undefined
    data = ma.masked_where(flags == 28, data)# Emissivity Information Missing
    data = ma.masked_where(flags == 156, data)# Viewing Angle Out of Range (EM Poor)
    data = ma.masked_where(flags == 284, data)# Viewing Angle Out of Range (EM Nominal)
    data = ma.masked_where(flags == 412, data)# Viewing Angle Out of Range (EM Excellent)
    data = ma.masked_where(flags == 668, data)# cwv information missing
    data = ma.masked_where(flags == 796, data)# cwv information missing
    data = ma.masked_where(flags == 924, data)# cwv information missing
##    data = ma.masked_where(flags == 5790, data)# Below Nominal (+ EM below nominal)
##    data = ma.masked_where(flags == 5918, data)# Below Nominal (+ EM nominal)
##    data = ma.masked_where(flags == 6046, data)# Below Nominal (+ EM above nominal)
##    data = ma.masked_where(flags == 10014, data)# Nominal (EM nominal)
##    data = ma.masked_where(flags == 10142, data)# Nominal (EM above nominal)
##    data = ma.masked_where(flags == 14238, data)# Above Nominal (EM above nominal)

    return data
