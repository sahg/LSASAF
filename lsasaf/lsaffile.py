"""Utilities for reading LSA-SAF product files.

This module contains some functions for reading the data products produced by
the Land Surface Analysis Satellite Applications Facility
(http://landsaf.meteo.pt). The data products are distributed in HDF5 format and
therefore the module requires the PyTables package. Numpy is also used.

"""
import bz2
import pathlib
import tempfile
from datetime import datetime

import h5py
import numpy as np

from .nav import geoloc_to_pixelloc

__all__ = ['LSAFFile', 'DSSFFile', 'DSLFFile', 'LSTFile']

class LSAFFile:
    """Base class for LSASAF file readers

    """
    def __init__(self, fname):
        self.fname = fname

        self.decomp_path = None
        if pathlib.Path(fname).suffix == '.bz2':
            self.compressed = True
        else:
            self.compressed = False

        self.metadata = {}  # Only evaluate on demand

    def read_metadata(self, dset_name='/'):
        if self.compressed:
            with tempfile.TemporaryDirectory() as tmpdir_name:
                self._decompress_bz2(tmpdir_name)
                self._read_metad(self.decomp_path, dset_name)
        else:
            self._read_metad(self.fname, dset_name)

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
        if self.compressed:
            with tempfile.TemporaryDirectory() as tmpdir_name:
                self._decompress_bz2(tmpdir_name)

                data = self._read_raw_dataset(self.decomp_path, dset_name)
        else:
            data = self._read_raw_dataset(self.fname, dset_name)

        return data

    def sample_raw_dataset(self, dset_name, lat, lon):
        """Sample a raw dataset as it appears on file

        Samples the requested dataset at the given geographic coordinates, as
        it appears in the HDF5 file. No shifting, scaling, masking or datatype
        changes are applied. The data returned by this function must be
        interpreted in conjunction with the dataset metadata and LSASAF
        documentation.

        Parameters
        ----------
        dset_name : string
            The name of the dataset to be read.
        lat : scalar or array
            Latitude/s of the requested sampling location/s.
        lon : scalar or array
            Longitude/s of the requested sampling location/s.

        Returns
        -------
        data : numpy ndarray
            A numpy ndarray object containing the dataset as stored on disk.

        """
        if self.compressed:
            with tempfile.TemporaryDirectory() as tmpdir_name:
                self._decompress_bz2(tmpdir_name)

                data = self._sample_raw_dataset(self.decomp_path,
                                                dset_name, lat, lon)
        else:
            data = self._sample_raw_dataset(self.fname, dset_name, lat, lon)

        return data

    def read_dataset(self, dset_name):
        if self.compressed:
            with tempfile.TemporaryDirectory() as tmpdir_name:
                self._decompress_bz2(tmpdir_name)

                data = self._read_dataset(self.decomp_path, dset_name)
        else:
            data = self._read_dataset(self.fname, dset_name)

        return data

    def sample_dataset(self, dset_name, lat, lon):
        if self.compressed:
            with tempfile.TemporaryDirectory() as tmpdir_name:
                self._decompress_bz2(tmpdir_name)

                data = self._sample_dataset(self.decomp_path,
                                            dset_name, lat, lon)
        else:
            data = self._sample_dataset(self.fname, dset_name, lat, lon)

        return data

    def _read_metad(self, fname, dset_name):
        """Write file metadata into class dict"""
        with h5py.File(fname) as h5file:
            self.metadata[dset_name] = {k: v for k, v in
                                        zip(h5file[dset_name].attrs.keys(),
                                            h5file[dset_name].attrs.values())
                                        }

    def _read_raw_dataset(self, fname, dset_name):
        """Read a raw dataset as it appears on file

        """
        with h5py.File(fname) as h5file:
            data = np.array(h5file[dset_name][...])

            return data

    def _sample_raw_dataset(self, fname, dset_name, lat, lon):
        """Sample a raw dataset as it appears on file

        """
        with h5py.File(fname) as h5file:
            data = h5file[dset_name][...]

            loff = h5file.attrs['LOFF'] - 1
            coff = h5file.attrs['COFF'] - 1

        row, col = geoloc_to_pixelloc(lat, lon, loff, coff)

        return data[row, col]

    def _read_dataset(self, fname, dset_name):
        with h5py.File(fname) as h5file:
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

    def _sample_dataset(self, fname, dset_name, lat, lon):
        with h5py.File(fname) as h5file:
            data = h5file[dset_name][...]

            loff = h5file.attrs['LOFF'] - 1
            coff = h5file.attrs['COFF'] - 1

            offset = h5file[dset_name].attrs['OFFSET']
            scale = h5file[dset_name].attrs['SCALING_FACTOR']
            missing = h5file[dset_name].attrs['MISSING_VALUE']

        row, col = geoloc_to_pixelloc(lat, lon, loff, coff)
        data = data[row, col]

        if (scale is not None) and (offset is not None):
            data = data/scale + offset

        if missing is not None:
            data[data == missing] = np.nan

        return data

    def _decompress_bz2(self, decomp_dir=None):
        comp_path = pathlib.Path(self.fname)

        if decomp_dir:
            # Decompression directory specified
            decomp_path = pathlib.Path(decomp_dir)
            decomp_path = decomp_path.joinpath(comp_path.stem)
        else:
            # Decompress in current working dir
            decomp_path = comp_path.stem

        with bz2.BZ2File(comp_path, 'rb') as bz2f:
            with open(decomp_path, 'wb') as f:
                f.write(bz2f.read())

        self.decomp_path = decomp_path

class DSSFFile(LSAFFile):
    """docstring for DSSFFile."""
    def __init__(self, fname):
        super().__init__(fname)

    def read_dataset(self):  # Override parent class method
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

    def sample_dataset(self, lat, lon):  # Override parent class method
        """Get a masked array containing the DSSF values.

        Sea, space and severly contaminated pixels are masked out. The mask is
        defined according to the bitfield specified in
        SAF_LAND_MF_PUM_DSSF_1.4.pdf. The masked array returned by this
        function contains the DSSF in W/m^2.

        """
        data = super().sample_dataset('/DSSF', lat, lon)
        flags = self.sample_raw_dataset('/DSSF_Q_Flag', lat, lon)

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

class LSTFile(LSAFFile):
    """docstring for LSTFile."""
    def __init__(self, fname):
        super().__init__(fname)

    def read_dataset(self):  # Override parent class version
        """Get a masked array containing the LST values.

        Sea, space and severly contaminated pixels are masked out. The masked
        array returned by this function contains the LST in degrees Centigrade.

        """
        data = super().read_dataset('/LST')
        flags = self.read_raw_dataset('/Q_FLAGS')

        # Mask based on the quality flags [THIS IS STILL IN-PROGRESS]
        data[flags == 0] = np.nan  # sea or space pixel

        data[flags == 0] = np.nan  # sea pixel
        data[flags == 4] = np.nan  # corrupted pixel
        data[flags == 12] = np.nan  # CMa - pixel non processed
        data[flags == 44] = np.nan  # CMa - pixel contaminated by clouds
        data[flags == 60] = np.nan  # CMa - Cloud filled
        data[flags == 76] = np.nan  # CMa - contaminated by snow/ice
        data[flags == 92] = np.nan  # CMa - Undefined
        data[flags == 28] = np.nan  # Emissivity Information Missing
        data[flags == 156] = np.nan  # View Angle Out of Range (EM Poor)
        data[flags == 284] = np.nan  # View Angle Out of Range (EM Nominal)
        data[flags == 412] = np.nan  # View Angle Out of Range (EM Excellent)
        data[flags == 668] = np.nan  # cwv information missing
        data[flags == 796] = np.nan  # cwv information missing
        data[flags == 924] = np.nan  # cwv information missing
        # data[flags == 5790] = np.nan  # Below Nominal (+ EM below nominal)
        # data[flags == 5918] = np.nan  # Below Nominal (+ EM nominal)
        # data[flags == 6046] = np.nan  # Below Nominal (+ EM above nominal)
        # data[flags == 10014] = np.nan  # Nominal (EM nominal)
        # data[flags == 10142] = np.nan  # Nominal (EM above nominal)
        # data[flags == 14238] = np.nan  # Above Nominal (EM above nominal)

        return data
