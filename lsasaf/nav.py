"""A simple module for converting between MSG image pixels and
Geographic co-ordinates. The module uses C code provided by EUMETSAT
(http://www.eumetsat.int/Home/Main/Access_to_Data/User_Support/SP_1117714787347?l=en),
under a freeware license.

This python package is released under the modified BSD license
(details in the included LICENSE.txt).

"""
import numpy as np

from . import _nav

__all__ = ['geoloc_to_pixelloc', 'pixelloc_to_geoloc']

def pixelloc_to_geoloc(row, col, loff, coff):
    """Convert from the pixel row and column to Geographic co-ordinates.

    Convert from the supplied MSG row and column to Geographic
    co-ordinates. The arguments are the MSG image pixel row and col,
    loff and coff are the row and column offsets respectively (usually
    both 1856). A tuple is returned containing the latitude and
    longitude of the pixel centre.

    """
    row = np.atleast_1d(row)
    col = np.atleast_1d(col)

    assert row.shape == col.shape

    lat, lon = _nav.pixelloc_to_geoloc(row, col, loff, coff)

    if lat.size == lon.size == 1:
        # Return scalars
        return lat[0], lon[0]
    else:
        return lat, lon

def geoloc_to_pixelloc(lat, lon, loff, coff):
    """Convert from the Geographic co-ordinates to pixel row and column.

    Convert from the supplied Geographic co-ordinates to MSG row and
    column. The arguments are the Geographic latitude and longitude,
    loff and coff are the row and column offsets respectively (usually
    both 1856). A tuple is returned containing the nearest pixel row
    and column.

    """
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)

    assert lat.shape == lon.shape

    row, col = _nav.geoloc_to_pixelloc(lat, lon, loff, coff)

    if row.size == col.size == 1:
        # Return scalars
        return row[0], col[0]
    else:
        return row, col
