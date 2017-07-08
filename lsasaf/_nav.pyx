import numpy as np

cdef extern int pixcoord2geocoord(int column, int row,
                                  int coff, int loff,
                                  double *latitude, double *longitude)

cdef extern int geocoord2pixcoord(double latitude, double longitude,
                                  int coff, int loff, int *column, int *row)

def pixelloc_to_geoloc(row, col, int loff, int coff):
    """Convert MSG grid locations to Geographic locations

    """
    cdef int k = 0
    cdef double lat = -999.0
    cdef double lon = -999.0

    cdef int nlocs = row.size

    lat_arr = np.zeros(nlocs, dtype=np.double)
    lon_arr = np.zeros(nlocs, dtype=np.double)

    for k in range(nlocs):
        pixcoord2geocoord(col[k], row[k], coff, loff, &lat, &lon)

        lat_arr[k] = lat
        lon_arr[k] = lon

    return lat_arr, lon_arr

def geoloc_to_pixelloc(lat, lon, int loff, int coff):
    """Convert Geographic locations to MSG grid locations

    """
    cdef int k = 0
    cdef int c = -999
    cdef int r = -999

    cdef int nlocs = lat.size

    row = np.zeros(nlocs, dtype=np.int)
    col = np.zeros(nlocs, dtype=np.int)

    for k in range(nlocs):
        geocoord2pixcoord(lat[k], lon[k], coff, loff, &c, &r)

        row[k] = r
        col[k] = c

    return row, col
