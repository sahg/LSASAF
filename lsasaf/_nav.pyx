cdef extern int pixcoord2geocoord(int column, int row,
                                  int coff, int loff,
                                  double *latitude, double *longitude)

cdef extern int geocoord2pixcoord(double latitude, double longitude,
                                  int coff, int loff, int *column, int *row)

def pixelloc_to_geoloc(row, col, loff, coff):
    cdef double lat = -999.0
    cdef double lon = -999.0

    pixcoord2geocoord(col, row, coff, loff, &lat, &lon)

    return lat, lon

def geoloc_to_pixelloc(lat, lon, loff, coff):
    cdef int col = -999
    cdef int row = -999

    geocoord2pixcoord(lat, lon, coff, loff, &col, &row)

    return row, col
