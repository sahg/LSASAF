import lsasaf

coff = 1856
loff = 1856

print('\n\n### Geo location supplied\n\n')

lat, lon = (-22.0, 14.97)
row, col = lsasaf.geoloc_to_pixelloc(lat, lon, loff, coff)

# Reverse the transform
lat1, lon1 = lsasaf.pixelloc_to_geoloc(row, col, loff, coff)

print('Requested location (lat, lon) = ({:6f}, {:6f})...'.format(lat, lon))

print('Calculated pixel is at (row, col) = ({:d}, {:d})...'.format(row, col))
print('Calculated location (lat, lon) = ({:6f}, {:6f})...'.format(lat1, lon1))

print('\n\n### Pixel location supplied\n\n')

row, col = (2631, 2356)
lat1, lon1 = lsasaf.pixelloc_to_geoloc(row, col, loff, coff)

# Reverse the transform
row1, col1 = lsasaf.geoloc_to_pixelloc(lat1, lon1, loff, coff)

print('Requested pixel is at (row, col) = ({:d}, {:d})...'.format(row, col))

print('Calculated location (lat, lon) = ({:6f}, {:6f})...'.format(lat1, lon1))
print('Calculated pixel is at (row, col) = ({:d}, {:d})...'.format(row1, col1))
