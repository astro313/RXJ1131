'''

Make CASA region file based on the X,Y positions supplied by user.

Since I want to get the number of beams within which I am calc. the flux from, and I don't know how to do that with go view in GILDAS. I noted the vertices of the polygon I used to extract flux and use CASA viewer to get that area.

Last Modified: 04 May 2016

History:
04 May 16:
    - created

'''

def make_region_file(vertices, outputfile):
    """
    Make a CASA region file for given vertices.
    poly[[x1, y1], [x2, y2], [x3, y3], ...]

    Parameters
    ----------
    vertices : list
        List of direction X and Y vertices in pixel position
    outputfile : str
        Name of output region file
    Returns
    -------
    region_filename : str
        Name of region file
    """

    lines = ['#CRTFv0\n\n']          # to be recognized as a CASA region file
    xylist = []
    X = vertices[0]
    Y = vertices[1]
    for x, y in zip(X, Y):
        xylist.append('[{0}pix, {1}pix]'.format(x, y))   # Units must always be included when defining a region.
    lines.append('poly[{0}]\n'.format(', '.join(xylist)))

    with open(outputfile, 'wb') as f:
        f.writelines(lines)

v = [[133, 118, 122, 153, 150], [150, 138, 104, 114, 141]]
make_region_file(v, 'Flux_CO32_37.2Jy_region.txt')


