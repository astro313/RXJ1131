#!/usr/bin/env python

'''

Organize RXJ1131 IR photometry data using astropy table:
wavelength [um], freq [GHz], flux density [mJy], flux_err [mJy], instrument
Save table as latex file


Last Modified: May 10 2016

History
-------
May 10 2016:
    add decomposed IRAC points, and the missing 3.6 um point
May 1 2016:
    update output error to 4 decimal places, one can truncate the extra sig. fig. in the LaTeX table as needed
Mar 28 2016:
    updated SPITZER/IRAC data to using 5.8" diameter aperture flux
Jan 09 2016:
    combine flux col with flux err col in latex table
    separate out column name with unit in latex table
    added CARAMA, PdBI, VLA C Band continuum data
Jan 08 2016:
    created, added 2MASS, WISE, Spitzer, Herschel, IRAS data

'''

from astropy.table import Table, Column
from astropy.io import ascii
from astropy.constants import c      # m/s
from astropy import units as u
import numpy as np

# (blah * u.Angstrom).to(u.m)
# u.microjansky
# u.microJansky
# u.micrometer = u.micron

# Utilize a row oriented table
#
tableColHead = ('Wavelength [micron]', 'Frequency [GHz]', 'Flux Density [mJy]', 'Flux Err [mJy]', 'Instrument')

data_rows = [(1.25, (1.25*u.micron).to(u.GHz, equivalencies=u.spectral()).value
              , 1.009, 0.090, '2MASS/J-Band'),
             (1.65, (1.65*u.micron).to(u.GHz, equivalencies=u.spectral()).value
              , 1.448, 0.1214, '2MASS/H-Band'),
             (2.17, (2.17*u.micron).to(u.GHz, equivalencies=u.spectral()).value
              , 2.064, 0.1597, '2MASS/Ks-Band'),
             (3.4, (3.4*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              7.027, 0.1424, 'WISE/W1'),
             (4.6, (4.6*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              8.872, 0.1634, 'WISE/W2'),
             (12, (12*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              21.96, 0.4247, 'WISE/W3'),
             (22, (22*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              55.11, 1.878, 'WISE/W4'),
             (3.6, (3.6*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              5.6183, 0.0021, 'Spitzer/IRAC(Extracted)'),
             (4.5, (4.5*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              7.803, 0.002082 , 'Spitzer/IRAC(Archive)'),
             (5.8, (5.8*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              10.72, 0.00512, 'Spitzer/IRAC(Archive)'),
             (8.0, (8.0*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              14.47, 0.004122, 'Spitzer/IRAC(Archive)'),
             (3.6, (3.6*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              5.0336, 0.0021, 'Spitzer/IRAC(Host)'),
             (4.5, (4.5*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              6.0092, 0.0017, 'Spitzer/IRAC(Host)'),
             (5.8, (5.8*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              7.5567, 0.0030, 'Spitzer/IRAC(Host)'),
             (8.0, (8.0*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              9.8811, 0.0039, 'Spitzer/IRAC(Host)'),
             (3.6, (3.6*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              0.5847, 0.00297, 'Spitzer/IRAC(Archive-Host)'),
             (4.5, (4.5*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              1.7938, 0.0027 , 'Spitzer/IRAC(Archive-Host)'),
             (5.8, (5.8*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              3.1633, 0.0059, 'Spitzer/IRAC(Archive-Host)'),
             (8.0, (8.0*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              4.5889, 0.00566, 'Spitzer/IRAC(Archive-Host)'),
             (24, (24*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              47.18, 0.02621, 'Spitzer/MIPS'),
             (250, (250*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              289.4271426, 9.554673, 'Herschel/SPIRE'),
             (350, (350*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              168.229, 8.6147, 'Herschel/SPIRE'),
             (500, (500*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              56.781866, 8.81188, 'Herschel/SPIRE')
             ]

data_IRAS = [(12, (12*u.micron).to(u.GHz, equivalencies=u.spectral()).value              , 0.4e3, None, 'IRAS'),
                 (25, (25*u.micron).to(u.GHz, equivalencies=u.spectral()).value
                  , 0.5e3, None, 'IRAS'),
                 (60, (60*u.micron).to(u.GHz, equivalencies=u.spectral()).value
                  , 0.6e3, None, 'IRAS'),
                 (100, (100*u.micron).to(u.GHz, equivalencies=u.spectral()).value, 1.0e3, None, 'IRAS')]


data_VLA = [((4.8815*u.GHz).to(u.micron, equivalencies=u.spectral()).value,
            4.8815, 1.273277, 4.208577829e-02, 'VLA/Cband-arc'),
            ((4.8815*u.GHz).to(u.micron, equivalencies=u.spectral()).value,
            4.8815,  8.662040e-01, 2.746052e-02, 'VLA/Cband-core')]


data_PdBI = [((139.256*u.GHz).to(u.micron, equivalencies=u.spectral()).value,
              139.256, 1.23, 0.22, 'PdBI-integrated'),
             ((139.256*u.GHz).to(u.micron, equivalencies=u.spectral()).value,
              139.256, 0.799, 1.3e-2, 'PdBI-peak')
            ]

rmsCARMA_cont = 8.30807E-01   # mJy
data_CARMA = [((216*u.GHz).to(u.micron, equivalencies=u.spectral()).value, 216, 3*rmsCARMA_cont, None, 'CARMA')]


# fill table
tbl = Table(rows=data_rows, names=tableColHead,
            dtype=['f4', 'f4', 'f4', 'f4', (str, 30)])
            # dtype=('f8, 'f8', 'f8', 'f8', 'S20')
            # dtype=[(str, 20), (str, 10), int, int, float, float, float, float, float, float, int]
# len(tbl)

def append_table(table, data):
    """
    append row to exiting astropy table

    Parameters
    ----------
    table: astropy.table.Table
        exisiting object

    data: list
        data to append to table

    """

    for i in np.arange(len(data)):
        tbl.add_row(data[i])

append_table(tbl, data_IRAS)
append_table(tbl, data_VLA)
append_table(tbl, data_PdBI)
append_table(tbl, data_CARMA)

formatsTbl = {
              tbl.colnames[0]: lambda x: '{0:.2f}'.format(x),
              tbl.colnames[1]: lambda x: '{0:.1f}'.format(x),
              tbl.colnames[2]: lambda x: '{0:.4f}'.format(x),
              tbl.colnames[3]: lambda x: '{0:.4f}'.format(x),
             }
tbl.write('RXJ1131photometry.dat', format='ascii.commented_header', formats=formatsTbl)
# Equvalently,
# ascii.write(tbl, 'RXJ1131photometry.dat', format='commented_header', formats=formatsTbl)

# ---------------------------------------------------------------------
# LaTexify table
from astropy.io.ascii.latex import latexdicts
from astropy.io.ascii.latex import AASTex

# sort table according to wavelength
import natsort
col_sort = tbl.colnames[0]
tbl = tbl[natsort.index_natsorted(tbl[col_sort])]

### avoid '_' --> subscript,
# same as using kwarg formats={row['Object Name']: lambda x: x.replace("_", "\_")}
# for row in tbl:
#    if "_" in row['Object Name']:
#        row['Object Name'] = row['Object Name'].replace("_", "\_")


### strip off unit in colhead, unit is a separate row in .tex
split1 = [x.find('[') for x in tbl.colnames]
for i, idx in enumerate(split1):
    # if '[' exists in column name
    if idx != -1:
        tbl.rename_column(tbl.colnames[i], tbl.colnames[i][:idx])


def combine_col(table, colname, errcolname, outcolname, unit=None, formatstr="{0:0.3f} $\pm$ {1:s}"):

    """
    combine flux and flux_err column -> with $\pm$

    Parameters
    ----------
    table: astropy.table.Table
        table containing data, and original columns

    colname: str
        of table

    errcolname: str
        of table

    outcolname: str
        column name for the combined column

    unit: astropy unit
        e.g. u.km/u.s

    formatstr: str
        how to format the data in the new column

    Returns
    -------
    astropy.table.Column
        data is str, with # of bits changes depending on data

    """

    # for upper limits, Err column = np.nan
    if np.isnan(table[errcolname]).any():
        tmpCol = table.columns[errcolname]

        # map the formatting function that changes np.nan to \nodata, else 2 decimal place
        # map(f, list)
        tmpCol = map(lambda x: r'\nodata' if np.isnan(x) else '{:.2}'.format(x), tmpCol)

    data = [formatstr.format(aa, err) for aa, err in zip(table[colname], tmpCol)]
    return Column(data=data, name=outcolname, unit=unit)

print tbl.colnames
idx = tbl.colnames.index('Flux Err ')
tableColHeadNew = tbl.colnames
tableColHeadNew.pop(idx)

newTbl = Table([tbl[tbl.colnames[0]],
                tbl[tbl.colnames[1]],
                combine_col(tbl, tbl.colnames[2], tbl.colnames[3], tbl.colnames[2]),
                tbl[tbl.colnames[-1]]])


### format table content
formats = {newTbl.colnames[1]: lambda x: '{0:0.1f}'.format(x)
          }


### setup table Class
demo = False       # testing out different formats

if demo:
    ### simplest, with \begin{table}
    simplest = latexdicts['AA']
    simplest['tablealign'] = 'tbpH'
    simplest['header_start'] = r'\label{tab:photometry}'
    simplest['units'] = {newTbl.colnames[0]: r'\micron', newTbl.colnames[1]: '  GHz', newTbl.colnames[2]: 'mJy', newTbl.colnames[3]: ' '}
    newTbl.write('table_photometry_AA.tex',
                 latexdict=simplest,
                 formats=formats)

    # more flexibility
    flex = latexdicts['template']
    flex['tablealign'] = 'tbpH'
    flex['preamble'] = r'\centering'
    # l['caption'] = ' '   # after preamble; uncomment if want caption above table
    flex['col_align'] = 'lccc \hline'       # default is all c
    flex['header_start'] = ' '
    flex['header_end'] = ' '
    flex['data_start'] = ' '
    flex['data_end'] = ' '
    flex['tablefoot'] = r'\caption{blah blah, \label{tab:blah}}'

    newTbl.colnames
    flex['units'] = {newTbl.colnames[0]: r'\micron', newTbl.colnames[1]: 'GHz',
                     newTbl.colnames[2]: 'mJy', newTbl.colnames[3]: ' '}
    newTbl.write('table_photometry_flex.tex',
                 latexdict=flex)

    # less flexibility, astropy API
    # only have the basics - begin, header, data, end
    newTbl.write('table_photometry_AASsimple.tex', format='ascii.aastex')


# deluxetable with customization
pream = r'\tabletypesize{\scriptsize}' + '\n' + r'\tablecolumns{' + str(len(newTbl.colnames)) + '}\n' + r'\tablecaption{Photometry data}'
cus = {'col_align': 'lccc',
       'tablealign': 'tbpH',
       'preamble': pream,
       'units': {newTbl.colnames[0]: 'micron',
                 newTbl.colnames[1]: 'GHz',
                 newTbl.colnames[2]: 'mJy',
                 newTbl.colnames[3]: ' '},
       'tablefoot': r'\label{tab:BLAH}' + '\n'
                    r'\tablecomments{blah}' + '\n %TablenotegoesBetween \n' +
                    r'\tablerefs{blah}'
       }
# stdout
ascii.write(newTbl, Writer=ascii.AASTex, latexdict=cus)
# save as a tex file
ascii.write(newTbl, 'table_photometry_AAScustomized.tex', Writer=ascii.AASTex, latexdict=cus)


