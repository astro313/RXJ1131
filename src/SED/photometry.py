'''

Organize RXJ1131 IR photometry data using astropy table
wavelength [µm], freq [GHz], flux density [mJy], flux_err [mJy], instrument


Last Modified: Jan 08 2016


TODO
- add VLA radio data
- add CARMA continuum as upper limit
- add PdBI continuum data
- add code to plot points on SED, put MIR and radio data on plot for background source but don't include in fit, maybe use IRAS 60, 100 µm to constrain SED peak


History
-------
Jan 08 2016:
    created, added 2MASS, WISE, Spitzer, Herschel, IRAS data

'''

from astropy.table import Table, Column
from astropy.constants import c      # m/s
from astropy import units as u
import numpy as np

# (blah * u.Angstrom).to(u.m)
# u.microjansky
# u.microJansky
# u.micrometer = u.micron

# Utilize a row oriented table
#
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
             (4.5, (4.5*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              6.241, 0.00207, 'Spitzer/IRAC'),
             (5.8, (5.8*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              9.354, 0.005694, 'Spitzer/IRAC'),
             (8.9, (8.9*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              13.56, 0.004518, 'Spitzer/IRAC'),
             (24, (24*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              47.18, 0.02621, 'Spitzer/MIPS'),
             (250, (250*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              289.4271426, 9.554673, 'Herschel/SPIRE'),
             (350, (350*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              168.229, 8.6147, 'Herschel/SPIRE'),
             (500, (500*u.micron).to(u.GHz, equivalencies=u.spectral()).value,
              56.781866, 8.81188, 'Herschel/SPIRE')
             ]

data_upperLim = [(12, (12*u.micron).to(u.GHz, equivalencies=u.spectral()).value              , 0.4e3, None, 'IRAS'),
                 (25, (25*u.micron).to(u.GHz, equivalencies=u.spectral()).value
                  , 0.5e3, None, 'IRAS'),
                 (60, (60*u.micron).to(u.GHz, equivalencies=u.spectral()).value
                  , 0.6e3, None, 'IRAS'),
                 (100, (100*u.micron).to(u.GHz, equivalencies=u.spectral()).value, 1.0e3, None, 'IRAS')]


# fill table
tbl = Table(rows=data_rows, names=('Wavelength [micron]', 'Frequency [GHz]', 'Flux density [mJy]', 'Flux_err [mJy]', 'Instrument'), dtype=['f4', 'f4', 'f4', 'f4', (str, 20)])
    #dtype=('f8, 'f8', 'f8', 'f8', 'S20')
    # dtype=[(str, 20), (str, 10), int, int, float, float, float, float, float, float, int]

# tbl.rename_column('a', 'A')
# tbl.colnames
# len(tbl)

# append IRAS data
for i in np.arange(len(data_upperLim)):
    tbl.add_row(data_upperLim[i])

tbl.write('IRphotometry.dat', format='ascii', overwrite=True)


# Latex table format
# perform some operations to conform to latex syntax
# sort according to wavelength?
# change np.nan to \nodata
# combine flux and flux_err column -> with $\pm$
# format to show {:.1f}
# strip off unit in colhead
#for row in tbl:
    # avoid _ --> subscript
#    if "_" in row['Object Name']:
#        row['Object Name'] = row['Object Name'].replace("_", "\_")

from astropy.io.ascii.latex import latexdicts
from astropy.io.ascii.latex import AASTex
from astropy.io import ascii

# setup table Class
# simplest, with \begin{table}
simplest = latexdicts['AA']
simplest['tablealign'] = 'tbpH'
simplest['header_start'] = r'\label{tab:photometry}'
tbl.write('table_photometry_AA.tex',
          latexdict=simplest)

# more flexibility
flex = latexdicts['template']
flex['tablealign'] = 'tbpH'
flex['preamble'] = r'\centering'
# l['caption'] = ' '   # after preamble; uncomment if want caption above table
flex['col_align'] = '|lllcc| \hline'       # default is all c
flex['header_start'] = ' '
flex['header_end'] = ' '
flex['data_start'] = ' '
flex['data_end'] = ' '
flex['tablefoot'] = r'\caption{blah blah, \label{tab:blah}}'

tbl.colnames
flex['units'] = {tbl.colnames[0]: 'micron', tbl.colnames[1]: 'GHz', tbl.colnames[2]: 'mJy', tbl.colnames[3]: 'mJy', tbl.colnames[4]: 'blah'}
tbl.write('table_photometry_flex.tex',
          latexdict=flex)

# less flexibility, astropy API
# only have the basics - begin, header, data, end
tbl.write('table_photometry_AASsimple.tex', format='ascii.aastex')


pream = r'\tabletypesize{\scriptsize}' + '\n' + r'\tablecolumns{' + str(len(tbl.colnames)) + '}\n' + r'\tablecaption{Photometry data}'

# deluxetable with customization
cus = {'col_align': '|rrccc|',
       'tablealign': 'tbpH',
       'preamble': pream,
       'units': {tbl.colnames[0]: 'micron',
                 tbl.colnames[1]: 'GHz',
                 tbl.colnames[2]: 'mJy',
                 tbl.colnames[3]: 'mJy',
                 tbl.colnames[4]: 'blah'},
       'tablefoot': r'\label{tab:BLAH}' + '\n'
                    r'\tablecomments{blah}' + '\n TablenotegoesBetween \n' +
                    r'\tablerefs{blah}'
       }
# stdout
ascii.write(tbl, Writer=ascii.AASTex, latexdict=cus)
# save as a tex file
ascii.write(tbl, 'table_photometry_AAScustomized.tex', Writer=ascii.AASTex, latexdict=cus)


