import matplotlib.pyplot as pl
import numpy as np
import pandas as pd

pl.style.use('defaults.mplstyle')

data = pd.read_excel('Global indicators paper - decadal means.xlsx', sheet_name='Sheet1', header=None)
data.columns = ['year_dec', 'year', 'Decadal mean', 'Annual mean']

# Decadal mean timepoint is one year out. 1850-1859 mean is the mean over the time
# period 1850.01.01 to 1859.12.31, basically the half open interval [1850, 1860).
# Therefore the midpoint is 1855, not 1854 or 1854.5.
data['year_dec'] = data['year_dec'] + 1

print(data['year_dec'])

pl.figure(figsize=(12/2.54, 9/2.54))
pl.plot(data.year, data['Annual mean'], lw=0.75, color='k')
pl.plot(data.year_dec, data['Decadal mean'], lw=2, color='k')
pl.ylabel('Temperature anomaly\nrelative to 1850-1900 (°C)')
pl.xlim(1850, 2025)
pl.ylim(-0.3, 1.35)
pl.tight_layout()
pl.savefig('annual_decadal_mean_gmst.pdf')
pl.savefig('annual_decadal_mean_gmst.png')
pl.show()
