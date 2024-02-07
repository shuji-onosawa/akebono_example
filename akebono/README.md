
## Akebono
The routines in this module can be used to load data from the Akebono mission.

### Instruments
- Plasma Waves and Sounder experiment (PWS)
- Radiation Moniter (RDM)
- Orbit data (orb)

### Examples
Get started by importing pyspedas and tplot; these are required to load and plot the data:

```python
import pyspedas
from pytplot import tplot
```

#### Plasma Waves and Sounder experiment (PWS)

```python
pws_vars = pyspedas.akebono.pws(trange=['2012-10-01', '2012-10-02'])

tplot(['akb_pws_RX1', 'akb_pws_RX2'])
```
Remote directory of electron density by 30s data estimated from UHR frequecny is https://www.darts.isas.jaxa.jp/stp/data/exosd/pws/NE/.

#### Radiation Moniter (RDM)

```python
rdm_vars = pyspedas.akebono.rdm(trange=['2012-10-01', '2012-10-02'])

tplot('akb_rdm_FEIO')
```


#### Orbit data (orb)

```python
orb_vars = pyspedas.akebono.orb(trange=['2012-10-01', '2012-10-02'])

tplot(['akb_orb_geo', 'akb_orb_MLT'])
```


### VLF/ Multi-channel Analyzer (MCA)

```python
akebono.vlf_mca(['2012-10-01', '2012-10-02'])
# akebono.vlf_mca(['2012-10-01', '2012-10-02'], datatype='pwr')
tplot(['akb_mca_Emax', 'akb_mca_Bmax'])
# tplot(['akb_mca_Emax_pwr', 'akb_mca_Bmax_pwr'])
```
if you want to change the directory which contains mca data cdf files, please change local_data_dir in config.py.
