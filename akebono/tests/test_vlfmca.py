import unittest
from pyspedas.utilities.data_exists import data_exists
import pyspedas
import pytplot

out = pyspedas.akebono.vlf_mca(datatype='dB')
pytplot.tplot(['akb_mca_Emax'])