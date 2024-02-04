import akebono
import pytplot

# Load the data
akebono.orb(['1990-02-11', '1990-02-12'])
akebono.vlf_mca(['1990-02-11', '1990-02-12'], datatype='pwr')

# Plot the data
pytplot.tlimit(['1990-02-11 18:00:00', '1990-02-11 18:15:00'])
pytplot.tplot(['akb_mca_Emax_pwr', 'akb_mca_Bmax_pwr'],
              var_label=['akb_orb_alt', 'akb_orb_mlat', 'akb_orb_mlt'])