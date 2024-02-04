from .load import load
import numpy as np
import pandas as pd
from pytplot import store_data, options, get_data, tplot_names, del_data
from pyspedas import time_double
from pyspedas.cotrans.xyz_to_polar import xyz_to_polar
from fnmatch import filter


def pws(trange=['2012-10-01', '2012-10-02'],
        datatype='ne',
        level='h1',
        suffix='',
        get_support_data=False,
        varformat=None,
        varnames=[],
        downloadonly=False,
        notplot=False,
        no_update=False,
        time_clip=False):
    """
    This function loads data from the Plasma Waves and Sounder experiment (PWS)

    Parameters
    ----------
        trange : list of str
            time range of interest [starttime, endtime] with the format
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        datatype: str
            Data type; Valid options:
                'ne', 'npw-ds', 'npw-py', 'spw'

        level: str
            Data level; options: 'h1' (default: h1)

        suffix: str
            The tplot variable names will be given this suffix.  By default,
            no suffix is added.

        get_support_data: bool
            Data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into tplot.  By default, only loads in data with a
            "VAR_TYPE" attribute of "data".

        varformat: str
            The file variable formats to load into tplot.  Wildcard character
            "*" is accepted.  By default, all variables are loaded in.

        varnames: list of str
            List of variable names to load (if not specified,
            all data variables are loaded)

        downloadonly: bool
            Set this flag to download the CDF files, but not load them into
            tplot variables

        notplot: bool
            Return the data in hash tables instead of creating tplot variables

        no_update: bool
            If set, only load data from your local cache

        time_clip: bool
            Time clip the variables to exactly the range specified in the trange keyword

    Returns
    ----------
        List of tplot variables created.

    """
    tvars = load(instrument='pws', trange=trange, level=level, datatype=datatype, suffix=suffix, get_support_data=get_support_data, varformat=varformat, varnames=varnames, downloadonly=downloadonly, notplot=notplot, time_clip=time_clip, no_update=no_update)

    if tvars is None or notplot or downloadonly:
        return tvars

    return pws_postprocessing(tvars)


def pws_postprocessing(variables):
    """
    Placeholder for PWS post-processing
    """
    return variables


def rdm(trange=['2012-10-01', '2012-10-02'],
        suffix='',
        get_support_data=False,
        varformat=None,
        varnames=[],
        downloadonly=False,
        notplot=False,
        no_update=False,
        time_clip=False):
    """
    This function loads data from the Radiation Moniter (RDM)

    Parameters
    ----------
        trange : list of str
            time range of interest [starttime, endtime] with the format
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        suffix: str
            The tplot variable names will be given this suffix.  By default,
            no suffix is added.

        get_support_data: bool
            Data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into tplot.  By default, only loads in data with a
            "VAR_TYPE" attribute of "data".

        varformat: str
            The file variable formats to load into tplot.  Wildcard character
            "*" is accepted.  By default, all variables are loaded in.

        varnames: list of str
            List of variable names to load (if not specified,
            all data variables are loaded)

        downloadonly: bool
            Set this flag to download the CDF files, but not load them into
            tplot variables

        notplot: bool
            Return the data in hash tables instead of creating tplot variables

        no_update: bool
            If set, only load data from your local cache

        time_clip: bool
            Time clip the variables to exactly the range specified in the trange keyword

    Returns
    ----------
        List of tplot variables created.

    """
    files = load(instrument='rdm', trange=trange, suffix=suffix, get_support_data=get_support_data, varformat=varformat, varnames=varnames, downloadonly=downloadonly, notplot=notplot, time_clip=time_clip, no_update=no_update)

    if files is None or notplot or downloadonly:
        return files

    return rdm_postprocessing(files)


def rdm_postprocessing(files):
    """
    Load the RDM ASCII files into tplot variables
    """
    data = load_csv_file(files)
    values = data.to_numpy()
    unix_times = time_double([ymd + '/' + hms for ymd, hms in zip(values[:, 0], values[:, 1])])

    L = np.float64(values[:, 2])
    INV = np.float64(values[:, 3])
    FMLAT = np.float64(values[:, 4])
    MLAT = np.float64(values[:, 5])
    MLT = np.float64(values[:, 6])
    ALT = np.float64(values[:, 7])
    GLAT = np.float64(values[:, 8])
    GLON = np.float64(values[:, 9])
    RDM_E3 = np.float64(values[:, 10])
    Energy = np.zeros(len(RDM_E3))
    Energy[:] = 2.5

    prefix_project = 'akb_'
    prefix_descriptor = 'rdm_'
    prefix = prefix_project + prefix_descriptor

    store_data(prefix_project+'L', data={'x': unix_times, 'y': L})
    store_data(prefix_project+'INV', data={'x': unix_times, 'y': INV})
    store_data(prefix_project+'FMLAT', data={'x': unix_times, 'y': FMLAT})
    store_data(prefix_project+'MLAT', data={'x': unix_times, 'y': MLAT})
    store_data(prefix_project+'MLT', data={'x': unix_times, 'y': MLT})
    store_data(prefix_project+'ALT', data={'x': unix_times, 'y': ALT})
    store_data(prefix_project+'GLAT', data={'x': unix_times, 'y': GLAT})
    store_data(prefix_project+'GLON', data={'x': unix_times, 'y': GLON})
    store_data(prefix+'FEIO', data={'x': unix_times, 'y': RDM_E3})
    store_data(prefix+'FEIO_Energy', data={'x': unix_times, 'y': Energy})

    options(prefix+'FEIO', 'spec', True)

    options(prefix_project+'L', 'ytitle', 'L-value')
    options(prefix_project+'INV', 'ytitle', 'Invariant Latitude [deg]')
    options(prefix_project+'FMLAT', 'ytitle', 'Footprint Latitude [deg]')
    options(prefix_project+'MLAT', 'ytitle', 'Magnetic Latitude [deg]')
    options(prefix_project+'MLT', 'ytitle', 'Magnetic Local Time [hour]')
    options(prefix_project+'ALT', 'ytitle', 'Altitude [km]')
    options(prefix_project+'GLAT', 'ytitle', 'Geographic Latitude [deg]')
    options(prefix_project+'GLON', 'ytitle', 'Geographic Longitude [deg]')
    options(prefix+'FEIO', 'ytitle', 'Omni-directional Integral Electron Flux')
    options(prefix+'FEIO', 'ysubtitle', '[/cm^22 sec str]')
    options(prefix+'FEIO_Energy', 'ytitle', 'Elctron energy [MeV]')

    return [prefix_project+'L',
            prefix_project+'INV',
            prefix_project+'FMLAT',
            prefix_project+'MLAT',
            prefix_project+'MLT',
            prefix_project+'ALT',
            prefix_project+'GLAT',
            prefix_project+'GLON',
            prefix+'FEIO',
            prefix+'FEIO_Energy']


def orb(trange=['2012-10-01', '2012-10-02'],
        suffix='',
        get_support_data=False,
        varformat=None,
        varnames=[],
        downloadonly=False,
        notplot=False,
        no_update=False,
        time_clip=False):
    """
    This function loads data from the Akebono orbit data (orb)

    Parameters
    ----------
        trange : list of str
            time range of interest [starttime, endtime] with the format
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        suffix: str
            The tplot variable names will be given this suffix.  By default,
            no suffix is added.

        get_support_data: bool
            Data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into tplot.  By default, only loads in data with a
            "VAR_TYPE" attribute of "data".

        varformat: str
            The file variable formats to load into tplot.  Wildcard character
            "*" is accepted.  By default, all variables are loaded in.

        varnames: list of str
            List of variable names to load (if not specified,
            all data variables are loaded)

        downloadonly: bool
            Set this flag to download the CDF files, but not load them into
            tplot variables

        notplot: bool
            Return the data in hash tables instead of creating tplot variables

        no_update: bool
            If set, only load data from your local cache

        time_clip: bool
            Time clip the variables to exactly the range specified in the trange keyword

    Returns
    ----------
        List of tplot variables created.
        ['pass','ut', 'ksc_azm', 'ksc_elv', 'ksc_dis', 'ksc_ang', 'syo_azm', 'syo_elv', 'syo_dis', 'syo_ang',
            'pra_azm', 'pra_elv', 'pra_dis', 'pra_ang', 'esr_azm', 'esr_elv', 'esr_dis', 'esr_ang', 'gclat','gclon',
            'inv', 'fmlat', 'mlat', 'mlt', 'bmdl_x', 'bmdl_y', 'bmdl_z', 'xxlon_sc', 'xxlat_sc', 'aheight','lsun',
            's_direc_x', 's_direc_y', 's_direc_z', 'sc_pos_x', 'sc_pos_y', 'sc_pos_z', 'sc_vel_x', 'sc_vel_y', 'sc_vel_z']
    """
    files = load(instrument='orb', trange=trange, suffix=suffix, get_support_data=get_support_data, varformat=varformat, varnames=varnames, downloadonly=downloadonly, notplot=notplot, time_clip=time_clip, no_update=no_update)

    if files is None or notplot or downloadonly:
        return files

    return orb_postprocessing(files)


def orb_postprocessing(files):
    """
    Load the orbit CSV files and create the tplot variables
    """
    prefix_project = 'akb_'
    prefix_descriptor = 'orb_'
    prefix = prefix_project + prefix_descriptor

    cols = ['pass','ut', 'ksc_azm', 'ksc_elv', 'ksc_dis', 'ksc_ang', 'syo_azm', 'syo_elv', 'syo_dis', 'syo_ang',
            'pra_azm', 'pra_elv', 'pra_dis', 'pra_ang', 'esr_azm', 'esr_elv', 'esr_dis', 'esr_ang', 'gclat','gclon',
            'inv', 'fmlat', 'mlat', 'mlt', 'bmdl_x', 'bmdl_y', 'bmdl_z', 'xxlon_sc', 'xxlat_sc', 'aheight','lsun',
            's_direc_x', 's_direc_y', 's_direc_z', 'sc_pos_x', 'sc_pos_y', 'sc_pos_z', 'sc_vel_x', 'sc_vel_y', 'sc_vel_z']

    data = load_csv_file(files, cols=cols)
    original_time = data['ut'].to_numpy()
    converted_time = convert_orb_time(original_time)
    unix_times = time_double(converted_time)

    km_in_re = 6374.4

    xyz = np.array([[data['sc_pos_x']], [data['sc_pos_y']], [data['sc_pos_z']]]).transpose([2, 0, 1]).squeeze()
    xyz = np.float64(xyz)
    xyz_re = xyz/km_in_re
    r_theta_phi = xyz_to_polar(xyz)
    rr = r_theta_phi[:, 0]
    th = r_theta_phi[:, 1]
    ph = r_theta_phi[:, 2]
    bmdl = np.array([[data['bmdl_x']], [data['bmdl_y']], [data['bmdl_z']]]).transpose([2, 0, 1]).squeeze()
    bmdl = np.float64(bmdl)
    bmdl_scaler = np.sqrt(np.sum(bmdl**2, axis=1))
    store_data(prefix + 'geo', data={'x': unix_times, 'y': xyz_re})
    store_data(prefix + 'gdlat', data={'x': unix_times, 'y': np.float64(data['gclat'])})
    store_data(prefix + 'gdlon', data={'x': unix_times, 'y': np.float64(data['gclon'])})
    store_data(prefix + 'inv', data={'x': unix_times, 'y': np.float64(data['inv'])})
    store_data(prefix + 'fmlat', data={'x': unix_times, 'y': np.float64(data['fmlat'])})
    store_data(prefix + 'mlat', data={'x': unix_times, 'y': np.float64(data['mlat'])})
    store_data(prefix + 'mlt', data={'x': unix_times, 'y': np.float64(data['mlt'])})
    store_data(prefix + 'alt', data={'x': unix_times, 'y': np.float64(data['aheight'])})
    store_data(prefix + 'gcalt', data={'x': unix_times, 'y': rr / km_in_re})
    store_data(prefix + 'gclat', data={'x': unix_times, 'y': th})
    store_data(prefix + 'gclon', data={'x': unix_times, 'y': ph})
    store_data(prefix + 'bmdl_scaler', data={'x': unix_times, 'y': bmdl_scaler})
    options(prefix + 'geo', 'ytitle', 'GEO')
    options(prefix + 'geo', 'ysubtitle', '[Re]')
    options(prefix + 'gdlat', 'ytitle', 'Geodetic latitude of the magnetic footprint')
    options(prefix + 'gdlat', 'ysubtitle', '(120km altitude) [deg]')
    options(prefix + 'gdlon', 'ytitle', 'Geodetic longitude of the magnetic footprint')
    options(prefix + 'gdlon', 'ysubtitle', '(120km altitude) [deg]')
    options(prefix + 'inv', 'ytitle', 'Invariant Latitude of the magnetic footprint')
    options(prefix + 'inv', 'ysubtitle', '(120km altitude) [deg]')
    options(prefix + 'fmlat', 'ytitle', 'Geomagnetic Latitude of the magnetic footprint')
    options(prefix + 'fmlat', 'ysubtitle', '(120km altitude) [deg]')
    options(prefix + 'mlat', 'ytitle', 'Magnetic Latitude of the spacecraft position')
    options(prefix + 'mlat', 'ysubtitle', '[deg]')
    options(prefix + 'mlt', 'ytitle', 'Magnetic Local Time')
    options(prefix + 'mlt', 'ysubtitle', '[hours]')
    options(prefix + 'alt', 'ytitle', 'Altitude')
    options(prefix + 'alt', 'ysubtitle', '[km]')
    options(prefix + 'gcalt', 'ytitle', 'Geocentric Altitude')
    options(prefix + 'gcalt', 'ysubtitle', '[Re]')
    options(prefix + 'gclat', 'ytitle', 'Geocentric Latitude')
    options(prefix + 'gclat', 'ysubtitle', '[deg]')
    options(prefix + 'gclon', 'ytitle', 'Geocentric Longitude')
    options(prefix + 'gclon', 'ysubtitle', '[deg]')

    stored_names = filter(tplot_names(quiet=True),
                          prefix_project + prefix_descriptor + '*')

    return stored_names


def convert_orb_time(original_time):
    """
    Convert the time format of the orb data
    original_time: str or numpy.ndarray of str
        The time in the original format. The format is 'yyMMddHHmmss'

    Returns
    ----------
        The time in the format 'yyyy-mm-dd/hh:mm:ss'
    """
    if isinstance(original_time, str):
        if 89 <= int(original_time[0:2]) <= 99:
            return '19' + original_time[0:2] + '-' + original_time[2:4] + '-' + original_time[4:6] + '/' + original_time[6:8] + ':' + original_time[8:10] + ':' + original_time[10:12]
        else:
            return '20' + original_time[0:2] + '-' + original_time[2:4] + '-' + original_time[4:6] + '/' + original_time[6:8] + ':' + original_time[8:10] + ':' + original_time[10:12]
    elif isinstance(original_time, np.ndarray):
        if 89 <= int(original_time[0][0:2]) <= 99:
            return np.array(['19' + time[0:2] + '-' + time[2:4] + '-' + time[4:6] + '/' + time[6:8] + ':' + time[8:10] + ':' + time[10:12] for time in original_time])
        else:
            return np.array(['20' + time[0:2] + '-' + time[2:4] + '-' + time[4:6] + '/' + time[6:8] + ':' + time[8:10] + ':' + time[10:12] for time in original_time])


def load_csv_file(filenames, cols=None):
    """
    Loads a list of CSV/txt files into pandas data frames
    """
    if not isinstance(filenames, list):
        filenames = [filenames]
    df = pd.concat((pd.read_csv(f, header=0, delim_whitespace=True, dtype=str, names=cols) for f in filenames), ignore_index=True)
    return df


def vlf_mca(trange=['2012-10-01', '2012-10-02'],
            datatype='dB',
            del_invalid_data=None,
            suffix='',
            get_support_data=False,
            varformat=None,
            varnames=[],
            downloadonly=False,
            notplot=False,
            no_update=False,
            time_clip=False):
    """
    This function loads data from the Plasma Waves and Sounder experiment (PWS)

    Parameters
    ----------
        trange : list of str
            time range of interest [starttime, endtime] with the format
            'YYYY-MM-DD','YYYY-MM-DD'] or to specify more or less than a day
            ['YYYY-MM-DD/hh:mm:ss','YYYY-MM-DD/hh:mm:ss']

        datatype: str
            Data type; Valid options:
                'dB': decibel, 0[dB]=10^-6[mV/m] for E-field and 0[dB]=10^-6[pT] for B-field,
                'amp': mV/m/Hz^1/2 or nT/Hz^1/2
                'pwr': (mV/m)^2/Hz or nT^2/Hz

        del_invalid_data :  list of string.
            mca cdf contain data from which the interference by BDR or SMS is *not* yet removed.
            You can remove data contaminated by interference by passing a list containing the following words.
            'off': mca is off
            'noisy': data is noisy
            'sms': SMS on
            'bdr': BDR on
            'bit rate m': Bit rate is medium. When the bit rate is medium, the data is not reliable.
            'pws': PWS sounder on
        suffix: str
            The tplot variable names will be given this suffix.  By default,
            no suffix is added.

        get_support_data: bool
            Data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into tplot.  By default, only loads in data with a
            "VAR_TYPE" attribute of "data".

        varformat: str
            The file variable formats to load into tplot.  Wildcard character
            "*" is accepted.  By default, all variables are loaded in.

        varnames: list of str
            List of variable names to load (if not specified,
            all data variables are loaded)

        downloadonly: bool
            Set this flag to download the CDF files, but not load them into
            tplot variables

        notplot: bool
            Return the data in hash tables instead of creating tplot variables

        no_update: bool
            If set, only load data from your local cache

        time_clip: bool
            Time clip the variables to exactly the range specified in the trange keyword

    Returns
    ----------
        List of tplot variables created.
    """

    tvars = load(instrument='mca', trange=trange, datatype=datatype, suffix=suffix, get_support_data=get_support_data, varformat=varformat, varnames=varnames, downloadonly=downloadonly, notplot=notplot, time_clip=time_clip, no_update=no_update)

    if tvars is None or notplot or downloadonly:
        return tvars

    return mca_postprocessing(datatype, del_invalid_data)


def mca_postprocessing(datatype, del_invalid_data):
    prefix = 'akb_mca_'
    if not del_invalid_data:
        pass
    else:
        Emax, Bmax, Eave, Bave = get_data(prefix+'Emax'), get_data(prefix+'Bmax'), get_data(prefix+'Eave'), get_data(prefix+'Bave')
        Emax_array, Bmax_array, Eave_array, Bave_array = Emax.y.astype(float), Bmax.y.astype(float), Eave.y.astype(float), Bave.y.astype(float)
        postgap = get_data(prefix+'PostGap')
        postgap_arr = postgap.y.T
        postgap_bin_arr = np.array([list(np.binary_repr(x, width=8)) for x in postgap_arr], dtype=int)
        for inst_name in del_invalid_data:
            inst_name = inst_name.lower()
            if inst_name in ['off', 'noisy', 'bdr', 'sms', 'bit rate m', 'pws']:
                pass
            else:
                raise Exception('del_invalid_data list must consist of either off, noisy, bdr, sms, bit rate m or pws')

            if inst_name == 'off':
                nodata_index_arr = np.nonzero(postgap_bin_arr[:, 7])[0]
                Emax_array[nodata_index_arr] = np.nan
                Bmax_array[nodata_index_arr] = np.nan
                Eave_array[nodata_index_arr] = np.nan
                Bave_array[nodata_index_arr] = np.nan
            if inst_name == 'noisy':
                noisy_index_arr = np.nonzero(postgap_bin_arr[:, 6])[0]
                Emax_array[noisy_index_arr] = np.nan
                Bmax_array[noisy_index_arr] = np.nan
                Eave_array[noisy_index_arr] = np.nan
                Bave_array[noisy_index_arr] = np.nan
            if inst_name == 'bdr':
                bdr_index_arr = np.nonzero(postgap_bin_arr[:, 3])[0]
                Emax_array[bdr_index_arr] = np.nan
                Bmax_array[bdr_index_arr] = np.nan
                Eave_array[bdr_index_arr] = np.nan
                Bave_array[bdr_index_arr] = np.nan
            if inst_name == 'sms':
                sms_index_arr = np.nonzero(postgap_bin_arr[:, 2])[0]
                Emax_array[sms_index_arr] = np.nan
                Bmax_array[sms_index_arr] = np.nan
                Eave_array[sms_index_arr] = np.nan
                Bave_array[sms_index_arr] = np.nan
            if inst_name == 'bit rate m':
                bitrate_index_arr = np.nonzero(postgap_bin_arr[:, 1])[0]
                Emax_array[bitrate_index_arr] = np.nan
                Bmax_array[bitrate_index_arr] = np.nan
                Eave_array[bitrate_index_arr] = np.nan
                Bave_array[bitrate_index_arr] = np.nan
            if inst_name == 'pws':
                pws_index_arr = np.nonzero(postgap_bin_arr[:, 0])[0]
                Emax_array[pws_index_arr] = np.nan
                Bmax_array[pws_index_arr] = np.nan
                Eave_array[pws_index_arr] = np.nan
                Bave_array[pws_index_arr] = np.nan

        store_data(prefix+'Emax',
                   data={'x': Emax.times, 'y': Emax_array, 'v': Emax.v})
        store_data(prefix+'Bmax',
                   data={'x': Bmax.times, 'y': Bmax_array, 'v': Bmax.v})
        store_data(prefix+'Eave',
                   data={'x': Eave.times, 'y': Eave_array, 'v': Eave.v})
        store_data(prefix+'Bave',
                   data={'x': Bave.times, 'y': Bave_array, 'v': Bave.v})

    if datatype == 'dB':
        options(prefix+'Emax',
                opt_dict={'spec': 1, 'ylog': 1, 'zlog': 1,
                          'yrange': [1, 2e4], 'ytitle': 'Emax', 'ysubtitle': 'Frecuency [Hz]',
                          'ztitle': 'dB'})
        options(prefix+'Bmax',
                opt_dict={'spec': 1, 'ylog': 1, 'zlog': 1,
                          'yrange': [1, 2e4], 'ytitle': 'Bmax', 'ysubtitle': 'Frecuency [Hz]',
                          'ztitle': 'dB'})
        options(prefix+'Eave',
                opt_dict={'spec': 1, 'ylog': 1, 'zlog': 1,
                          'yrange': [1, 2e4], 'ytitle': 'Eave', 'ysubtitle': 'Frecuency [Hz]',
                          'ztitle': 'dB'})
        options(prefix+'Bave',
                opt_dict={'spec': 1, 'ylog': 1, 'zlog': 1,
                          'yrange': [1, 2e4], 'ytitle': 'Bave', 'ysubtitle': 'Frecuency [Hz]',
                          'ztitle': 'dB'})
        return [prefix+'Emax',
                prefix+'Bmax',
                prefix+'Eave',
                prefix+'Bave']
    elif datatype == 'pwr':
        mca_h1cdf_dB_to_absolute('pwr')
        del_data(prefix+'Emax')
        del_data(prefix+'Bmax')
        del_data(prefix+'Eave')
        del_data(prefix+'Bave')
        return [prefix+'Emax_pwr',
                prefix+'Bmax_pwr',
                prefix+'Eave_pwr',
                prefix+'Bave_pwr']
    elif datatype == 'amp':
        mca_h1cdf_dB_to_absolute('amp')
        del_data(prefix+'Emax')
        del_data(prefix+'Bmax')
        del_data(prefix+'Eave')
        del_data(prefix+'Bave')
        return [prefix+'Emax_amp',
                prefix+'Bmax_amp',
                prefix+'Eave_amp',
                prefix+'Bave_amp']


def dB_to_absolute(dB_value, reference_value):
    return reference_value * 10**(dB_value/10)


def mca_h1cdf_dB_to_absolute(spec_type: str):
    prefix = 'akb_mca_'
    if spec_type == 'pwr':
        tvar_names = [prefix+'Emax', prefix+'Eave', prefix+'Bmax', prefix+'Bave']
        for i in range(4):
            tvar = get_data(tvar_names[i])
            tvar_pwr = dB_to_absolute((tvar.y).astype(float), 1e-12)
            # (mV/m)^2/Hz or pT^2/Hz
            store_data(tvar_names[i] + '_pwr',
                       data={'x': tvar.times, 'y': tvar_pwr, 'v': tvar.v})
            if tvar_names[i] == prefix + 'Emax' or tvar_names[i] == prefix + 'Eave':
                opt_dict = {'spec': 1, 'ylog': 1, 'zlog': 1,
                            'yrange': [1, 2e4], 'ysubtitle': 'Frequency [Hz]',
                            'zrange': [1e-8, 1e2], 'ztitle': '$[(mV/m)^2/Hz]$'}
                options(tvar_names[i] + '_pwr', opt_dict=opt_dict)
            if tvar_names[i] == prefix + 'Bmax' or tvar_names[i] == prefix + 'Bave':
                opt_dict = {'spec': 1, 'ylog': 1, 'zlog': 1,
                            'yrange': [1, 2e4], 'ysubtitle': 'Frecuency [Hz]',
                            'zrange': [1e-5, 1e9], 'ztitle': '$[pT^2/Hz]$'}
                options(tvar_names[i] + '_pwr', opt_dict=opt_dict)

    if spec_type == 'amp':
        tvar_names = [prefix+'Emax', prefix+'Eave', prefix+'Bmax', prefix+'Bave']
        for i in range(4):
            tvar = get_data(tvar_names[i])
            tvar_pwr = dB_to_absolute((tvar.y).astype(float), 1e-12)
            tvar_amp = np.sqrt(tvar_pwr)
            # mV/m/Hz^0.5 or pT/Hz^0.5
            store_data(tvar_names[i] + '_amp',
                       data={'x': tvar.times, 'y': tvar_amp, 'v': tvar.v})
            if tvar_names[i] == prefix + 'Emax' or tvar_names[i] == prefix + 'Eave':
                opt_dict = {'spec': 1, 'ylog': 1, 'zlog': 1,
                            'yrange': [1, 2e4], 'ysubtitle': 'Frequency [Hz]',
                            'zrange': [1e-4, 10], 'ztitle': '$[mV/m/Hz^{0.5}]$'}
                options(tvar_names[i] + '_amp', opt_dict=opt_dict)
            if tvar_names[i] == prefix + 'Bmax' or tvar_names[i] == prefix + 'Bave':
                opt_dict = {'spec': 1, 'ylog': 1, 'zlog': 1,
                            'yrange': [1, 2e4], 'ysubtitle': 'Frequency [Hz]',
                            'zrange': [1e-4, 1e3], 'ztitle': '$[pT/Hz^{0.5}]$'}
                options(tvar_names[i] + '_amp', opt_dict=opt_dict)
