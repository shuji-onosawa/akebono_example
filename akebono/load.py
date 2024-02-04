import logging
from pyspedas.utilities.dailynames import dailynames
from pyspedas.utilities.download import download
from pyspedas.analysis.time_clip import time_clip as tclip
from pytplot import cdf_to_tplot
import os
from .config import CONFIG
import urllib.request


def load(trange=['2012-10-01', '2012-10-02'],
         instrument='pws',
         datatype='epd',
         level='l2',
         suffix='',
         get_support_data=False,
         varformat=None,
         varnames=[],
         downloadonly=False,
         notplot=False,
         no_update=False,
         time_clip=False):
    """
    This function loads data from the Akebono mission; this function is not meant 
    to be called directly; instead, see the wrappers:

        pyspedas.akebono.pws
        pyspedas.akebono.rdm
        pyspedas.akebono.orb

    """
    prefix = ''

    #  PWS and MCA data are available in CDF files
    if instrument == 'pws':
        prefix = 'akb_pws_'
        pathformat = instrument + '/NPW-DS/%Y/ak_h1_pws_%Y%m%d_v??.cdf'
    elif instrument == 'rdm':
        prefix = 'akb_rdm_'
        pathformat = instrument + '/%Y/sf%y%m%d'
    elif instrument == 'orb':
        prefix = 'akb_orb_'
        pathformat = 'orbit/daily/%Y%m/ED%y%m%d.txt'
    elif instrument == 'mca':
        prefix = 'akb_mca_'
        pathformat = 'https://akebono-vlf.db.kanazawa-u.ac.jp/permalink.php?keyword=ak_h1_mca_%Y%m%d_v02.cdf'
    else:
        logging.error('Unknown instrument: ' + instrument)
        return

    # find the full remote path names using the trange
    remote_names = dailynames(file_format=pathformat, trange=trange)
    out_files = []

    # fiels are local file names
    if instrument != 'mca':
        files = download(remote_file=remote_names, remote_path=CONFIG['remote_data_dir'], local_path=CONFIG['local_data_dir'], no_download=no_update)
    if instrument == 'mca':
        local_path = CONFIG['local_data_dir']+'vlf/mca/h1/ave8s/'+trange[0][0:4]+'/'
        files = download_mca_files(remote_names, local_path)

    if files is not None:
        for file in files:
            out_files.append(file)

    out_files = sorted(out_files)
    if downloadonly or (instrument != 'pws' and instrument != 'mca'):
        return out_files

    tvars = cdf_to_tplot(out_files, prefix=prefix, suffix=suffix, get_support_data=get_support_data, varformat=varformat, varnames=varnames, notplot=notplot)
    
    if notplot:
        return tvars

    if time_clip:
        for new_var in tvars:
            tclip(new_var, trange[0], trange[1], suffix='')

    return tvars


def download_mca_files(remote_names, local_path):
    # make an empty list of local file names
    local_names = []
    for remote_name in remote_names:
        # example
        # remote_name = 'https://akebono-vlf.db.kanazawa-u.ac.jp/permalink.php?keyword=ak_h1_mca_20121001_v02.cdf'
        # local_name = 'vlf/mca/h1/ave8s/ak_h1_mca_20121001_v02.cdf'
        local_name = local_path + remote_name.split('=')[-1]
        # If local_path does not exist, create it
        if not os.path.exists(local_path):
            os.makedirs(local_path)
        
        # If local_name does not exist, download
        if not os.path.exists(local_name):
            logging.info('Downloading remote index: ' + remote_name)
            logging.info("Downloading remote file: " + remote_name + " to "+ local_name)
            data = urllib.request.urlopen(remote_name).read()
            with open(local_name, mode="wb") as f:
                f.write(data)            
        if os.path.exists(local_name):
            logging.info("Local file: " + local_name + " exists.")
        
        local_names.append(local_name)

    return local_names

