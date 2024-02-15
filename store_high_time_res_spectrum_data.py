import pytplot
import akebono
from datetime import datetime

def store_mca_high_time_res_data(date: str = '1990-02-25',
                                 datatype: str = 'pwr',
                                 del_invalid_data: list = ['off', 'bit rate m']):
    """
    高時間分解能のMCAスペクトルデータをtplot変数にする
    Parameters
    ----------
    date : str, optional
    datatype : str, optional
        データタイプ, by default 'pwr'
    del_invalid_data : list of string \n
        mca cdf contain data from which the interference by BDR or SMS is *not* yet removed. \n
        You can remove data contaminated by interference by passing a list containing the following words.\n
        'off': mca is off\n
        'noisy': data is noisy\n
        'sms': SMS is on\n
        'bdr': BDR is on\n
        'bit rate m': Bit rate is medium. When the bit rate is medium, the data is not reliable.\n
        'pws': PWS sounder on\n
    """
    # dateを'yyyymmdd'に変換する
    year = date[:4]
    date = datetime.strptime(date, '%Y-%m-%d').strftime('%Y%m%d')

    mca_cdf_name = '../akebono_data/vlf/mca/h1/ave0.5s/'+year+'/ak_h1_mca_'+date+'_v02.cdf'
    pytplot.cdf_to_tplot(mca_cdf_name, prefix='akb_mca_', get_metadata=True)
    akebono.mca_postprocessing(datatype=datatype, del_invalid_data=del_invalid_data)
