import pytplot
import numpy as np
import pandas as pd
import pytplot
import akebono
import matplotlib.pyplot as plt
'''
MCAスペクトルデータを読み込み、それぞれの周波数ごとに強度vs.時間のグラフを作成する
各プロットには、周波数ごとの強度(Todo:とその1時間平均値、1 sigmaの値を表示する)
プロットにはpytplotモジュールを使用する
'''


def calcMovingAverage(xarray, startTime, endTime, timeWindow):
    '''
    時間移動平均を計算する
    :param xarray: xarray.DataArray
    :param startTime: str (yyyy-mm-dd HH:MM:SS)
    :param endTime: str (yyyy-mm-dd HH:MM:SS)
    :param timeWindow: froat (min)
    :return: xarray.DataArray
    '''
    # startTimeの30分前からendTimeの30分後までの時間範囲を指定
    aveStartTime = pd.to_datetime(startTime) - pd.Timedelta(minutes=30)
    aveEndTime = pd.to_datetime(endTime) + pd.Timedelta(minutes=30)
    xarray = xarray.sel(time=slice(aveStartTime, aveEndTime))
    # データの時間幅を判別
    time = xarray.coords['time'].values
    timeDelta = pd.to_datetime(time[1]) - pd.to_datetime(time[0])

    # 時間移動平均を計算
    timeWindow = int(timeWindow*60/timeDelta.total_seconds())   # timeWindowをデータ数に変換

    aveXry = xarray.rolling(time=timeWindow, center=True, min_periods=1).mean(skipna=True)
    aveXry = aveXry.sel(time=slice(startTime, endTime))
    return aveXry


def calcMovingStd(xarray, startTime, endTime, timeWindow):
    '''
    時間移動標準偏差を計算する
    :param xarray: xarray.DataArray
    :param startTime: str (yyyy-mm-dd HH:MM:SS)
    :param endTime: str (yyyy-mm-dd HH:MM:SS)
    :param timeWindow: froat (min)
    :return: xarray.DataArray
    '''
    # startTimeの30分前からendTimeの30分後までの時間範囲を指定
    aveStartTime = pd.to_datetime(startTime) - pd.Timedelta(minutes=30)
    aveEndTime = pd.to_datetime(endTime) + pd.Timedelta(minutes=30)
    xarray = xarray.sel(time=slice(aveStartTime, aveEndTime))
    # データの時間幅を判別
    time = xarray.coords['time'].values
    timeDelta = pd.to_datetime(time[1]) - pd.to_datetime(time[0])

    # 時間移動標準偏差を計算
    timeWindow = int(timeWindow*60/timeDelta.total_seconds())   # timeWindowをデータ数に変換

    stdXry = xarray.rolling(time=timeWindow, center=True, min_periods=1).std(skipna=True)
    stdXry = stdXry.sel(time=slice(startTime, endTime))
    return stdXry


def storeEpwrLines(xarray, startTime, endTime):
    '''
    MCA電場データのxarray.DataArrayを受け取り、移動平均と移動標準偏差と併せてpytplotに格納する
    :param xarray: xarray.DataArray
    :param startTime: str (yyyy-mm-dd HH:MM:SS)
    :param endTime: str (yyyy-mm-dd HH:MM:SS)
    :return: None
    '''
    # 移動平均を計算
    maXry = calcMovingAverage(xarray, startTime, endTime, 60)
    maMatrix = maXry.values
    '''
    # 移動標準偏差を計算
    print('calc moving std')
    stdXry = calcMovingStd(xarray, startTime, endTime, 60)
    stdMatrix = stdXry.values
    '''
    # 時間範囲を指定
    xarray = xarray.sel(time=slice(startTime, endTime))
    pwrMatrix = xarray.values
    timePwr = xarray.coords['time'].values
    specBin = xarray.coords['spec_bins'].values
    # 各周波数の最大値、最小値を取得
    pwrMax = np.nanmax(pwrMatrix, axis=0)
    pwrMin = np.nanmin(pwrMatrix, axis=0)
    for i in range(pwrMatrix.shape[1]):
        storeName = 'akb_mca_Epwr_ch{}'.format(i+1)
        chPwrMatrix = np.array([pwrMatrix[:, i].T, maMatrix[:, i].T])
        pytplot.store_data(storeName, data={'x': timePwr, 'y': chPwrMatrix.T})
        pytplot.options(storeName, 'color', ['black', 'red'])
        pytplot.options(storeName, 'ytitle', '{} Hz'.format(str(specBin[i])))
        pytplot.options(storeName, 'ysubtitle', '[(mV/m^2)/Hz]')
        pytplot.options(storeName, 'yrange', [pwrMin[i], pwrMax[i]])
        pytplot.options(storeName, 'ylog', 1)
        pytplot.options(storeName, 'legend_names', ['pwr', 'ma'])


def storeBpwrLines(xarray, startTime, endTime):
    '''
    MCA磁場データのxarray.DataArrayを受け取り、移動平均と移動標準偏差と併せてpytplotに格納する
    :param xarray: xarray.DataArray
    :param startTime: str (yyyy-mm-dd HH:MM:SS)
    :param endTime: str (yyyy-mm-dd HH:MM:SS)
    :return: None
    '''
    # 移動平均を計算
    maXry = calcMovingAverage(xarray, startTime, endTime, 60)
    maMatrix = maXry.values
    '''
    # 移動標準偏差を計算
    stdXry = calcMovingStd(xarray, startTime, endTime, 60)
    stdMatrix = stdXry.values
    '''
    # 時間範囲を指定
    xarray = xarray.sel(time=slice(startTime, endTime))
    pwrMatrix = xarray.values
    timePwr = xarray.coords['time'].values
    specBin = xarray.coords['spec_bins'].values
    # 各周波数の最大値、最小値を取得
    pwrMax = np.nanmax(pwrMatrix, axis=0)
    pwrMin = np.nanmin(pwrMatrix, axis=0)
    for i in range(pwrMatrix.shape[1]):
        storeName = 'akb_mca_Bpwr_ch{}'.format(i+1)
        chPwrMatrix = np.array([pwrMatrix[:, i].T, maMatrix[:, i].T])
        pytplot.store_data(storeName, data={'x': timePwr, 'y': chPwrMatrix.T})
        pytplot.options(storeName, 'color', ['black', 'red'])
        pytplot.options(storeName, 'ytitle', '{} Hz'.format(str(specBin[i])))
        pytplot.options(storeName, 'ysubtitle', '[(nT^2/Hz)]')
        pytplot.options(storeName, 'yrange', [pwrMin[i], pwrMax[i]])
        pytplot.options(storeName, 'ylog', 1)
        pytplot.options(storeName, 'legend_names', ['pwr', 'ma'])


def plotEpwrLines(startTime, endTime):
    '''
    MCA電場データのpytplotを受け取り、それぞれの周波数ごとに強度vs.時間のグラフを作成する
    :param startTime: str (yyyy-mm-dd HH:MM:SS)
    :param endTime: str (yyyy-mm-dd HH:MM:SS)
    :return: None
    '''
    # mca epwrのデータを取得
    pwrList = []
    for ch in range(16):
        pwrXry = pytplot.get_data('akb_mca_Epwr_ch{}'.format(ch+1), xarray=True)
        pwrXrySel = pwrXry.sel(time=slice(startTime, endTime))
        pwrList.append([[pwrXrySel['time'].values], [pwrXrySel.values[:,0]], [pwrXrySel.values[:,1]]])
        # pwrList = [time, pwr, ma]
    # angle_b0_Eyのデータを取得
    angleXry = pytplot.get_data('angle_b0_Ey', xarray=True)
    angleXrySel = angleXry.sel(time=slice(startTime, endTime))
    angleAry = angleXrySel.values

    # plot
    # 2*6のサブプロットを作成
    fig, ax = plt.subplots(6, 2, figsize=(16, 18))
    ax = ax.flatten()
    for i in range(12):
        # データを取得
        time = pwrList[i][0]
        pwr = pwrList[i][1]
        MovingAve = pwrList[i][2]
        # MAplus2Sigma = pwrList[i][3]
        # プロット
        ax[i].scatter(time, pwr, s=1, c=angleAry, cmap='jet')
        ax[i].plot(time, MovingAve, color='k', label='ma', linestyle='dashed')
        # ax[i].plot(time, MAplus2Sigma, color='k', label='ma+2sigma', linestyle='dotted')
        ax[i].set_yscale('log')
        ax[i].set_xlabel('time')
        ax[i].set_ylabel('Ch{} \n [(mV/m)^2/Hz]'.format(i+1))
    plt.legend()
    plt.tight_layout()
    fig.colorbar(ax[0].collections[0], ax=ax, orientation='horizontal', label='angle_b0_Ey')
    plt.savefig('test.png')