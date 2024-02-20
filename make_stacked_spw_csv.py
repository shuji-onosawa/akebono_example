# 指定フォルダ内のspwデータのcdfファイルを読み込み、指定したdelaytimeのデータをcsvファイルに書き出す
# csvはdelaytimeごとにファイルを分ける
# 複数delaytimeのデータの場合は、それぞれのdelaytimeのデータを結合したcsvファイルを作成する

import os
import re
import cdflib
import numpy as np

delayTimeIdx = 6 # cdfに保存されているdelaytimeの配列のインデックス
dataDir = 'data/spw/'

timeList = []
heightList = []
gmlatList = []
gmltList = []
fceList = []
rxMatrix = []

cdfList = os.listdir(dataDir) # cdfファイルの名前はak_h1_pws_spw_YYYYMMDD_HHMMSS.cdf
cdfList.sort() # 時刻順にソート
for cdfName in cdfList:
    if re.match(r'.*\.cdf', cdfName): # .cdfファイルのみ
        xry = cdflib.cdf_to_xarray(dataDir+cdfName)
        delayTimeAry = xry['Delay'].values
        rxXrySel = xry['RX'].sel(Delay=delayTimeAry[delayTimeIdx], method='nearest')
        rxAry = rxXrySel.values
        rxMatrix.append(rxAry.tolist())

        # 時刻、高度、地磁気緯度、地磁気地方時、電子サイクロトロン周波数は1行目のみ取得
        epoch = xry['Epoch'].values[0]
        # epoch(s) to numpy datetime64
        epochDatetime = cdflib.cdfepoch.to_datetime(epoch)
        epochStr = epochDatetime[0].strftime('%Y-%m-%d %H:%M:%S')
        timeList.append(epochStr)

        heightList.append(xry['Height'].values[0]) # 高度 (m)
        gmlatList.append(xry['GMLat'].values[0]) # 地磁気緯度 (deg)
        gmltList.append(xry['GMLT'].values[0]) # 地磁気地方時 (hr)
        fceList.append(xry['Fce'].values[0]) # 電子サイクロトロン周波数 (Hz)

# データをcsvに書き出す
# ヘッダー用に周波数のリストをcdfから取得
freqList = xry['Frequency'].values.tolist()
freqList = [str(int(f)) for f in freqList]
# ヘッダーはtime、hight、地磁気緯度、地磁気地方時、電子サイクロトロン周波数、掃引周波数
header = ['time'] + ['height'] + ['gmlat'] + ['gmlt'] + ['fce'] + freqList
# データは時刻、高度、地磁気緯度、地磁気地方時、電子サイクロトロン周波数、RXのdB値
data = np.c_[timeList, heightList, gmlatList, gmltList, fceList, rxMatrix]
data = np.insert(data, 0, header, axis=0)
os.makedirs('execute', exist_ok=True)
np.savetxt('execute/delaytime_'+str(delayTimeAry[delayTimeIdx]*1000)+'ms.csv', data, delimiter=',', fmt='%s')
