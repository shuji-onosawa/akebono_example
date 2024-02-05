import numpy as np
import xarray as xr
import os
import pytplot
import akebono
import datetime

def get_plot_trange_list(ilat_mlt_ds, mlt_range, ilat_range, mlat_range=[-90, 90]):
    '''
    時刻データの配列とMLT、ILATの範囲を指定して、
    プロットする時刻の範囲を取得する関数
    args:
        ilat_mlt_ds: 時刻, ilat, mltのdataset. 1日分のデータを想定
                    orbit dataの時刻は30秒間隔であることを想定
        mlt_range: プロットするMLTの範囲 [下限値, 上限値] list
        ilat_range: プロットするILATの範囲 [下限値, 上限値] list
    return:
        プロットする時刻の範囲のリスト [[開始時刻, 終了時刻], ...] list
    '''
    # ilat_mlt_dsが'akb_orb_inv'と'akb_orb_MLT'の2つの変数を持っていることを確認
    assert 'akb_orb_inv' in ilat_mlt_ds.variables
    assert 'akb_orb_mlt' in ilat_mlt_ds.variables
    # 条件を指定してデータをフィルタリング
    condition = ((ilat_mlt_ds['akb_orb_inv'] >= ilat_range[0]) &
                 (ilat_mlt_ds['akb_orb_inv'] <= ilat_range[1]) &
                 (ilat_mlt_ds['akb_orb_mlt'] >= mlt_range[0]) &
                 (ilat_mlt_ds['akb_orb_mlt'] <= mlt_range[1]))
    if mlat_range == [-90, 90]:
        pass
    else:
        condition = condition & ((ilat_mlt_ds['akb_orb_mlat'] >= mlat_range[0]) &
                                 (ilat_mlt_ds['akb_orb_mlat'] <= mlat_range[1]))
    filtered_data = ilat_mlt_ds.where(condition, drop=True)
    # filtered_dataの時刻データをnumpy配列として取得
    time_ary = filtered_data['time'].values
    # time_arrayで時間差が30秒以上のインデックスを取得
    break_indices = np.where(np.diff(time_ary) > np.timedelta64(30, 's'))[0] + 1

    # break_indicesが空の場合は、データが存在しないことを示すために空のリストを返す
    if len(break_indices) == 0:
        return []

    # 時間差が30秒以内のブロックの開始時刻と終了時刻のリストを作成
    time_range_list = []
    for i in range(len(break_indices)+1):
        if i == 0:
            start_index = 0
        else:
            start_index = break_indices[i-1]
        if i == len(break_indices):
            end_index = len(time_ary)
        else:
            end_index = break_indices[i]
        time_range_list.append([time_ary[start_index], time_ary[end_index-1]])

    # 各ブロックの開始時刻から4分後の時刻、さらに4分後の時刻...と計算していき、
    # 各ブロックの終了時刻を超えるまでの時刻の範囲を取得
    plot_trange_list = []
    for time_range in time_range_list:
        start_time = time_range[0]
        end_time = time_range[1]
        while True:
            next_time = start_time + np.timedelta64(4, 'm')
            if next_time > end_time:
                break
            plot_trange_list.append([start_time, next_time])
            start_time = next_time

    # plot_time_range_listの中にはnumpy.datetime64型のデータが入っているので、
    # 'YYYY-MM-DD HH:MM:SS'の形式に変換
    for i in range(len(plot_trange_list)):
        plot_trange_list[i][0] = plot_trange_list[i][0].astype('M8[ns]').astype(str)
        plot_trange_list[i][1] = plot_trange_list[i][1].astype('M8[ns]').astype(str)
    return plot_trange_list


def plot_mca_w_1day(date: str):
    '''
    MCAのデータとMGFのデータをプロットする関数
    dateの日のデータをプロットする
    1プロットの時間範囲は4分(general_plot.pyのplot_trange_listの設定に依存)
    '''
    # データを取得、前処理
    ilat_ds = pytplot.get_data('akb_orb_inv', xarray=True)
    mlt_ds = pytplot.get_data('akb_orb_mlt', xarray=True)
    mlat_ds = pytplot.get_data('akb_orb_mlat', xarray=True)
    ilat_mlt_ds = xr.merge([ilat_ds, mlt_ds, mlat_ds])

    # MLTとILATの範囲を指定
    mlt_range = [10, 14]  # MLTの開始値と終了値を指定
    ilat_range = [70, 80]  # ILATの開始値と終了値を指定
    mlat_range = [-90, 90]  # MLATの開始値と終了値を指定

    # プロットする時刻の範囲を取得
    plot_trange_list = get_plot_trange_list(ilat_mlt_ds,
                                            mlt_range, ilat_range, mlat_range)

    # プロット
    # save path の設定
    condition_str = 'mlt_'+str(mlt_range[0])+'_'+str(mlt_range[1])+'_ilat_'+str(ilat_range[0])+'_'+str(ilat_range[1])
    save_dir = '../plots/mca_w_mgf/'+date[:4]+'/'+condition_str+'/'
    os.makedirs(save_dir, exist_ok=True)
    # plot_trange_listの中身が空の場合は、データが存在しないことを表示して終了
    if len(plot_trange_list) == 0:
        print("No data in the specified range")
        return

    for plot_trange in plot_trange_list:
        pytplot.tlimit(plot_trange)
        save_path = save_dir+plot_trange[0][:10]+'_'+plot_trange[0][11:13]+plot_trange[0][14:16]
        pytplot.tplot(['akb_mca_Emax_pwr',
                       'akb_mca_Bmax_pwr'],
                      var_label=['akb_orb_alt', 'akb_orb_inv', 'akb_orb_mlt'],
                      xsize=14, ysize=14, save_png=save_path, display=False)

akebono.orb(['1990-02-11', '1990-02-12'])
akebono.vlf_mca(['1990-02-11', '1990-02-12'],
                datatype='pwr',
                del_invalid_data=['off', 'sms'])
plot_mca_w_1day('1990-02-11')