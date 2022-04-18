### 处理有自旋的情况
import numpy as np
import pandas as pd
from pymatgen.io.vasp.outputs import Vasprun
fermi_energy = 7.8847


def GetDictBandData():
    '''
    读取计算结果中的能带数据，vasprun.xml文件和计算能带所使用的k点文件是必须的
    '''
    my_vasprun = Vasprun(filename="vasprun.xml", separate_spins=True)
    my_bd_data = my_vasprun.get_band_structure(kpoints_filename='KPOINTS')
    dt_band_data = my_bd_data.as_dict()
    return dt_band_data

dt_band_data = GetDictBandData()

### 获取一些重要计算能带结果中的一些重要数据
my_dt_bands_up = dt_band_data["bands"]['1']
my_dt_bands_dn = dt_band_data["bands"]['-1']
my_dt_kpnts = dt_band_data["kpoints"]
my_dt_brchs = dt_band_data["branches"]
my_dt_reclat = dt_band_data["lattice_rec"]['matrix']


kpts_count = len(my_dt_kpnts)
bds_count = len(my_dt_bands)


def GetKptDistance(ls_kpt1, ls_kpt2, rec_lat):
    '''
    根据两个k点的坐标得到它们在倒空间中的距离
    '''
    mat_rec_lat = np.array(rec_lat)
    vec_k1 = ls_kpt1[0]*mat_rec_lat[0]+ls_kpt1[1]*mat_rec_lat[1]+ls_kpt1[2]*mat_rec_lat[2]
    vec_k2 = ls_kpt2[0]*mat_rec_lat[0]+ls_kpt2[1]*mat_rec_lat[1]+ls_kpt2[2]*mat_rec_lat[2]
    k_distance = np.linalg.norm(vec_k2-vec_k1)
    return k_distance


def GetLabelPosiInAxisX(ls_brchs_info, rec_lat, ls_kpnts):
    '''
    获取每个k点在X轴上的位置，即是将三维k点坐标转化到一维空间坐标系
    '''
    brchs_count = len(ls_brchs_info)
    sum_distance = 0.0
    df_klabels = pd.DataFrame(data=ls_brchs_info)
    ls_labels_distance = []
    ls_kpnts_distance = []
    for brchs_index in range(brchs_count):
        dict_brch = ls_brchs_info[brchs_index]
        st_index = dict_brch['start_index']
        ed_index = dict_brch['end_index']
        brch_distance = GetKptDistance(ls_kpt1=ls_kpnts[st_index], ls_kpt2=ls_kpnts[ed_index], rec_lat=rec_lat)
        ls_labels_distance.append(sum_distance+brch_distance)
        for i in range(st_index, ed_index+1):
            kpt_distance = sum_distance+(i-st_index)*brch_distance/(ed_index-st_index)
            ls_kpnts_distance.append(kpt_distance)
        sum_distance = sum_distance+brch_distance
    df_klabels = df_klabels.assign(label_sum_distance=ls_labels_distance)
    df_klabels.to_csv("klable.csv")           ### klabel.csv文件存储了每个k点的符号和在X轴上的位置
    return ls_kpnts_distance

     
def WriteBandsEnergyToTxt(filename, dt_bands_up, dt_bands_dn, ls_kpts, fermi_energy):
    '''
    把能带的能量数据写入到文本呢文件中,其中filename是文件名,dt_bands_up和dt_bands_dn是读取到的自旋向上和向下的能带文件，fermi_energy是费米能
    '''
    bds_count = len(dt_bands_up)
    kps_count = len(ls_kpts)
    with open(file=filename, mode='w') as bd_file:
        bd_file.write("kpos\tenergy_up(eV)\tenergy_down(ev)\tenergy_with_zero_fermi_energy_up\tenergy_with_zero_fermi_energy_down\n")
        for i in range(bds_count):
            bd_energy_up = dt_bands_up[i]
            bd_energy_dn = dt_bands_dn[i]
            for j in range(kps_count):
                bd_energy_kj_up = bd_energy_up[j]
                bd_energy_kj_zero_efermi_up = bd_energy_kj_up-fermi_energy
                bd_energy_kj_dn = bd_energy_dn[j]
                bd_energy_kj_zero_efermi_dn = bd_energy_kj_dn-fermi_energy
                distance_kj = ls_kpts[j]
                str_oneline = '{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(distance_kj, bd_energy_kj_up,bd_energy_kj_dn, bd_energy_kj_zero_efermi_up, bd_energy_kj_zero_efermi_dn)
                bd_file.write(str_oneline)
            bd_file.write("\n")


ls_kd = GetLabelPosiInAxisX(my_dt_brchs, my_dt_reclat, my_dt_kpnts)
WriteBandsEnergyToTxt(filename='bd.txt', dt_bands_up=my_dt_bands_up, dt_bands_dn=my_dt_bands_dn, ls_kpts=ls_kd, fermi_energy=fermi_energy)
