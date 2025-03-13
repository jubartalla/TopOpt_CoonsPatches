# -*- coding: utf-8 -*-
"""
Created on Fri Mar 7 09:05:26 2025

@author: Ahmed Mohamed Jubartalla Ali

FUNCTIONS: 
surf_squr       Evaluates Coons patch inside a square surface
surf_tria       Evaluates Coons patch inside a triangle surface
volm_cube       Evaluates Coons patch inside a cube volume
volm_trip       Evaluates Coons patch inside a trianglular prism volume

INPUT:
var_xxx, var_yyy, var_zzz   Cartesian coordinates
crnr_ij_arr                 Elastcity tensor for the material at x=i & y=j
curv_xj_arr, derv_xj_arr    Interpolated property and its slope along x-axis at y=j
curv_iy_arr, derv_iy_arr    Interpolated property and its slope along y-axis at x=i

crnr_ijk_arr                Elastcity tensor for the material at x=i & y=j & z=k
curv_xjk_arr, derv_xjk_arr  Interpolated property and its slope along x-axis at y=j & z=k
curv_iyk_arr, derv_1yk_arr  Interpolated property and its slope along y-axis at x=i & z=k
curv_ijz_arr, derv_ijz_arr  Interpolated property and its slope along z-axis at x=i & y=j

var_uuu, var_vvv                        Barycentric coordinates
crnr_u_arr, crnr_v_arr, crnr_w_arr      Elastcity tensor for the material at u, v, w
curv_kl_arr, derv_kl_arr                Interpolated property and its slope in the direction from node_k to node_l

crnr_ui_arr, crnr_vi_arr, crnr_wi_arr   Elastcity tensor for the material at u, v, w at z=i
curv_kli_arr, derv_kli_arr              Interpolated property and its slope in the direction from node_k to node_l at z=i
"""


def proj_var(var_uuu, var_vvv):
    #
    var_wu = var_uuu / (1 - var_vvv)
    var_wv = var_vvv / (1 - var_uuu)
    var_uv = var_vvv / (var_uuu + var_vvv)
    
    dvar_wu_du = 1 / (1 - var_vvv)
    dvar_wv_dv = 1 / (1 - var_uuu)
    
    dvar_wu_dv = var_uuu / (1 - var_vvv)**2
    dvar_wv_du = var_vvv / (1 - var_uuu)**2
    
    dvar_uv_du = -var_vvv / (var_uuu + var_vvv)**2
    dvar_uv_dv =  var_uuu / (var_uuu + var_vvv)**2
    
    #
    list_wu = [var_wu, dvar_wu_du, dvar_wu_dv]
    list_wv = [var_wv, dvar_wv_du, dvar_wv_dv]
    list_uv = [var_uv, dvar_uv_du, dvar_uv_dv]
    
    return list_wu, list_wv, list_uv


def map_drv(drvu_arr, drvv_arr):
    # read
    pnt_x_vec = lib.domn.dmn_x_list.copy()
    pnt_y_vec = lib.domn.dmn_y_list.copy()
    
    #    
    du_dx_val = 1 / (pnt_x_vec[0] - pnt_x_vec[2])
    du_dy_val = 1 / (pnt_y_vec[0] - pnt_y_vec[2])
    
    dv_dx_val = 1 / (pnt_x_vec[1] - pnt_x_vec[2])
    dv_dy_val = 1 / (pnt_y_vec[1] - pnt_y_vec[2])    
        
    # map
    drvx_arr = drvu_arr * du_dx_val + drvv_arr * dv_dx_val
    drvy_arr = drvu_arr * du_dy_val + drvv_arr * dv_dy_val
    
    # return
    return drvx_arr, drvy_arr


def surf_squr(var_list, crnr_arr_list, curv_arr_list, derv_arr_list):    
    # read
    var_xxx, var_yyy = var_list
    
    crnr_00_arr, crnr_10_arr, crnr_11_arr, crnr_01_arr = crnr_arr_list
    curv_x0_arr, curv_1y_arr, curv_x1_arr, curv_0y_arr = curv_arr_list
    derv_x0_arr, derv_1y_arr, derv_x1_arr, derv_0y_arr = derv_arr_list
    
    # ruled surfaces
    surf_x0_arr = (1 - var_yyy) * curv_x0_arr + var_yyy * curv_x1_arr
    surf_0y_arr = (1 - var_xxx) * curv_0y_arr + var_xxx * curv_1y_arr
    
    surf_xy_arr = (1 - var_xxx) * (1 - var_yyy) * crnr_00_arr + var_xxx * (1 - var_yyy) * crnr_10_arr + var_xxx * var_yyy * crnr_11_arr + (1 - var_xxx) * var_yyy * crnr_01_arr
    
    # derivatives
    drvx_0y_arr = curv_1y_arr - curv_0y_arr
    drvy_x0_arr = curv_x1_arr - curv_x0_arr
    
    drvx_x0_arr = (1 - var_yyy) * derv_x0_arr + var_yyy * derv_x1_arr
    drvy_0y_arr = (1 - var_xxx) * derv_0y_arr + var_xxx * derv_1y_arr
    
    drvx_xy_arr = -(1 - var_yyy) * crnr_00_arr + (1 - var_yyy) * crnr_10_arr + var_yyy * crnr_11_arr - var_yyy  * crnr_01_arr
    drvy_xy_arr = -(1 - var_xxx) * crnr_00_arr + (1 - var_xxx) * crnr_01_arr + var_xxx * crnr_11_arr - var_xxx  * crnr_10_arr
    
    # combine
    surf_arr = surf_x0_arr + surf_0y_arr - surf_xy_arr
    drvx_arr = drvx_x0_arr + drvx_0y_arr - drvx_xy_arr
    drvy_arr = drvy_x0_arr + drvy_0y_arr - drvy_xy_arr
    
    # return
    return surf_arr, drvx_arr, drvy_arr


def surf_tria(var_uuu, var_vvv, crnr_arr_list, curv_arr_list, derv_arr_list):
    # read
    crnr_u_arr, crnr_v_arr, crnr_w_arr = crnr_arr_list
    curv_wu_arr, curv_uv_arr, curv_vu_arr, curv_wv_arr, curv_vw_arr, curv_uw_arr = curv_arr_list
    derv_wu_arr, derv_uv_arr, derv_vu_arr, derv_wv_arr, derv_vw_arr, derv_uw_arr = derv_arr_list
       
    # normailize
    list_wu, list_wv, list_uv = proj_var(var_uuu, var_vvv)
    
    var_wu, dvar_wu_du, dvar_wu_dv = list_wu
    var_wv, dvar_wv_du, dvar_wv_dv = list_wv
    var_uv, dvar_uv_du, dvar_uv_dv = list_uv
        
    # ruled surfaces    
    surf_u_arr = (1 - var_wu) * curv_wv_arr + var_wu * curv_uv_arr
    surf_v_arr = (1 - var_wv) * curv_wu_arr + var_wv * curv_vu_arr
    surf_w_arr = (1 - var_uv) * curv_uw_arr + var_uv * curv_vw_arr
    
    surf_0_arr = var_uuu * crnr_u_arr + var_vvv * crnr_v_arr + (1 - var_uuu - var_vvv) * crnr_w_arr
    
    # derivatives
    drvu_u_arr = -dvar_wu_du * curv_wv_arr + dvar_wu_du * curv_uv_arr
    drvv_v_arr = -dvar_wv_dv * curv_wu_arr + dvar_wv_dv * curv_vu_arr
    
    drvv_u_arr = -(1 - var_wu) * derv_wv_arr - dvar_wu_dv * curv_wv_arr + var_wu * derv_uv_arr + dvar_wu_dv * curv_uv_arr
    drvu_v_arr =  (1 - var_wv) * derv_wu_arr - dvar_wv_du * curv_wu_arr - var_wv * derv_vu_arr + dvar_wv_du * curv_vu_arr
	
    drvu_w_arr =  (1 - var_uv) * derv_uw_arr - dvar_uv_du * curv_uw_arr - var_uv * derv_vw_arr + dvar_uv_du * curv_vw_arr
    drvv_w_arr =  (1 - var_uv) * derv_uw_arr - dvar_uv_dv * curv_uw_arr - var_uv * derv_vw_arr + dvar_uv_dv * curv_vw_arr

    drvu_0_arr = crnr_u_arr - crnr_w_arr
    drvv_0_arr = crnr_v_arr - crnr_w_arr
    
    # combine
    surf_arr = surf_u_arr + surf_v_arr + surf_w_arr - surf_0_arr
    drvu_arr = drvu_u_arr + drvu_v_arr + drvu_w_arr - drvu_0_arr
    drvv_arr = drvv_u_arr + drvv_v_arr + drvv_w_arr - drvv_0_arr
    
    # correct
    surf_arr *= 1/2
    drvu_arr *= 1/2
    drvv_arr *= 1/2
    
    # return
    return surf_arr, drvu_arr, drvv_arr


def volm_cube(var_list, crnr_arr_list, curv_arr_list_list, derv_arr_list_list):    
    # read
    x_val, y_val, z_val = var_list
    
    crvx_arr_list, crvy_arr_list, crvz_arr_list = curv_arr_list_list
    drvx_arr_list, drvy_arr_list, drvz_arr_list = derv_arr_list_list
    
    crnr_iii_arr, crnr_jii_arr, crnr_jji_arr, crnr_iji_arr, crnr_iij_arr, crnr_jij_arr, crnr_jjj_arr, crnr_ijj_arr = crnr_arr_list
    
    curv_xii_arr, curv_xji_arr, curv_xjj_arr, curv_xij_arr = crvx_arr_list
    curv_iyi_arr, curv_jyi_arr, curv_jyj_arr, curv_iyj_arr = crvy_arr_list
    curv_iiz_arr, curv_jiz_arr, curv_jjz_arr, curv_ijz_arr = crvz_arr_list
    
    derv_xii_arr, derv_xji_arr, derv_xjj_arr, derv_xij_arr = drvx_arr_list
    derv_iyi_arr, derv_jyi_arr, derv_jyj_arr, derv_iyj_arr = drvy_arr_list
    derv_iiz_arr, derv_jiz_arr, derv_jjz_arr, derv_ijz_arr = drvz_arr_list
    
    # prepare
    var_xii = (1 - y_val) * (1 - z_val)
    var_xij = (1 - y_val) * (    z_val)
    var_xji = (    y_val) * (1 - z_val)
    var_xjj = (    y_val) * (    z_val)
    
    dry_xii = -(1 - z_val)
    dry_xij = -(    z_val)
    dry_xji =  (1 - z_val)
    dry_xjj =  (    z_val)
    
    drz_xii = -(1 - y_val)
    drz_xij =  (1 - y_val)
    drz_xji = -(    y_val)
    drz_xjj =  (    y_val)
    
    var_iyi = (1 - x_val) * (1 - z_val)
    var_iyj = (1 - x_val) * (    z_val)
    var_jyi = (    x_val) * (1 - z_val)
    var_jyj = (    x_val) * (    z_val)
    
    drx_iyi = -(1 - z_val)
    drx_iyj = -(    z_val)
    drx_jyi =  (1 - z_val)
    drx_jyj =  (    z_val)
    
    drz_iyi = -(1 - x_val)
    drz_iyj =  (1 - x_val)
    drz_jyi = -(    x_val)
    drz_jyj =  (    x_val)
    
    var_iiz = (1 - x_val) * (1 - y_val)
    var_ijz = (1 - x_val) * (    y_val)
    var_jiz = (    x_val) * (1 - y_val)
    var_jjz = (    x_val) * (    y_val)
    
    drx_iiz = -(1 - y_val)
    drx_ijz = -(    y_val)
    drx_jiz =  (1 - y_val)
    drx_jjz =  (    y_val)
    
    dry_iiz = -(1 - x_val)
    dry_ijz =  (1 - x_val)
    dry_jiz = -(    x_val)
    dry_jjz =  (    x_val)
    
    var_iii = (1 - x_val) * (1 - y_val) * (1 - z_val)
    var_iij = (1 - x_val) * (1 - y_val) * (    z_val)
    var_iji = (1 - x_val) * (    y_val) * (1 - z_val)
    var_ijj = (1 - x_val) * (    y_val) * (    z_val)
    var_jii = (    x_val) * (1 - y_val) * (1 - z_val)
    var_jij = (    x_val) * (1 - y_val) * (    z_val)
    var_jji = (    x_val) * (    y_val) * (1 - z_val)
    var_jjj = (    x_val) * (    y_val) * (    z_val)
    
    drx_iii = -(1 - y_val) * (1 - z_val)
    drx_iij = -(1 - y_val) * (    z_val)
    drx_iji = -(    y_val) * (1 - z_val)
    drx_ijj = -(    y_val) * (    z_val)
    drx_jii =  (1 - y_val) * (1 - z_val)
    drx_jij =  (1 - y_val) * (    z_val)
    drx_jji =  (    y_val) * (1 - z_val)
    drx_jjj =  (    y_val) * (    z_val)
    
    dry_iii = -(1 - x_val) * (1 - z_val)
    dry_iij = -(1 - x_val) * (    z_val)
    dry_iji =  (1 - x_val) * (1 - z_val)
    dry_ijj =  (1 - x_val) * (    z_val)
    dry_jii = -(    x_val) * (1 - z_val)
    dry_jij = -(    x_val) * (    z_val)
    dry_jji =  (    x_val) * (1 - z_val)
    dry_jjj =  (    x_val) * (    z_val)
    
    drz_iii = -(1 - x_val) * (1 - y_val)
    drz_iij =  (1 - x_val) * (1 - y_val)
    drz_iji = -(1 - x_val) * (    y_val)
    drz_ijj =  (1 - x_val) * (    y_val)
    drz_jii = -(    x_val) * (1 - y_val)
    drz_jij =  (    x_val) * (1 - y_val)
    drz_jji = -(    x_val) * (    y_val)
    drz_jjj =  (    x_val) * (    y_val)
    
    # ruled volumes
    volm_xy0_arr = var_iiz * curv_iiz_arr + var_ijz * curv_ijz_arr + var_jiz * curv_jiz_arr + var_jjz * curv_jjz_arr
    volm_0yz_arr = var_xjj * curv_xjj_arr + var_xii * curv_xii_arr + var_xij * curv_xij_arr + var_xji * curv_xji_arr
    volm_x0z_arr = var_iyi * curv_iyi_arr + var_iyj * curv_iyj_arr + var_jyi * curv_jyi_arr + var_jyj * curv_jyj_arr
    
    volm_xyz_arr = (var_iii * crnr_iii_arr + var_iji * crnr_iji_arr + var_ijj * crnr_ijj_arr + var_iij * crnr_iij_arr + 
    var_jii * crnr_jii_arr + var_jji * crnr_jji_arr + var_jjj * crnr_jjj_arr + var_jij * crnr_jij_arr)
    
	# derivatives
    drvz_xy0_arr = var_iiz * derv_iiz_arr + var_ijz * derv_ijz_arr + var_jiz * derv_jiz_arr + var_jjz * derv_jjz_arr
    drvx_0yz_arr = var_xjj * derv_xjj_arr + var_xii * derv_xii_arr + var_xij * derv_xij_arr + var_xji * derv_xji_arr
    drvy_x0z_arr = var_iyi * derv_iyi_arr + var_iyj * derv_iyj_arr + var_jyi * derv_jyi_arr + var_jyj * derv_jyj_arr
    
    drvx_xy0_arr = drx_iiz * curv_iiz_arr + drx_ijz * curv_ijz_arr + drx_jiz * curv_jiz_arr + drx_jjz * curv_jjz_arr
    drvy_xy0_arr = dry_iiz * curv_iiz_arr + dry_ijz * curv_ijz_arr + dry_jiz * curv_jiz_arr + dry_jjz * curv_jjz_arr
    
    drvy_0yz_arr = dry_xjj * curv_xjj_arr + dry_xii * curv_xii_arr + dry_xij * curv_xij_arr + dry_xji * curv_xji_arr
    drvz_0yz_arr = drz_xjj * curv_xjj_arr + drz_xii * curv_xii_arr + drz_xij * curv_xij_arr + drz_xji * curv_xji_arr
    
    drvx_x0z_arr = drx_iyi * curv_iyi_arr + drx_iyj * curv_iyj_arr + drx_jyi * curv_jyi_arr + drx_jyj * curv_jyj_arr
    drvz_x0z_arr = drz_iyi * curv_iyi_arr + drz_iyj * curv_iyj_arr + drz_jyi * curv_jyi_arr + drz_jyj * curv_jyj_arr
    
    drvx_xyz_arr = (  drx_iii * crnr_iii_arr + drx_iji * crnr_iji_arr + drx_ijj * crnr_ijj_arr + drx_iij * crnr_iij_arr 
                    + drx_jii * crnr_jii_arr + drx_jji * crnr_jji_arr + drx_jjj * crnr_jjj_arr + drx_jij * crnr_jij_arr)
 
    drvy_xyz_arr = (  dry_iii * crnr_iii_arr + dry_iji * crnr_iji_arr + dry_ijj * crnr_ijj_arr + dry_iij * crnr_iij_arr 
                    + dry_jii * crnr_jii_arr + dry_jji * crnr_jji_arr + dry_jjj * crnr_jjj_arr + dry_jij * crnr_jij_arr)
    
    drvz_xyz_arr = (  drz_iii * crnr_iii_arr + drz_iji * crnr_iji_arr + drz_ijj * crnr_ijj_arr + drz_iij * crnr_iij_arr 
                    + drz_jii * crnr_jii_arr + drz_jji * crnr_jji_arr + drz_jjj * crnr_jjj_arr + drz_jij * crnr_jij_arr)
    
    # combine
    volm_arr = volm_xy0_arr + volm_0yz_arr + volm_x0z_arr - 2 * volm_xyz_arr
    drvx_arr = drvx_xy0_arr + drvx_0yz_arr + drvx_x0z_arr - 2 * drvx_xyz_arr
    drvy_arr = drvy_xy0_arr + drvy_0yz_arr + drvy_x0z_arr - 2 * drvy_xyz_arr
    drvz_arr = drvz_xy0_arr + drvz_0yz_arr + drvz_x0z_arr - 2 * drvz_xyz_arr
    
    # return
    return volm_arr, drvx_arr, drvy_arr, drvz_arr

    
def volm_trip(var_list, crnr_arr_list_list, curv_arr_list_list, derv_arr_list_list):    
    # read
    var_uuu, var_vvv, var_zzz = var_list
    
    crni_arr_list, crnj_arr_list = crnr_arr_list_list
    
    crvi_arr_list, crvj_arr_list, crvz_arr_list = curv_arr_list_list
    drvi_arr_list, drvj_arr_list, drvz_arr_list = derv_arr_list_list
    
    crnr_ui_arr, crnr_vi_arr, crnr_wi_arr = crni_arr_list
    crnr_uj_arr, crnr_vj_arr, crnr_wj_arr = crnj_arr_list
    
    curv_wui_arr, curv_uvi_arr, curv_vui_arr, curv_wvi_arr, curv_vwi_arr, curv_uwi_arr = crvi_arr_list
    curv_wuj_arr, curv_uvj_arr, curv_vuj_arr, curv_wvj_arr, curv_vwj_arr, curv_uwj_arr = crvj_arr_list
    
    curv_uu_arr, curv_vv_arr, curv_ww_arr = crvz_arr_list
    derv_uu_arr, derv_vv_arr, derv_ww_arr = drvz_arr_list
    
    derv_wui_arr, derv_uvi_arr, derv_vui_arr, derv_wvi_arr, derv_vwi_arr, derv_uwi_arr = drvi_arr_list
    derv_wuj_arr, derv_uvj_arr, derv_vuj_arr, derv_wvj_arr, derv_vwj_arr, derv_uwj_arr = drvj_arr_list
    
    # surfaces
    surf_t0_arr, drvu_t0_arr, drvv_t0_arr = surf_tria(var_uuu, var_vvv, crni_arr_list, crvi_arr_list, drvi_arr_list)
    surf_t1_arr, drvu_t1_arr, drvv_t1_arr = surf_tria(var_uuu, var_vvv, crnj_arr_list, crvj_arr_list, drvj_arr_list)
    
    # volumes    
    volm_zz_arr = (1 - var_zzz) * surf_t0_arr + var_zzz * surf_t1_arr    
    volm_uv_arr = var_uuu * curv_uu_arr + var_vvv * curv_vv_arr + (1 - var_uuu - var_vvv) * curv_ww_arr
    volm_oo_arr = (1 - var_zzz) * (var_uuu * crnr_ui_arr + var_vvv * crnr_vi_arr + (1 - var_uuu - var_vvv) * crnr_wi_arr) +\
                       var_zzz  * (var_uuu * crnr_uj_arr + var_vvv * crnr_vj_arr + (1 - var_uuu - var_vvv) * crnr_wj_arr)
            
    # derivatives
    drvu_zz_arr = (1 - var_zzz) * drvu_t0_arr + var_zzz * drvu_t1_arr    
    drvv_zz_arr = (1 - var_zzz) * drvv_t0_arr + var_zzz * drvv_t1_arr    
    drvz_zz_arr = -surf_t0_arr + surf_t1_arr
    
    drvu_uv_arr = curv_uu_arr - curv_ww_arr
    drvv_uv_arr = curv_vv_arr - curv_ww_arr
    drvz_uv_arr = var_uuu * derv_uu_arr + var_vvv * derv_vv_arr + (1 - var_uuu - var_vvv) * derv_ww_arr
    
    drvu_oo_arr = (1 - var_zzz) * (crnr_ui_arr - crnr_wi_arr) + var_zzz * (crnr_uj_arr - crnr_wj_arr)
    drvv_oo_arr = (1 - var_zzz) * (crnr_vi_arr - crnr_wi_arr) + var_zzz * (crnr_vj_arr - crnr_wj_arr)
    drvz_oo_arr = -(var_uuu * crnr_ui_arr + var_vvv * crnr_vi_arr + (1 - var_uuu - var_vvv) * crnr_wi_arr) \
                  +(var_uuu * crnr_uj_arr + var_vvv * crnr_vj_arr + (1 - var_uuu - var_vvv) * crnr_wj_arr)
    
    # combine
    volm_arr = volm_zz_arr + volm_uv_arr / 2 - volm_oo_arr / 2
    drvu_arr = drvu_zz_arr + drvu_uv_arr / 2 - drvu_oo_arr / 2
    drvv_arr = drvv_zz_arr + drvv_uv_arr / 2 - drvv_oo_arr / 2
    drvz_arr = drvz_zz_arr + drvz_uv_arr / 2 - drvz_oo_arr / 2
    
    # map
    drvx_arr, drvy_arr = map_drv(drvu_arr, drvv_arr)
             
    # return
    return volm_arr, drvx_arr, drvy_arr, drvz_arr