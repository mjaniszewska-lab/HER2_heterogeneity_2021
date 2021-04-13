#####
# author:  engje, grael, bue,
# date: 2020-04-07
# license: GPLv3
#####

# library
import os
import pandas as pd
import shutil
from mplex_image import analyze
import numpy as np

# const input output
s_path = "/home/groups/graylab_share/OMERO.rdsStore/engje/Data/20200706_NP001-36"

ls_file = ["20201021_NP001-36_ManualPositive"]  #_edge31
s_out = '20201022_NP001-36' #'20201012_NP001-36_edge31'


df_data = pd.DataFrame()

for s_file in ls_file:
    print(f'Loading {s_file}')
    df_tt = pd.read_csv(f'./{s_file}.csv', index_col=0)
    df_data=df_data.append(df_tt, sort=True)

# do we want slide or scene?
df_data['slide'] = [item.split('_')[0] for item in df_data.index]
df_data['scene'] = [item.split('_')[1] for item in df_data.index]
#call slide_scene "slide"
df_data['slide'] = df_data.slide + '_' + df_data.scene
df_data.drop('scene',axis=1,inplace=True)

#########################
#### parameters #########
#########################

ls_endothelial = ['CD31_Ring']
ls_immune = ['CD45_Ring','CD68_Ring','CD3_Ring','CD4_Ring'] 
ls_stromal = []
ls_tumor = ['CK7_Ring','CK19_Ring','CK8_Ring','CK5_Ring','CK14_Ring'] #,'Ecad_Ring'
ls_prolif = ['Ki67_Nuclei']

#tcell/myeloid
ls_tcell = ['CD3_Ring','CD4_Ring']
s_bcell = 'CD20_Ring'
s_myeloid = 'CD68_Ring'
ls_immune_functional = ['PD1_Ring','FoxP3_Nuclei','GRNZB_Nuclei','prolif'] #,
df_data.rename(dict(zip(ls_immune_functional,[item.split('_')[0] for item in ls_immune_functional])),axis=1,inplace=True)

#luminal/basal/mesenchymal
ls_luminal = ['CK19_Ring','CK7_Ring','CK8_Ring'] #
ls_basal = ['CK5_Ring','CK14_Ring'] 
ls_mes = ['Vim_Ring'] #
ls_tumor_functional = ['EGFR_Ring']#['Ecad_Ring'] #
ls_stromal_function =  ['Vim_Ring','aSMA_Ring'] # _Ring,'ColIV_Ring''PDPN_Ring',,'CD44_Ring' (missing from K169 set
ls_tumor_signal = ['pS6RP_Ring','PCNA_Nuclei','Ki67_Nuclei','pRB_Nuclei']#'pRB_Nuclei','LamAC_Nuclei']'H3K27_Nuclei','H3K4_Nuclei',

ls_cellline = []
'''
['JE-TMA-22_scene02', 'JE-TMA-22_scene03', 'JE-TMA-22_scene05', 'JE-TMA-22_scene06',
 'JE-TMA-22_scene08', 'JE-TMA-22_scene09', 'JE-TMA-22_scene10', 'JE-TMA-22_scene11',
 'JE-TMA-22_scene13', 'JE-TMA-22_scene14','JE-TMA-32_scene02', 'JE-TMA-32_scene03', 'JE-TMA-32_scene05', 'JE-TMA-32_scene06',
 'JE-TMA-32_scene08', 'JE-TMA-32_scene09', 'JE-TMA-32_scene10', 'JE-TMA-32_scene11',
 'JE-TMA-32_scene13', 'JE-TMA-32_scene14','JE-TMA-35_scene02', 'JE-TMA-35_scene03', 'JE-TMA-35_scene05', 'JE-TMA-35_scene06',
 'JE-TMA-35_scene08', 'JE-TMA-35_scene09', 'JE-TMA-35_scene10', 'JE-TMA-35_scene11',
 'JE-TMA-35_scene13', 'JE-TMA-35_scene14']
'''
#####################################################################
### I. exclusive: ['endothelial',  'tumor', 'stromal',  'immune'] ###
#define 4 main celltype, and proliferation of those 4 types
######################################################################

#celltpye
#1 endothelial
df_data['endothelial'] = df_data.loc[:,ls_endothelial].any(axis=1)

#3 tumor
ls_exclude =  ls_endothelial 
df_data['tumor'] = df_data.loc[:,ls_tumor].any(axis=1) & ~df_data.loc[:,ls_exclude].any(axis=1) 

#2 immune
ls_exclude = ls_endothelial + ls_tumor
df_data['immune'] = df_data.loc[:,ls_immune].any(axis=1) & ~df_data.loc[:,ls_exclude].any(axis=1)

#4 stromal
ls_exclude = ls_immune + ls_endothelial + ls_tumor
df_data['stromal'] = ~df_data.loc[:,ls_exclude].any(axis=1)
#add celltype
ls_cell_names = ['stromal','endothelial','tumor','immune']
s_type_name = 'celltype'
analyze.add_celltype(df_data, ls_cell_names, s_type_name)

#fix cell lines (all tumor!)
for s_cellline in ls_cellline:
    ls_index = df_data[df_data.slide==s_cellline].index
    df_data.loc[ls_index,'celltype'] = 'tumor'
df_data['immune'] = df_data.loc[:,'celltype'] == 'immune'
df_data['stromal'] = df_data.loc[:,'celltype'] == 'stromal'
df_data['endothelial'] = df_data.loc[:,'celltype'] == 'endothelial'

#proliferation
df_data['prolif'] = df_data.loc[:,ls_prolif].any(axis=1)
df_data['nonprolif'] = ~df_data.loc[:,ls_prolif].any(axis=1)
#add proliferation
ls_cell_names = ['prolif','nonprolif']
s_type_name = 'proliferation'
analyze.add_celltype(df_data, ls_cell_names, s_type_name)




#####################################################################################################################################
##### II.  immune status marker ['CD3_Ring', 'CD68_Ring','FoxP3_Nuclei', 'GRNZB_Nuclei','CD8_Ring','CD4_Ring','PD1_Ring','CD45_Ring']
#####################################################################################################################################

## T cell, B cell or myeloid
df_data['MyeloidImmune'] = df_data.loc[:,[s_myeloid,'immune']].all(axis=1) 
df_data['CD20Bcell'] = df_data.loc[:,[s_bcell,'immune']].all(axis=1) & ~df_data.loc[:,['MyeloidImmune']+ls_tcell].any(axis=1)
df_data['TcellImmune'] = df_data.loc[:,['immune']].all(axis=1) & df_data.loc[:,ls_tcell].any(axis=1) & ~df_data.loc[:,['CD20Bcell','MyeloidImmune']].any(axis=1)
df_data['UnspecifiedImmune'] = df_data.loc[:,'immune'] & ~df_data.loc[:,['CD20Bcell','TcellImmune','MyeloidImmune']].any(axis=1)
## CD4 and CD8 
if df_data.columns.isin(['CD8_Ring','CD4_Ring']).sum()==2:
    print('CD4 AND CD8')
    df_data['CD8Tcell'] = df_data.loc[: ,['CD8_Ring','TcellImmune']].all(axis=1)
    df_data['CD4Tcell'] = df_data.loc[: ,['CD4_Ring','TcellImmune']].all(axis=1) & ~df_data.loc[: ,'CD8Tcell']
    df_data['UnspecifiedTcell'] = df_data.TcellImmune & ~df_data.loc[:,['CD8Tcell','CD4Tcell']].any(axis=1) #if cd4 or 8 then sum = 2
    ## check
    ls_immune = df_data[df_data.loc[:,'TcellImmune']].index.tolist()
    if ((df_data.loc[ls_immune,['CD8Tcell','CD4Tcell','UnspecifiedTcell']].sum(axis=1)!=1)).any():
        print('Error in Tcell cell types')
    ls_immuntype = ['MyeloidImmune','CD20Bcell','UnspecifiedImmune','CD8Tcell','CD4Tcell','UnspecifiedTcell'] #'TcellImmune',
#add Immunetype
ls_cell_names = ls_immuntype
s_type_name = 'ImmuneType'
analyze.add_celltype(df_data, ls_cell_names, s_type_name)

#get rid of unspecfied immune cells (make them stroma)
ls_index = df_data[df_data.ImmuneType.fillna('x').str.contains('Unspecified')].index
df_data.loc[ls_index,'celltype'] = 'stromal'
df_data.loc[ls_index,'ImmuneType'] = np.nan
df_data.loc[ls_index,'stromal'] = True
df_data.loc[ls_index,'immune'] = False


#Immune functional states (don't forget to merge!)
df_func = analyze.combinations(df_data,[item.split('_')[0] for item in ls_immune_functional])
df_data = df_data.merge(df_func,how='left', left_index=True, right_index=True, suffixes = ('_all',''))
#gated combinations: immune type plus fuctional status
ls_gate = ['CD4Tcell', 'CD8Tcell']#sorted(df_data[~df_data.ImmuneType.isna()].loc[:,'ImmuneType'].unique()) make fewer immune types
ls_marker = df_func.columns.tolist()
df_gate_counts = analyze.gated_combinations(df_data,ls_gate,ls_marker)
df_data = df_data.merge(df_gate_counts, how='left', left_index=True, right_index=True,suffixes = ('_all',''))
#add TcellImmune
ls_cell_names = df_gate_counts.columns.tolist()
s_type_name ='TcellImmune'
analyze.add_celltype(df_data, ls_cell_names, s_type_name)

#make non Tcell nas
df_data.loc[df_data.TcellImmune==False,'TcellImmune'] = np.nan
df_data.loc[df_data.TcellImmune==True,'TcellImmune'] = np.nan
########################################
#CellProlif combinations, main cell types and proliferation
######################################

ls_gate =  ['endothelial',  'tumor', 'stromal', 'immune']
ls_combo =['prolif','nonprolif']
df_gate_counts2 = analyze.gated_combinations(df_data,ls_gate,ls_combo)
df_data = df_data.merge(df_gate_counts2, how='left', left_index=True, right_index=True,suffixes = ('_all',''))

#add CellProlif
ls_cell_names = ['endothelial_prolif','endothelial_nonprolif', 'tumor_prolif', 'tumor_nonprolif',
       'stromal_prolif', 'stromal_nonprolif', 'immune_prolif','immune_nonprolif']
s_type_name = 'CellProlif'
analyze.add_celltype(df_data, ls_cell_names, s_type_name)


#########################################
#### III tumor status ####################
###########################################

print('differentiation')
df_data['Lum'] = df_data.loc[:,ls_luminal].any(axis=1) & df_data.tumor
df_data['Bas'] = df_data.loc[:,ls_basal].any(axis=1)  & df_data.tumor
df_data['Mes'] = df_data.loc[:,ls_mes].any(axis=1) & df_data.tumor

print('hormonal status')
df_data['ER'] = df_data.loc[:,['tumor','ER_Nuclei']].all(axis=1)
df_data['HER2'] = df_data.loc[:,['tumor','HER2_Ring']].all(axis=1)
if df_data.columns.isin(['PgR_Nuclei']).any():
    df_data['PR'] = df_data.loc[:,['tumor','PgR_Nuclei']].all(axis=1)

#df_data['HR'] = df_data.loc[:,['ER','PR']].any(axis=1) & df_data.tumor
df_data['ER'] = df_data.loc[:,['ER']].any(axis=1) & df_data.tumor

#Diffenrentiation state combinations (256 total)
#1 luminal/basal/ mesenchymal/ functional!
df_data.rename(dict(zip(ls_tumor_functional,[item.split('_')[0] for item in ls_tumor_functional])),axis=1,inplace=True)
ls_marker = ['Lum','Bas','Mes'] +  [item.split('_')[0] for item in ls_tumor_functional]
df_diff = analyze.combinations(df_data,ls_marker)
df_data = df_data.merge(df_diff,how='left', left_index=True, right_index=True, suffixes = ('_all',''))

#add DiffState
ls_cell_names = df_diff.columns.tolist()
s_type_name = 'DiffState'
analyze.add_celltype(df_data, ls_cell_names, s_type_name)
#change non-tumor to NA (works!)
df_data.loc[df_data[df_data.celltype != 'tumor'].index,s_type_name] = np.nan

#2 ER/PR/HER2
ls_marker =  ['ER','HER2']#,['ER','HER2','PR']
df_hr = analyze.combinations(df_data,ls_marker)
df_hr.rename({'__':'TN'},axis=1,inplace=True)
df_data = df_data.merge(df_hr,how='left', left_index=True, right_index=True,suffixes = ('_all',''))
ls_cell_names = df_hr.columns.tolist()
s_type_name = 'HRStatus'
analyze.add_celltype(df_data, ls_cell_names, s_type_name)
#change non-tumor to NA (works!)
df_data.loc[df_data[df_data.celltype != 'tumor'].index,s_type_name] = np.nan

#3 combinations: differentiation and HR status
ls_gate = df_diff.columns.tolist()
ls_marker = df_hr.columns.tolist()
df_gate_counts = analyze.gated_combinations(df_data,ls_gate,ls_marker)
df_data = df_data.merge(df_gate_counts, how='left', left_index=True, right_index=True,suffixes = ('_all',''))

# make Tumor Diff plus HR Status object column
ls_cell_names =  df_gate_counts.columns.tolist()
s_type_name = 'DiffStateHRStatus'
analyze.add_celltype(df_data, ls_cell_names, s_type_name)
#change non-tumor to NA (works!)
df_data.loc[df_data[df_data.celltype != 'tumor'].index,s_type_name] = np.nan

#tumor signaling and proliferation
#rename
df_data.rename(dict(zip(ls_tumor_signal,[item.split('_')[0] for item in ls_tumor_signal])),axis=1,inplace=True)
ls_marker = [item.split('_')[0] for item in ls_tumor_signal]
df_sign = analyze.combinations(df_data,ls_marker)
df_data = df_data.merge(df_sign,how='left', left_index=True, right_index=True,suffixes = ('_all',''))
ls_cell_names = df_sign.columns.tolist()
s_type_name = 'TumorProlif'
analyze.add_celltype(df_data, ls_cell_names, s_type_name)
#change non-tumor to NA (works!)
df_data.loc[df_data[df_data.celltype != 'tumor'].index,s_type_name] = np.nan

#########################################################################################
################ stromal combos #########################################################
############################################################################################

#careful, run this last because of diff state using same marker names (i.e. Vim_Ring)
df_data.rename(dict(zip(ls_stromal_function,[item.split('_')[0] for item in ls_stromal_function])),axis=1,inplace=True)

#functional states (stromal) (don't forget to merge!)
df_func = analyze.combinations(df_data,[item.split('_')[0] for item in ls_stromal_function])
df_data = df_data.merge(df_func,how='left', left_index=True, right_index=True, suffixes = ('_all',''))

#gated combinations: immune type plus fuctional status
ls_gate = ['stromal']
ls_marker = df_func.columns.tolist()
df_gate_counts = analyze.gated_combinations(df_data,ls_gate,ls_marker)
df_data = df_data.merge(df_gate_counts, how='left', left_index=True, right_index=True,suffixes = ('_all',''))

# make Functional Stromal Status object column
ls_cell_names = df_gate_counts.columns.tolist()
s_type_name = 'StromalType'
analyze.add_celltype(df_data, ls_cell_names, s_type_name)
df_data.loc[df_data[df_data.celltype != 'stromal'].index,s_type_name] = np.nan

## Lineage Stroma Definitions
df_data['Pericyte'] = df_data.loc[:,['stromal']].all(axis=1)& df_data.loc[:,['aSMA']].any(axis=1) & ~df_data.loc[:,['Vim']].any(axis=1)
df_data['ActivatedFB'] = df_data.loc[:,['stromal']].all(axis=1)& df_data.loc[:,['Vim','aSMA']].any(axis=1) & ~df_data.loc[:,['Pericyte']].any(axis=1)
# aSMA
df_data['Other_NonTumor']= df_data.loc[:,['stromal']].all(axis=1) & ~df_data.loc[:,['aSMA','Vim']].any(axis=1)

#add Strmaltype
ls_cell_names_str = ['Pericyte','ActivatedFB','Other_NonTumor']
s_type_name = 'MainStromal'
analyze.add_celltype(df_data, ls_cell_names_str, s_type_name)

#####################################################################################################################################

#one more column: all non-tumor cells
index_endothelial = df_data[df_data.celltype=='endothelial'].index
index_immune = df_data[df_data.celltype=='immune'].index
index_stroma = df_data[df_data.celltype=='stromal'].index
index_tumor = df_data[df_data.celltype=='tumor'].index

#more cell types
'''
df_data.loc[index_endothelial,'NonTumorFunc'] = df_data.loc[index_endothelial,'CellProlif']
df_data.loc[index_immune,'NonTumorFunc'] = df_data.loc[index_immune,'FuncImmune']
df_data.loc[index_stroma,'NonTumorFunc'] = df_data.loc[index_stroma,'StromalType']
df_data.loc[index_tumor,'NonTumorFunc'] = np.nan
'''
#fewer cell tpyes
df_data.loc[index_endothelial,'NonTumor'] = 'endothelial'
df_data.loc[index_immune,'NonTumor'] = df_data.loc[index_immune,'ImmuneType']
df_data.loc[index_stroma,'NonTumor'] = df_data.loc[index_stroma,'MainStromal']
df_data.loc[index_tumor,'NonTumor'] = np.nan

###############################################################################################################3
#Final Cell
#add a stromal/immune/endothelial column
for s_marker in ['StromalType','TcellImmune','DiffStateHRStatus']:
    df_data.loc[df_data[~df_data.loc[:,s_marker].isna()].index,'FinalCell'] = df_data[~df_data.loc[:,s_marker].isna()].loc[:,s_marker]

#add immune
for s_cell in  ['CD20Bcell', 'MyeloidImmune']:
    print(s_cell)
    print((df_data.ImmuneType==s_cell).sum())
    df_data.loc[df_data.ImmuneType==s_cell,'FinalCell'] = s_cell
    print((df_data.FinalCell == s_cell).sum())

# add endothelial
df_data['FinalCell'] = df_data.FinalCell.fillna('endothelial')
df_data['GRNZB_immune'] = df_data.loc[:,['GRNZB','immune']].all(axis=1)

#drop rare
ls_drop = df_data.groupby('FinalCell').tumor.count()[df_data.groupby('FinalCell').tumor.count() < 1000].index
if ls_drop.isin(['GRNZB_immune']).sum()!=0:
    print('Drop')
    ls_drop = ls_drop.drop('GRNZB_immune')
#add 
#df_filter_data = df_data.loc[(~(df_data.loc[:,'FinalCell'].isin(ls_drop.tolist()))),:]
ls_keep_index = df_data[~df_data.loc[:,'FinalCell'].isin(ls_drop)].index
df_filter_data = df_data.loc[ls_keep_index,:]

#add back GRNZ to Final cell, immune type, non tumor
df_filter_data.loc[df_filter_data.GRNZB_immune==True,'FinalCell'] = 'GRNZB_immune'
df_filter_data.loc[df_filter_data.GRNZB_immune==True,'ImmuneType'] = 'GRNZB_immune'
df_filter_data.loc[df_filter_data.GRNZB_immune==True,'NonTumor'] = 'GRNZB_immune'

print(sorted(set(df_filter_data.FinalCell)))

#collapse a couple states
df_filter_data.loc[:,'FinalCell'] = df_filter_data.FinalCell.replace({'Bas_EGFR_TN':'Bas_TN','Lum_Bas_EGFR_TN':'Lum_Bas_TN',
 'CD4Tcell___':'CD4Tcell','CD8Tcell___':'CD8Tcell','stromal___':'Other_NonTumor'})
 
print(sorted(set(df_filter_data.FinalCell)))

####################################################################################
############### LAST STEPS #########################################################
####################################################################################

#drop extra colums
df_gate = df_filter_data.loc[:,df_filter_data.dtypes!='bool']
df_gate.to_csv(f'{s_path}/{s_out}_GatedPositiveCellNames.csv')

df_filter_edge = pd.read_csv('20201021_NP001-36_ManualPositive_edge31.csv',index_col=0)
df_gate[df_gate.index.isin(df_filter_edge.index)].to_csv(f'{s_path}/{s_out}_edge31_GatedPositiveCellNames.csv')
