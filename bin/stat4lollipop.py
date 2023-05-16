import glob
import pandas as pd
import os
import shutil

go_filenames = glob.glob("GO*.csv")
for go_filename in go_filenames:
    df = pd.read_csv(go_filename)
    # 将表格按照 "Condition" 列和 "ONTOLOGY" 列分组
    grouped = df.groupby(['Condition', 'ONTOLOGY'])
    # 遍历分组后的每一组
    for name, group in grouped:
        # 获取 Condition 列的值和 ONTOLOGY 列的值
        condition_value, ontology_value = name
        # 创建文件夹
        original_name = condition_value
        if not os.path.exists(original_name):
            os.makedirs(original_name)
        # 生成新的表格文件名
        new_filename = '{}_{}_{}'.format(condition_value, ontology_value,go_filename,)
        # 将分组后的数据写入新的表格文件中
        group.to_csv(new_filename, index = False)
        shutil.move(new_filename, original_name)

kegg_filenames = glob.glob("KEGG*.csv")        
for kegg_filename in kegg_filenames:
    df = pd.read_csv(kegg_filename)
    # 将表格按照 "Condition" 列和 "ONTOLOGY" 列分组
    grouped = df.groupby('Condition')
    # 遍历分组后的每一组
    for name, group in grouped:
        # 获取 Condition 列的值和 ONTOLOGY 列的值
        condition_value = name
        # 创建文件夹
        original_name = condition_value
        if not os.path.exists(original_name):
            os.makedirs(original_name)
        # 生成新的表格文件名
        new_filename = '{}_{}'.format(condition_value,kegg_filename,)
        # 将分组后的数据写入新的表格文件中
        group.to_csv(new_filename, index = False)
        shutil.move(new_filename, original_name)