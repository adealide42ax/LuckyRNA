import glob
import pandas as pd
import re
import sys

filenames = glob.glob('*/3*.csv')

columns = ['GO_diff', 'GO_Up', 'GO_Down', 'KEGG_diff', 'KEGG_Up', 'KEGG_Down']
data = {column: [] for column in columns}
for filename in filenames:
    for category1 in ['GO', 'KEGG']:
        for category2 in ['DESeq2', 'edgeR', 'limma', 'merge']:
            for category3 in ['diff', 'Up', 'Down']:
                if f'{category1}_{category2}_{category3}' in filename:
                    data[f'{category1}_{category3}'].append(filename)
                    break
df = pd.DataFrame(data, columns=columns)

for col in df.columns:
    data_list = []
    for path in df[col]:
        data = pd.read_csv(str(path))
        if 'ONTOLOGY' in data.columns:
            for group, sub_df in data.groupby('ONTOLOGY'):
                sub_df = sub_df.sort_values('p.adjust', ascending=False).tail(10)
                sub_df.insert(0, 'path', path)
                data_list.append(sub_df)
        else:
            data = data.sort_values('p.adjust', ascending=False).tail(10)
            data.insert(0, 'path', path)
            data_list.append(data)
        combined_data = pd.concat(data_list)
        combined_data.rename(columns={combined_data.columns[1]: 'Number'}, inplace=True)
        combined_data['Condition'] = ''
        combined_data['FC'] = ''
        combined_data['FDR'] = ''
        combined_data['Method'] = ''
        for i, row in combined_data.iterrows():
            combined_data['path_split'] = combined_data['path'].apply(lambda x: re.split(r'[/_]', x))
            # 将分割后的倒数第九和倒数第八个值传递给condition列
            combined_data['Condition'] = combined_data['path_split'].apply(lambda x: x[-9] + '_vs_' + x[-8])
            # 将分割后的倒数第七个值传递给FC列
            combined_data['FC'] = combined_data['path_split'].apply(lambda x: x[-7])
            # 将分割后的倒数第六个值传递给FDR列
            combined_data['FDR'] = combined_data['path_split'].apply(lambda x: x[-6])
            # 将分割后的倒数第三个值传递给Method列
            combined_data['Method'] = combined_data['path_split'].apply(lambda x: x[-3])

        combined_data = combined_data.drop(columns=['path_split'])
        # 获取 DataFrame 中的后四列
        columns_to_move = combined_data.columns[-4:]
        # 使用 reindex() 函数将后四列移动到前四列
        combined_data = combined_data.reindex(columns=columns_to_move.tolist() + combined_data.columns[:-4].tolist())
        combined_data.to_csv(f'{col}.csv', index=False, float_format='%.2e', header=True)
