import pandas as pd

# 读取 CSV 文件
file_path = 'prefecture_forecast.csv'  # 替换为你的文件路径

# 读取数据
df = pd.read_csv(file_path)
print(df)

# 定义新的顺序
new_order = [23, 5, 2, 12, 38, 18, 40, 7, 21, 10, 34, 1, 28, 8, 17, 3, 37, 46, 14, 39, 43, 26, 24, 4, 45, 20, 42, 29, 15, 44, 33, 27, 41, 11, 25, 32, 22, 9, 36, 13, 31, 16, 30, 6, 35, 19]

# 创建一个新的 DataFrame 按照新顺序重新排列
df_reordered = df.iloc[[i-1 for i in new_order]]  # 减1是因为Python的索引从0开始

# 将重新排序的数据写回 CSV
output_file_path = 'prefecture_forecast_new.csv'  # 替换为你想要保存的文件路径
df_reordered.to_csv(output_file_path, index=False)

print('数据已按照新顺序重新排序并保存到', output_file_path)