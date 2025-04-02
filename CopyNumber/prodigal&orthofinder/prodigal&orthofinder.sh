# 先prodigal
cd ./CompleteGenome_noplasmid  # 存储.fna文件的目录
mkdir prodigal_res  # prodigal计算结果目录
for fna_file in *.fna; do
    res_file_name="${fna_file%.fna}"
    prodigal -i "$fna_file" -a "./prodigal_res/$res_file_name.faa"
done
# 再orthofinder
orthofinder -f ./prodigal_res -t 8 -a 8