cd ./CompleteGenome_noplasmid  # 存储.fna文件的目录
mkdir prokka_res  # prokka计算结果目录
for fna_file in *.fna; do
    res_file_name="${fna_file%.fna}"
    prokka --outdir "./prokka_res/$res_file_name" --prefix "$res_file_name" "$fna_file"
done