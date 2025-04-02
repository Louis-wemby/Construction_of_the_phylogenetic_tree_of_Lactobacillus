cd ./CompleteGenome_noplasmid  # 存储.fna文件的目录
cd prokka_res
# ffn
target_dir="ffn_file"
mkdir -p "$target_dir"
find . -type f -name "*.ffn" -print0 | while IFS= read -r -d '' file; do
    filename=$(basename "$file")
    cp -v --backup=numbered "$file" "$target_dir/$filename"
done
# tsv
target_dir="tsv_file"
mkdir -p "$target_dir"
find . -type f -name "*.tsv" -print0 | while IFS= read -r -d '' file; do
    filename=$(basename "$file")
    cp -v --backup=numbered "$file" "$target_dir/$filename"
done