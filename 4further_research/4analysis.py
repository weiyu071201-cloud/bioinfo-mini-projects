from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os
import pandas as pd

def find_longest_cds(sequence):
    """寻找并提取序列中的最长 CDS"""
    longest_cds_seq = ""
    max_len = 0
    start_pos = sequence.find("ATG")
    while start_pos != -1:
        # translate(to_stop=True) 模拟翻译
        temp_protein = sequence[start_pos:].translate(to_stop=True)
        current_len = len(temp_protein) * 3
        if current_len > max_len:
            max_len = current_len
            longest_cds_seq = sequence[start_pos : start_pos + current_len]
        start_pos = sequence.find("ATG", start_pos + 1)
    return longest_cds_seq, max_len

def get_species_name(description):
    """
    鲁棒地提取物种名：支持 PREDICTED 标签及多单词学名
    示例输入: "XM_0123.1 PREDICTED: Canis lupus familiaris interleukin 6..."
    """
    # 1. 统一转为大写处理后再移除 PREDICTED: 标签
    desc_clean = description.replace("PREDICTED:", "").replace("predicted:", "").strip()
    parts = desc_clean.split()
    
    if len(parts) < 2:
        return "Unknown"

    # parts[0] 是编号（如 NM_000600.5），物种名从 parts[1] 开始
    species_parts = []
    # 定义基因关键词列表，用于标识物种名结束的位置
    stop_words = ["interleukin", "insulin", "preproinsulin"]
    
    for word in parts[1:]:
        # 如果遇到基因相关的词汇，停止提取
        if word.lower().strip(",()") in stop_words:
            break
        species_parts.append(word)
    
    # 重新组合物种名，如 "Canis lupus familiaris" 或 "Homo sapiens"
    return " ".join(species_parts).strip(",")

def get_gene_name(description):
    """从FASTA表头识别基因名称"""
    desc_lower = description.lower()
    if "interleukin" in desc_lower or "il-6" in desc_lower or "il6" in desc_lower:
        return "IL6"
    elif "insulin" in desc_lower or "ins" in desc_lower:
        return "INS"
    else:
        return "Unknown"


def extract_and_sort_data(file_path):
    """遍历文件夹，提取 CDS 并在内部完成物种排序"""
    data_list = []
    for file in os.listdir(file_path):
        if file.endswith((".fasta", ".fa", ".fas")):
            full_path = os.path.join(file_path, file)
            for record in SeqIO.parse(full_path, "fasta"):
                cds_seq, cds_len = find_longest_cds(record.seq)
                if cds_len > 0:
                    # 提取精确的物种名
                    species = get_species_name(record.description)
                    gene = get_gene_name(record.description) 
                    
                    data_list.append({
                        "Gene": gene,
                        "Species": species,
                        "CDS_Length": cds_len,
                        "GC_Content":  f"{(gc_fraction(record.seq) * 100):.2f}%",
                    })
    
    df = pd.DataFrame(data_list)
    # 按物种名排序，确保 IL-6 和 INS 组的物种顺序完全一致
    df_sorted = df.sort_values(by=["Species", "Gene"], ascending=True).reset_index(drop=True)
    return df_sorted

if __name__ == "__main__":
    # 使用你之前定义的路径
    input_dir = r"D:\biocode\3cytokine"
    output_csv = r"D:\biocode\4further_research\IL6_INS_summary.csv"
    
    result_df = extract_and_sort_data(input_dir)
    result_df.to_csv(output_csv, index=False, encoding='utf-8',mode='a')
    
    print("数据提取并按物种排序完成：")
    print(result_df[["Gene","Species", "GC_Content"]].head())