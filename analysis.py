from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os

def analysis_all_files_CDS(file_path):
    for file in os.listdir(file_path):
        if file.endswith(".fasta"):
            for record in SeqIO.parse(os.path.join(file_path, file), "fasta"):
                longest_CDS = ""
                longest_length = 0
                
                # 遍历序列，寻找所有潜在的起始密码子 ATG
                start_pos = record.seq.find("ATG")
                while start_pos != -1:
                    # 从当前 ATG 翻译，遇到终止密码子停止
                    temp_protein = record.seq[start_pos:].translate(to_stop=True)
                    current_cds_length = len(temp_protein) * 3
                    
                    # 如果当前阅读框比之前记录的最长阅读框更长，则更新结果
                    if current_cds_length > longest_length:
                        longest_length = current_cds_length
                        longest_CDS = record.seq[start_pos : start_pos + current_cds_length]
                        
                    # 从当前 ATG 的下一个碱基开始，继续寻找下一个 ATG
                    start_pos = record.seq.find("ATG", start_pos + 1)
                
                # 确保找到了至少一个有效的 CDS
                if longest_length > 0:
                    gc_content = f"{(gc_fraction(longest_CDS)*100):.2f}%"
                    yield file, longest_length, gc_content

def calculate_fasta_stats(target_dir="."):
    """
    遍历指定文件夹下的所有 FASTA 文件，计算文件名、整条序列的总长度和总 GC 含量。
    """
    for file in os.listdir(target_dir):
        # 兼容常见的 fasta 后缀名
        if file.endswith((".fasta", ".fa", ".fas")):
            file_path = os.path.join(target_dir, file)
            
            # 解析 FASTA 文件（支持单个文件内包含多条序列的情况）
            for record in SeqIO.parse(file_path, "fasta"):
                # 获取整条序列（mRNA）的总长度
                total_length = len(record.seq)
                
                # 计算整条序列的总 GC 含量，保留两位小数并转为百分比字符串
                gc_content = f"{(gc_fraction(record.seq) * 100):.2f}%"
                
                # 返回：文件名, 序列长度, GC含量
                yield file, total_length, gc_content
            
def in_csv(file, CDS_length, gc_content):
    csv_path = r"D:\biocode\3cytokine\cytokine_summary.csv"
    # 仅首次运行时写入表头（用a+模式判断文件是否为空）
    with open(csv_path, "a+", encoding='utf-8') as f:
        f.seek(0)#将文件指针移动到开头，以便检查文件是否为空
        if not f.read():
            f.write("File_Name,CDS_Length,GC_Content\n")
        # 修复编码参数错误（原把encoding写进内容里）+ 移除多余逗号
        f.write(f"{file},{CDS_length},{gc_content}\n")

if __name__ == "__main__":
    file_path = "D:\\biocode\\3cytokine"
    for result in analysis_all_files_CDS(file_path):
        in_csv(result[0], result[1], result[2])
