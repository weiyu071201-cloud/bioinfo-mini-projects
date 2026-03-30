from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os

def analysis_all_files(file_path):
    for file in os.listdir(file_path):
        if file.endswith(".fasta"):
            for record in SeqIO.parse(os.path.join(file_path, file), "fasta"):
                start_pos = record.seq.find("ATG")#给人看的时候加1就行了，在代码里就从0开始算
                if start_pos != -1:
                    real_protein = record.seq[start_pos:].translate(to_stop=True)
                    stop_pos = start_pos + (len(real_protein) * 3)#编码区不需要加3
                    CDS_seq = record.seq[start_pos:stop_pos]
                CDS_length = len(CDS_seq)
                gc_content = f"{(gc_fraction(record.seq)*100):.2f}%"#gc含量两位百分数实质上是字符串
                yield file, CDS_length, gc_content#yield能返回多个结果
            
def in_csv(file, CDS_length, gc_content):
    csv_path = r"D:\biocode\2project_insulin_analysis\insulin_summary.csv"
    # 仅首次运行时写入表头（用a+模式判断文件是否为空）
    with open(csv_path, "a+", encoding='utf-8') as f:
        f.seek(0)#将文件指针移动到开头，以便检查文件是否为空
        if not f.read():
            f.write("File_Name,CDS_Length,GC_Content\n")
        # 修复编码参数错误（原把encoding写进内容里）+ 移除多余逗号
        f.write(f"{file},{CDS_length},{gc_content}\n")

if __name__ == "__main__":
    file_path = "D:\\biocode\\2project_insulin_analysis"
    for result in analysis_all_files(file_path):
        in_csv(result[0], result[1], result[2])

# #原来的代码中，in_csv函数的参数传递有误，导致写入CSV文件时出现了错误。现在已经修复了这个问题，并且添加了注释以便更好地理解代码的功能。
# with open(csv_path, "w") as f:
#         f.write("File_Name,CDS_Length,GC_Content\n")   
#     with open(csv_path, "a") as f:
#         f.write(f"{file},{CDS_length},{gc_content} , encoding='utf-8'\n")
