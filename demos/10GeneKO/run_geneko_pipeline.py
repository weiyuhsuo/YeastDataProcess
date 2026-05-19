#!/usr/bin/env python3
"""
GeneKO数据处理完整流程运行脚本

执行顺序：
1. 数据预处理 (process_geneko_data.py)
2. 构建numpy文件 (build_geneko_numpy.py)

使用方法：
python run_geneko_pipeline.py
"""

import os
import sys
import subprocess
import time

def run_command(command, description):
    """运行命令并显示进度"""
    print(f"\n{'='*60}")
    print(f"开始执行: {description}")
    print(f"{'='*60}")
    print(f"命令: {command}")
    
    start_time = time.time()
    
    try:
        # 运行命令
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        
        # 显示输出
        if result.stdout:
            print("标准输出:")
            print(result.stdout)
        
        if result.stderr:
            print("标准错误:")
            print(result.stderr)
        
        end_time = time.time()
        print(f"\n✅ {description} 执行成功！耗时: {end_time - start_time:.2f} 秒")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ {description} 执行失败！")
        print(f"错误代码: {e.returncode}")
        if e.stdout:
            print("标准输出:")
            print(e.stdout)
        if e.stderr:
            print("标准错误:")
            print(e.stderr)
        return False

def check_requirements():
    """检查运行环境"""
    print("检查运行环境...")
    
    # 检查Python包 - 使用更可靠的导入方式
    required_packages = [
        ('pandas', 'pandas'),
        ('numpy', 'numpy'), 
        ('sklearn', 'scikit-learn'),
        ('tqdm', 'tqdm'),
        ('matplotlib', 'matplotlib'),
        ('seaborn', 'seaborn')
    ]
    
    missing_packages = []
    for import_name, package_name in required_packages:
        try:
            if import_name == 'sklearn':
                import sklearn
            else:
                __import__(import_name)
            print(f"  ✅ {package_name}")
        except ImportError:
            print(f"  ❌ {package_name}")
            missing_packages.append(package_name)
    
    if missing_packages:
        print(f"\n❌ 缺少必要的Python包: {missing_packages}")
        print("请使用以下命令安装:")
        print(f"pip install {' '.join(missing_packages)}")
        return False
    
    print("✅ Python包检查通过")
    
    # 检查数据文件
    required_files = [
        'data/Saccharomyces_cerevisiae.gene_info',
        'data/ncbiRefSeqCurated.txt',
        'data/ATAC1_matrix.csv',
        'data/ATAC1_peaks.narrowPeak'
    ]
    
    missing_files = []
    for file_path in required_files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
            print(f"  ❌ {file_path}")
        else:
            print(f"  ✅ {file_path}")
    
    if missing_files:
        print(f"\n❌ 缺少必要的数据文件:")
        for file_path in missing_files:
            print(f"  {file_path}")
        return False
    
    print("✅ 数据文件检查通过")
    
    # 检查GeneKO数据目录
    geneko_data_dir = 'data/20250801data'
    if not os.path.exists(geneko_data_dir):
        print(f"❌ GeneKO数据目录不存在: {geneko_data_dir}")
        return False
    
    gse_folders = [d for d in os.listdir(geneko_data_dir) if d.startswith('GSE')]
    if not gse_folders:
        print(f"❌ 在 {geneko_data_dir} 中未找到GSE文件夹")
        return False
    
    print(f"✅ 找到 {len(gse_folders)} 个GSE数据集: {gse_folders}")
    
    return True

def main():
    """主函数"""
    print("GeneKO数据处理完整流程")
    print("="*60)
    
    # 检查运行环境
    if not check_requirements():
        print("\n❌ 环境检查失败，请解决上述问题后重试")
        sys.exit(1)
    
    # 记录开始时间
    total_start_time = time.time()
    
    # 步骤1: 数据预处理
    if not run_command("python process_geneko_data.py", "GeneKO数据预处理"):
        print("\n❌ 数据预处理失败，流程终止")
        sys.exit(1)
    
    # 步骤2: 构建numpy文件
    if not run_command("python build_geneko_numpy.py", "构建GeneKO numpy文件"):
        print("\n❌ numpy文件构建失败，流程终止")
        sys.exit(1)
    
    # 计算总耗时
    total_end_time = time.time()
    total_time = total_end_time - total_start_time
    
    # 输出总结
    print(f"\n{'='*60}")
    print("🎉 GeneKO数据处理完整流程执行成功！")
    print(f"{'='*60}")
    print(f"总耗时: {total_time:.2f} 秒 ({total_time/60:.2f} 分钟)")
    print(f"输出目录: NumpyFileOutput/")
    
    # 检查输出文件
    output_dir = "NumpyFileOutput"
    if os.path.exists(output_dir):
        output_files = [f for f in os.listdir(output_dir) if f.endswith('.npy')]
        print(f"生成的numpy文件数量: {len(output_files)}")
        for file_name in sorted(output_files):
            file_path = os.path.join(output_dir, file_name)
            file_size = os.path.getsize(file_path) / (1024 * 1024)  # MB
            print(f"  {file_name} ({file_size:.2f} MB)")
    
    print(f"\n下一步：使用生成的numpy文件进行基因KO预测模型训练和测试")

if __name__ == "__main__":
    main()
