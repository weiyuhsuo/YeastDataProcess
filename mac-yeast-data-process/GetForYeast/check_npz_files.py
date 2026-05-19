#!/usr/bin/env python3
"""
检查npz文件是否完整和有效
"""
import numpy as np
import os
import sys

def check_npz_file(file_path):
    """检查单个npz文件"""
    print(f"\n检查文件: {file_path}")
    print("=" * 60)
    
    # 检查文件是否存在
    if not os.path.exists(file_path):
        print(f"❌ 文件不存在")
        return False
    
    # 检查文件大小
    file_size = os.path.getsize(file_path)
    print(f"文件大小: {file_size / (1024*1024):.2f} MB ({file_size:,} bytes)")
    
    if file_size == 0:
        print(f"❌ 文件为空")
        return False
    
    # 尝试加载文件
    try:
        npz_file = np.load(file_path, mmap_mode='r', allow_pickle=True)
        print(f"✅ 文件可以加载")
        
        # 列出所有键
        keys = list(npz_file.keys())
        print(f"包含的键: {keys}")
        
        # 检查关键数据
        if 'data' in keys:
            data = npz_file['data']
            print(f"  data形状: {data.shape}")
            print(f"  data类型: {data.dtype}")
            print(f"  data大小: {data.nbytes / (1024*1024):.2f} MB")
        else:
            print(f"⚠️ 警告: 未找到'data'键")
        
        if 'peak_ids' in keys:
            peak_ids = npz_file['peak_ids']
            print(f"  peak_ids数量: {len(peak_ids)}")
        else:
            print(f"⚠️ 警告: 未找到'peak_ids'键")
        
        if 'sample_ids' in keys:
            sample_ids = npz_file['sample_ids']
            print(f"  sample_ids数量: {len(sample_ids)}")
        else:
            print(f"⚠️ 警告: 未找到'sample_ids'键")
        
        npz_file.close()
        print(f"✅ 文件检查通过")
        return True
        
    except EOFError as e:
        print(f"❌ 文件损坏或格式错误 (EOFError)")
        print(f"   错误详情: {e}")
        print(f"   建议: 重新生成npz文件")
        return False
    except Exception as e:
        print(f"❌ 无法加载文件")
        print(f"   错误类型: {type(e).__name__}")
        print(f"   错误详情: {e}")
        return False

if __name__ == '__main__':
    # 默认检查KM的4个文件
    base_dir = os.path.dirname(os.path.abspath(__file__))
    files_to_check = [
        os.path.join(base_dir, 'input/KM/matrix_C1.csv.npz'),
        os.path.join(base_dir, 'input/KM/matrix_C3.csv.npz'),
        os.path.join(base_dir, 'input/KM/matrix_O2.csv.npz'),
        os.path.join(base_dir, 'input/KM/matrix_O3.csv.npz'),
    ]
    
    # 如果提供了命令行参数，使用命令行参数
    if len(sys.argv) > 1:
        files_to_check = sys.argv[1:]
    
    print("=" * 60)
    print("NPZ文件完整性检查")
    print("=" * 60)
    
    results = []
    for file_path in files_to_check:
        result = check_npz_file(file_path)
        results.append((file_path, result))
    
    print("\n" + "=" * 60)
    print("检查结果汇总")
    print("=" * 60)
    for file_path, result in results:
        status = "✅ 正常" if result else "❌ 异常"
        print(f"{status}: {os.path.basename(file_path)}")
    
    # 如果有异常文件，退出码为1
    if not all(result for _, result in results):
        sys.exit(1)
