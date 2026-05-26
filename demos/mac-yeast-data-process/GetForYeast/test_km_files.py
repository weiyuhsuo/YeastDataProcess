#!/usr/bin/env python3
"""测试KM npz文件是否能正常加载（模拟推理脚本的加载过程）"""
import numpy as np
import os

# 使用与推理脚本相同的路径
base_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(base_dir, 'input/KM')

files_to_check = [
    'matrix_C1.csv.npz',
    'matrix_C3.csv.npz',
    'matrix_O2.csv.npz',
    'matrix_O3.csv.npz',
]

print("=" * 70)
print("测试 KM NPZ 文件（模拟推理脚本加载）")
print("=" * 70)
print(f"数据目录: {data_dir}")
print()

all_success = True

for filename in files_to_check:
    file_path = os.path.join(data_dir, filename)
    print(f"测试文件: {filename}")
    print("-" * 70)
    
    if not os.path.exists(file_path):
        print(f"❌ 文件不存在: {file_path}")
        all_success = False
        continue
    
    file_size = os.path.getsize(file_path)
    print(f"文件大小: {file_size / (1024*1024):.2f} MB")
    
    # 使用与推理脚本相同的加载方式（先标准模式，失败再mmap）
    try:
        # 方法1: 标准模式（不使用mmap）
        npz_file = np.load(file_path, allow_pickle=True)
        print(f"✅ 标准模式加载成功")
        
        # 检查关键数据
        if 'data' not in npz_file:
            print(f"❌ 缺少'data'键")
            all_success = False
            continue
        
        data = npz_file['data']
        num_samples, num_peaks, num_features = data.shape
        print(f"  数据形状: ({num_samples}, {num_peaks}, {num_features})")
        
        # 检查特征维度（应该是545特征 + 2标签 = 547）
        if num_features == 547:
            print(f"  ✅ 特征维度正确: 545特征 + 2标签 = 547")
        elif num_features == 545:
            print(f"  ✅ 特征维度: 545（仅特征，无标签，适合推理）")
        else:
            print(f"  ⚠️ 特征维度异常: {num_features}（期望547或545）")
        
        # 检查peak_ids
        if 'peak_ids' in npz_file:
            peak_ids = npz_file['peak_ids']
            if hasattr(peak_ids, 'tolist'):
                peak_ids = peak_ids.tolist()
            print(f"  Peak IDs数量: {len(peak_ids)}")
        else:
            print(f"  ⚠️ 缺少'peak_ids'键")
        
        # 检查sample_ids
        if 'sample_ids' in npz_file:
            sample_ids = npz_file['sample_ids']
            if hasattr(sample_ids, 'tolist'):
                sample_ids = sample_ids.tolist()
            print(f"  Sample IDs数量: {len(sample_ids)}")
        else:
            print(f"  ⚠️ 缺少'sample_ids'键")
        
        # 尝试访问数据（模拟推理过程）
        try:
            sample_features = data[0, 0, :545]  # 第一个样本，第一个peak，前545个特征
            print(f"  ✅ 可以访问数据: 特征向量长度={len(sample_features)}")
        except Exception as e:
            print(f"  ❌ 无法访问数据: {e}")
            all_success = False
        
        npz_file.close()
        print(f"✅ 文件测试通过\n")
        
    except EOFError as e:
        print(f"❌ EOFError: {e}")
        print(f"   尝试mmap模式...")
        try:
            npz_file = np.load(file_path, mmap_mode='r', allow_pickle=True)
            print(f"✅ mmap模式加载成功")
            data = npz_file['data']
            print(f"  数据形状: {data.shape}")
            npz_file.close()
            print(f"✅ 文件测试通过（使用mmap模式）\n")
        except Exception as e2:
            print(f"❌ mmap模式也失败: {e2}")
            all_success = False
            print()
    except Exception as e:
        print(f"❌ 加载失败: {type(e).__name__}: {e}")
        all_success = False
        print()

print("=" * 70)
if all_success:
    print("✅ 所有文件测试通过，可以用于推理")
else:
    print("❌ 部分文件测试失败，请检查")
print("=" * 70)
