#!/usr/bin/env python3
"""检查KM npz文件"""
import numpy as np
import os
import sys

base_dir = 'KM/numpy/output/run_20260129_152130'
files_to_check = [
    'matrix_C1.csv.npz',
    'matrix_C3.csv.npz',
    'matrix_O2.csv.npz',
    'matrix_O3.csv.npz',
]

print("=" * 70)
print("检查 NPZ 文件完整性")
print("=" * 70)

results = []
for filename in files_to_check:
    file_path = os.path.join(base_dir, filename)
    print(f"\n检查文件: {filename}")
    print("-" * 70)
    
    if not os.path.exists(file_path):
        print(f"❌ 文件不存在")
        results.append((filename, False, "文件不存在"))
        continue
    
    file_size = os.path.getsize(file_path)
    print(f"文件大小: {file_size / (1024*1024):.2f} MB ({file_size:,} bytes)")
    
    if file_size == 0:
        print(f"❌ 文件为空")
        results.append((filename, False, "文件为空"))
        continue
    
    # 尝试不同的加载方式
    success = False
    error_msg = None
    
    # 方法1: 使用mmap
    try:
        npz_file = np.load(file_path, mmap_mode='r', allow_pickle=True)
        print(f"✅ 方法1 (mmap) 成功")
        success = True
    except EOFError as e:
        print(f"❌ 方法1 (mmap) 失败: EOFError - {e}")
        error_msg = f"EOFError: {e}"
        
        # 方法2: 不使用mmap
        try:
            npz_file = np.load(file_path, allow_pickle=True)
            print(f"✅ 方法2 (无mmap) 成功")
            success = True
        except Exception as e2:
            print(f"❌ 方法2 (无mmap) 也失败: {type(e2).__name__} - {e2}")
            error_msg = f"{type(e2).__name__}: {e2}"
            results.append((filename, False, error_msg))
            continue
    except Exception as e:
        print(f"❌ 方法1 (mmap) 失败: {type(e).__name__} - {e}")
        error_msg = f"{type(e).__name__}: {e}"
        
        # 方法2: 不使用mmap
        try:
            npz_file = np.load(file_path, allow_pickle=True)
            print(f"✅ 方法2 (无mmap) 成功")
            success = True
        except Exception as e2:
            print(f"❌ 方法2 (无mmap) 也失败: {type(e2).__name__} - {e2}")
            error_msg = f"{type(e2).__name__}: {e2}"
            results.append((filename, False, error_msg))
            continue
    
    if success:
        try:
            keys = list(npz_file.keys())
            print(f"包含的键: {keys}")
            
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
            results.append((filename, True, None))
        except Exception as e:
            print(f"❌ 读取数据时出错: {type(e).__name__} - {e}")
            results.append((filename, False, f"读取错误: {e}"))

print("\n" + "=" * 70)
print("检查结果汇总")
print("=" * 70)
all_ok = True
for filename, result, error in results:
    status = "✅ 正常" if result else "❌ 异常"
    print(f"{status}: {filename}")
    if error:
        print(f"   错误: {error}")
        all_ok = False

if not all_ok:
    print("\n⚠️ 发现异常文件，建议重新生成npz文件")
    sys.exit(1)
else:
    print("\n✅ 所有文件检查通过")
