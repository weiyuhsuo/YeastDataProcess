#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
快速构建质粒npy文件的脚本
支持命令行参数和配置文件选择
"""

import sys
import argparse
from build_plasmid_npy_general import build_plasmid_npy, PlasmidConfig
from plasmid_configs import PLASMID_CONFIGS, get_config_by_name, list_available_configs

def build_from_config_name(name, promoter_name="", copy_number=None):
    """根据配置名称构建npy文件"""
    config = get_config_by_name(name, promoter_name, copy_number)
    if not config:
        print(f"错误: 未找到配置 {name}")
        if promoter_name:
            print(f"启动子: {promoter_name}")
        if copy_number:
            print(f"拷贝数: {copy_number}")
        return False
    
    print(f"使用配置: {config.name}")
    return build_plasmid_npy(config)

def build_from_custom_config(fimo_file, expression_file, sample_info_file, 
                           target_gene, copy_number, name="Custom", promoter_name=""):
    """使用自定义配置构建npy文件"""
    config = PlasmidConfig(
        name=name,
        fimo_file=fimo_file,
        expression_file=expression_file,
        sample_info_file=sample_info_file,
        target_gene=target_gene,
        copy_number=copy_number,
        promoter_name=promoter_name
    )
    
    print(f"使用自定义配置: {config.name}")
    return build_plasmid_npy(config)

def main():
    parser = argparse.ArgumentParser(description='快速构建质粒npy文件')
    parser.add_argument('--list', action='store_true', help='列出所有可用配置')
    parser.add_argument('--name', type=str, help='质粒名称')
    parser.add_argument('--promoter', type=str, help='启动子名称')
    parser.add_argument('--copy', type=int, help='拷贝数')
    parser.add_argument('--fimo', type=str, help='Fimo文件路径')
    parser.add_argument('--expression', type=str, help='表达数据文件路径')
    parser.add_argument('--sample-info', type=str, help='样本信息文件路径')
    parser.add_argument('--target-gene', type=str, help='目标基因')
    parser.add_argument('--interactive', action='store_true', help='交互式选择配置')
    
    args = parser.parse_args()
    
    if args.list:
        list_available_configs()
        return
    
    if args.interactive:
        # 交互式选择
        list_available_configs()
        try:
            choice = int(input("\n请选择要构建的质粒配置 (输入序号): ")) - 1
            if 0 <= choice < len(PLASMID_CONFIGS):
                selected_config = PLASMID_CONFIGS[choice]
                print(f"\n已选择: {selected_config.name}")
                success = build_plasmid_npy(selected_config)
                if success:
                    print(f"\n{selected_config.name}质粒npy文件构建成功！")
                else:
                    print(f"\n{selected_config.name}质粒npy文件构建失败！")
            else:
                print("无效选择！")
        except ValueError:
            print("请输入有效的数字！")
        except KeyboardInterrupt:
            print("\n用户取消操作")
        return
    
    # 检查是否有足够的参数进行自定义构建
    if all([args.fimo, args.expression, args.sample_info, args.target_gene, args.copy]):
        success = build_from_custom_config(
            fimo_file=args.fimo,
            expression_file=args.expression,
            sample_info_file=args.sample_info,
            target_gene=args.target_gene,
            copy_number=args.copy,
            name=args.name or "Custom",
            promoter_name=args.promoter or ""
        )
        if success:
            print("自定义配置构建成功！")
        else:
            print("自定义配置构建失败！")
        return
    
    # 使用预定义配置
    if args.name:
        success = build_from_config_name(args.name, args.promoter, args.copy)
        if success:
            print(f"{args.name}质粒npy文件构建成功！")
        else:
            print(f"{args.name}质粒npy文件构建失败！")
        return
    
    # 如果没有参数，显示帮助信息
    print("请使用以下方式之一:")
    print("1. --list: 列出所有可用配置")
    print("2. --interactive: 交互式选择配置")
    print("3. --name <质粒名> [--promoter <启动子>] [--copy <拷贝数>]: 使用预定义配置")
    print("4. 提供完整的自定义参数: --fimo, --expression, --sample-info, --target-gene, --copy")
    print("\n示例:")
    print("  python build_quick.py --list")
    print("  python build_quick.py --interactive")
    print("  python build_quick.py --name STE12 --promoter promoter1 --copy 10")
    print("  python build_quick.py --fimo data/New/fimo.tsv --expression data/New/expr.csv --sample-info data/New/info.csv --target-gene YNEW001 --copy 20")

if __name__ == "__main__":
    main()
