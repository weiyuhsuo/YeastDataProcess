#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
质粒配置文件
包含所有质粒的配置信息，方便快速添加新的质粒
"""

from build_plasmid_npy_general import PlasmidConfig

# 质粒配置列表
PLASMID_CONFIGS = [
    # Cup1质粒配置
    PlasmidConfig(
        name="Cup1",
        fimo_file="data/FimoofCup1/fimo_Cup1promoter.tsv",
        expression_file="data/OE数据/OE_测试数据_Rip.csv",
        sample_info_file="data/OE数据/OE_样本信息_Rip.csv",
        target_gene="YEL024W",
        copy_number=100
    ),
    
    # STE12质粒配置 - promoter1, copy10
    PlasmidConfig(
        name="STE12",
        fimo_file="data/STE12/promoter1/fimo.tsv",
        expression_file="data/STE12/STE12表达矩阵.csv",
        sample_info_file="data/STE12/STE12样品信息.csv",
        target_gene="YHR084W",
        copy_number=10,
        promoter_name="promoter1"
    ),
    
    # STE12质粒配置 - promoter2, copy5
    PlasmidConfig(
        name="STE12",
        fimo_file="data/STE12/promoter2/fimo.tsv",
        expression_file="data/STE12/STE12表达矩阵.csv",
        sample_info_file="data/STE12/STE12样品信息.csv",
        target_gene="YHR084W",
        copy_number=5,
        promoter_name="promoter2"
    ),
    
    # 可以在这里添加更多质粒配置
    # 例如：
    # PlasmidConfig(
    #     name="NewPlasmid",
    #     fimo_file="data/NewPlasmid/fimo.tsv",
    #     expression_file="data/NewPlasmid/表达矩阵.csv",
    #     sample_info_file="data/NewPlasmid/样品信息.csv",
    #     target_gene="YNEW001",
    #     copy_number=20,
    #     promoter_name="promoter1"
    # ),
]

def get_config_by_name(name, promoter_name="", copy_number=None):
    """根据名称、启动子和拷贝数获取配置"""
    for config in PLASMID_CONFIGS:
        if (config.name == name and 
            (not promoter_name or config.promoter_name == promoter_name) and
            (copy_number is None or config.copy_number == copy_number)):
            return config
    return None

def list_available_configs():
    """列出所有可用的配置"""
    print("可用的质粒配置:")
    for i, config in enumerate(PLASMID_CONFIGS):
        print(f"{i+1}. {config.name} - {config.target_gene} - 拷贝数{config.copy_number}")
        if config.promoter_name:
            print(f"   启动子: {config.promoter_name}")
        print(f"   Fimo文件: {config.fimo_file}")
        print(f"   表达数据: {config.expression_file}")
        print(f"   样本信息: {config.sample_info_file}")
        print()

def add_new_config(name, fimo_file, expression_file, sample_info_file, 
                  target_gene, copy_number, promoter_name=""):
    """添加新的质粒配置"""
    new_config = PlasmidConfig(
        name=name,
        fimo_file=fimo_file,
        expression_file=expression_file,
        sample_info_file=sample_info_file,
        target_gene=target_gene,
        copy_number=copy_number,
        promoter_name=promoter_name
    )
    
    # 检查是否已存在相同配置
    existing = get_config_by_name(name, promoter_name, copy_number)
    if existing:
        print(f"警告: 已存在相同配置 {name}")
        return False
    
    PLASMID_CONFIGS.append(new_config)
    print(f"成功添加新配置: {name}")
    return True

if __name__ == "__main__":
    list_available_configs()
