#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
修复 p.csv 中的顶点索引
从 plane.obj 中比对原坐标，把原坐标替换成索引（0-based）
格式：原坐标x,原坐标y,原坐标z,变换后x,变换后y,变换后z
输出：索引,变换后x,变换后y,变换后z
"""

def read_obj_vertices(obj_file):
    """读取 OBJ 文件中的所有顶点，返回顶点列表（0-based 索引）"""
    vertices = []
    with open(obj_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('v '):
                parts = line.split()
                if len(parts) >= 4:
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    vertices.append((x, y, z))
    return vertices

def read_csv_points(csv_file):
    """读取 CSV 文件中的控制点
    格式：原坐标x,原坐标y,原坐标z,变换后x,变换后y,变换后z
    """
    points = []
    with open(csv_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            if len(parts) >= 6:
                orig_x = float(parts[0])
                orig_y = float(parts[1])
                orig_z = float(parts[2])
                new_x = float(parts[3])
                new_y = float(parts[4])
                new_z = float(parts[5])
                points.append((orig_x, orig_y, orig_z, new_x, new_y, new_z))
    return points

def find_matching_vertex(vertices, target_x, target_y, target_z, tolerance=1e-5):
    """在顶点列表中查找匹配的顶点，返回 0-based 索引"""
    for i, (vx, vy, vz) in enumerate(vertices):
        if abs(vx - target_x) < tolerance and abs(vy - target_y) < tolerance and abs(vz - target_z) < tolerance:
            return i
    return None

def fix_csv_indices(obj_file, csv_file, output_file):
    """修复 CSV 文件中的索引"""
    print(f"读取 OBJ 文件: {obj_file}")
    vertices = read_obj_vertices(obj_file)
    print(f"找到 {len(vertices)} 个顶点")
    
    print(f"读取 CSV 文件: {csv_file}")
    points = read_csv_points(csv_file)
    print(f"找到 {len(points)} 个控制点")
    
    fixed_points = []
    not_found = []
    
    for orig_x, orig_y, orig_z, new_x, new_y, new_z in points:
        idx = find_matching_vertex(vertices, orig_x, orig_y, orig_z)
        if idx is not None:
            fixed_points.append((idx, new_x, new_y, new_z))
            print(f"找到匹配: 原坐标 ({orig_x}, {orig_y}, {orig_z}) -> 索引 {idx}, 变换后 ({new_x}, {new_y}, {new_z})")
        else:
            not_found.append((orig_x, orig_y, orig_z, new_x, new_y, new_z))
            print(f"警告: 未找到匹配的顶点 (原坐标: {orig_x}, {orig_y}, {orig_z})")
    
    if not_found:
        print(f"\n警告: {len(not_found)} 个点未找到匹配的顶点")
        print("尝试使用更大的容差...")
        for orig_x, orig_y, orig_z, new_x, new_y, new_z in not_found:
            idx = find_matching_vertex(vertices, orig_x, orig_y, orig_z, tolerance=1e-3)
            if idx is not None:
                fixed_points.append((idx, new_x, new_y, new_z))
                print(f"找到匹配 (容差 1e-3): 原坐标 ({orig_x}, {orig_y}, {orig_z}) -> 索引 {idx}")
            else:
                print(f"错误: 仍然未找到匹配的顶点 (原坐标: {orig_x}, {orig_y}, {orig_z})")
    
    # 写入修复后的 CSV 文件
    # 格式：索引,变换后x,变换后y,变换后z
    print(f"\n写入修复后的文件: {output_file}")
    with open(output_file, 'w') as f:
        for idx, new_x, new_y, new_z in fixed_points:
            f.write(f"{idx},{new_x},{new_y},{new_z}\n")
    
    print(f"完成! 成功处理 {len(fixed_points)} 个点的索引")

if __name__ == '__main__':
    import sys
    obj_file = 'plane.obj'
    csv_file = 'p.csv'
    output_file = 'p.csv'
    
    if len(sys.argv) > 1:
        obj_file = sys.argv[1]
    if len(sys.argv) > 2:
        csv_file = sys.argv[2]
    if len(sys.argv) > 3:
        output_file = sys.argv[3]
    
    fix_csv_indices(obj_file, csv_file, output_file)

