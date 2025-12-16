import trimesh
import pandas as pd


mesh = trimesh.load("plane.obj", force='mesh')
verts = mesh.vertices       
faces = mesh.faces           
print(mesh)

df = pd.read_csv("point.csv")
print(df.head())

# 约束点
src_points = df[['x', 'y', 'z']].values
tar_points = df[['x1', 'y1', 'z1']].values
print("Source Points:\n", src_points)