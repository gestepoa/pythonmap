#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
甘肃省地级市 GeoJSON 转单笔闭合欧拉回路 SVG 工具（精简版）
====================================================================
核心流程：
1. 经纬度余弦投影与画布自适应缩放。
2. 构建无向平面图并自动合并相邻地级市公共边界。
3. 拓扑保持抽稀（可调节精度）。
4. 最小生成树（MST）自动连接平凉内部飞地/孤岛。
5. 中国邮递员算法（最小权完美匹配）消除奇度顶点并提取单一闭合欧拉回路。
6. 导出仅含单个 <path> 标签的标准 SVG 文件。

依赖库安装：
    pip install networkx numpy
"""

import json
import math
import os
from typing import Dict, List, Tuple
import networkx as nx
import numpy as np

# 坐标节点类型定义：(量化后的经度, 量化后的纬度)
NodeKey = Tuple[float, float]
Point2D = Tuple[float, float]


# ==============================================================================
# 第一部分：几何预处理、投影与拓扑保持抽稀
# ==============================================================================


def load_geojson(filepath: str) -> dict:
  """读取 GeoJSON 文件。"""
  if not os.path.exists(filepath):
    raise FileNotFoundError(
        f"未找到文件 '{filepath}'。请将您提供的甘肃省 GeoJSON 内容保存为"
        f" '{filepath}' 并放在同一目录下。"
    )
  with open(filepath, "r", encoding="utf-8") as f:
    return json.load(f)


def extract_rings_from_geojson(geojson_data: dict) -> List[List[NodeKey]]:
  """从 GeoJSON 的 FeatureCollection 中提取所有 Polygon/MultiPolygon 的闭合环（含外环与内环孔洞）。

  对经纬度保留 6 位小数，确保相邻地级市（如兰州与白银、张掖与山丹马场）的公共边界顶点完全一致。
  """
  rings: List[List[NodeKey]] = []
  features = geojson_data.get("features", [])

  for feat in features:
    geom = feat.get("geometry")
    if not geom:
      continue
    gtype = geom.get("type")
    coords = geom.get("coordinates", [])

    if gtype == "Polygon":
      polys = [coords]
    elif gtype == "MultiPolygon":
      polys = coords
    else:
      continue

    for poly in polys:
      for ring in poly:
        # 量化到 6 位小数以消除浮点误差并对齐公共边界顶点
        clean_ring: List[NodeKey] = []
        for pt in ring:
          node = (round(float(pt[0]), 6), round(float(pt[1]), 6))
          if not clean_ring or clean_ring[-1] != node:
            clean_ring.append(node)
        if len(clean_ring) >= 3:
          # 确保首尾闭合
          if clean_ring[0] != clean_ring[-1]:
            clean_ring.append(clean_ring[0])
          rings.append(clean_ring)

  return rings


def build_planar_graph(
    rings: List[List[NodeKey]],
    canvas_width: float = 1000.0,
    canvas_height: float = 1000.0,
    padding: float = 50.0,
) -> nx.Graph:
  """将所有多边形环构建为无向平面图 G，并完成经纬度投影与公共边自动去重。

  返回带有节点属性 'pos': (x, y) 和边属性 'weight': 欧氏距离 的无向图。
  """
  all_lons = [pt[0] for ring in rings for pt in ring]
  all_lats = [pt[1] for ring in rings for pt in ring]

  min_lon, max_lon = min(all_lons), max(all_lons)
  min_lat, max_lat = min(all_lats), max(all_lats)

  # 甘肃省跨纬度较大（约 32.5°N ~ 42.8°N），取中心纬度做余弦投影校正，防止地图横向拉伸
  mean_lat_rad = math.radians((min_lat + max_lat) / 2.0)
  cos_lat = math.cos(mean_lat_rad)

  geo_width = (max_lon - min_lon) * cos_lat
  geo_height = max_lat - min_lat

  avail_w = canvas_width - 2 * padding
  avail_h = canvas_height - 2 * padding
  scale = min(avail_w / geo_width, avail_h / geo_height)

  # 居中偏移量计算
  offset_x = padding + (avail_w - geo_width * scale) / 2.0
  offset_y = padding + (avail_h - geo_height * scale) / 2.0

  def project_point(lon: float, lat: float) -> Point2D:
    # SVG 坐标系中 X 向右为正，Y 向下为正（因此纬度需要用 max_lat - lat 翻转）
    x = offset_x + (lon - min_lon) * cos_lat * scale
    y = offset_y + (max_lat - lat) * scale
    return (x, y)

  G = nx.Graph()

  for ring in rings:
    for i in range(len(ring) - 1):
      u, v = ring[i], ring[i + 1]
      if u == v:
        continue
      if u not in G:
        G.add_node(u, pos=project_point(u[0], u[1]))
      if v not in G:
        G.add_node(v, pos=project_point(v[0], v[1]))

      # 无向图 G.add_edge 会自动把两个相邻市共用的边界合并为 1 条边
      if not G.has_edge(u, v):
        pos_u = G.nodes[u]["pos"]
        pos_v = G.nodes[v]["pos"]
        dist = math.hypot(pos_u[0] - pos_v[0], pos_u[1] - pos_v[1])
        G.add_edge(u, v, weight=dist, edge_type="normal")

  return G


def rdp_simplify_chain(
    chain: List[NodeKey], pos_dict: Dict[NodeKey, Point2D], epsilon: float
) -> List[NodeKey]:
  """Ramer-Douglas-Peucker (RDP) 折线抽稀算法，保持链条首尾交界点严格不变。"""
  if len(chain) <= 2 or epsilon <= 0:
    return chain

  p_start = np.array(pos_dict[chain[0]])
  p_end = np.array(pos_dict[chain[-1]])
  line_vec = p_end - p_start
  line_len = np.linalg.norm(line_vec)

  max_dist = -1.0
  index = -1

  for i in range(1, len(chain) - 1):
    p = np.array(pos_dict[chain[i]])
    if line_len == 0:
      dist = float(np.linalg.norm(p - p_start))
    else:
      # 计算点到线段的垂直距离
      t = np.clip(np.dot(p - p_start, line_vec) / (line_len**2), 0.0, 1.0)
      proj = p_start + t * line_vec
      dist = float(np.linalg.norm(p - proj))

    if dist > max_dist:
      max_dist = dist
      index = i

  if max_dist > epsilon:
    left = rdp_simplify_chain(chain[: index + 1], pos_dict, epsilon)
    right = rdp_simplify_chain(chain[index:], pos_dict, epsilon)
    return left[:-1] + right
  else:
    return [chain[0], chain[-1]]


def simplify_planar_graph_topology(G: nx.Graph, epsilon: float = 0.8) -> nx.Graph:
  """在不破坏地级市公共边界拓扑的前提下，对图中度为 2 的边界链条进行抽稀。

  避免直接对各市 Polygon 单独抽稀导致的边界裂缝问题。
  """
  if epsilon <= 0:
    return G

  pos_dict = nx.get_node_attributes(G, "pos")
  junctions = {n for n, d in G.degree() if d != 2}

  visited_edges = set()
  chains: List[List[NodeKey]] = []

  # 1. 提取以交界点（degree != 2）为端点的所有简单链
  for u in junctions:
    for v in G.neighbors(u):
      edge_key = tuple(sorted([u, v]))
      if edge_key in visited_edges:
        continue
      chain = [u, v]
      visited_edges.add(edge_key)
      curr = v
      prev = u
      while curr not in junctions:
        neighbors = list(G.neighbors(curr))
        nxt = neighbors[0] if neighbors[1] == prev else neighbors[1]
        edge_key = tuple(sorted([curr, nxt]))
        if edge_key in visited_edges:
          break
        visited_edges.add(edge_key)
        chain.append(nxt)
        prev, curr = curr, nxt
      chains.append(chain)

  # 2. 提取完全不经过交界点的纯闭合环（例如平凉市内部的孤立飞地孔洞）
  for u, v in G.edges():
    edge_key = tuple(sorted([u, v]))
    if edge_key in visited_edges:
      continue
    chain = [u, v]
    visited_edges.add(edge_key)
    curr = v
    prev = u
    while True:
      neighbors = list(G.neighbors(curr))
      nxt = neighbors[0] if neighbors[1] == prev else neighbors[1]
      edge_key = tuple(sorted([curr, nxt]))
      if edge_key in visited_edges:
        chain.append(nxt)
        break
      visited_edges.add(edge_key)
      chain.append(nxt)
      prev, curr = curr, nxt
    chains.append(chain)

  # 3. 对每条链执行 RDP 抽稀并重建平面图
  G_simp = nx.Graph()
  for chain in chains:
    if chain[0] == chain[-1] and len(chain) > 4:
      # 闭合环拆成两半分别抽稀，防止首尾退化为单条线段
      mid = len(chain) // 2
      c1 = rdp_simplify_chain(chain[: mid + 1], pos_dict, epsilon)
      c2 = rdp_simplify_chain(chain[mid:], pos_dict, epsilon)
      simp_chain = c1[:-1] + c2
    else:
      simp_chain = rdp_simplify_chain(chain, pos_dict, epsilon)

    for i in range(len(simp_chain) - 1):
      a, b = simp_chain[i], simp_chain[i + 1]
      if a == b:
        continue
      if a not in G_simp:
        G_simp.add_node(a, pos=pos_dict[a])
      if b not in G_simp:
        G_simp.add_node(b, pos=pos_dict[b])
      pa, pb = pos_dict[a], pos_dict[b]
      G_simp.add_edge(
          a, b, weight=math.hypot(pa[0] - pb[0], pa[1] - pb[1]), edge_type="normal"
      )

  return G_simp


# ==============================================================================
# 第二部分：连通孤岛（搭桥）与中国邮递员算法（构造欧拉回路）
# ==============================================================================


def connect_isolated_components(G: nx.Graph) -> Tuple[nx.Graph, List[Tuple[NodeKey, NodeKey]]]:
  """检测图中的孤立连通分量（如平凉市内部的飞地孔洞或独立保护区），

  使用最小生成树（MST）在最近顶点之间自动架设“桥梁线”（Bridge Lines），使全图连通。
  """
  comps = [list(c) for c in nx.connected_components(G)]
  num_comps = len(comps)
  bridges: List[Tuple[NodeKey, NodeKey]] = []

  if num_comps <= 1:
    return G, bridges

  print(f"[拓扑检测] 发现 {num_comps} 个独立连通分量（含内部飞地/孔洞），正在计算最短搭桥方案...")
  pos_dict = nx.get_node_attributes(G, "pos")

  # 预提取每个连通分量的坐标矩阵，加速最近点对搜索
  comp_coords = []
  for comp in comps:
    pts = np.array([pos_dict[n] for n in comp])
    comp_coords.append(pts)

  comp_graph = nx.Graph()
  for i in range(num_comps):
    for j in range(i + 1, num_comps):
      pts_i = comp_coords[i]
      pts_j = comp_coords[j]
      # 计算两个分量所有顶点对之间的欧氏距离矩阵
      diff = pts_i[:, np.newaxis, :] - pts_j[np.newaxis, :, :]
      dist_matrix = np.hypot(diff[:, :, 0], diff[:, :, 1])
      min_idx = np.unravel_index(np.argmin(dist_matrix), dist_matrix.shape)
      min_dist = float(dist_matrix[min_idx])
      u = comps[i][min_idx[0]]
      v = comps[j][min_idx[1]]
      comp_graph.add_edge(i, j, weight=min_dist, bridge=(u, v))

  # 求连通分量间的最小生成树，保证搭桥总长度最短
  mst = nx.minimum_spanning_tree(comp_graph, weight="weight")
  for _, _, data in mst.edges(data=True):
    u, v = data["bridge"]
    w = data["weight"]
    G.add_edge(u, v, weight=w, edge_type="bridge")
    bridges.append((u, v))
    print(f"  -> 已架设桥梁连接孤岛: {u} <---> {v} (画布距离: {w:.2f} px)")

  return G, bridges


def build_eulerian_circuit(G: nx.Graph) -> Tuple[List[Point2D], nx.MultiGraph]:
  """解决加权中国邮递员问题（Chinese Postman Problem）：

  1. 找出所有奇度顶点（度为奇数的交叉路口）。
  2. 计算奇度顶点间的最短路径距离，构造完全图并求最小权完美匹配。
  3. 沿着匹配路径复制边（使画笔在这些公共州界/市界上原路折返或重走一次），令全图所有顶点度均为偶数。
  4. 使用 Hierholzer 算法提取单一闭合欧拉回路。
  """
  if not nx.is_connected(G):
    G, _ = connect_isolated_components(G)

  M = nx.MultiGraph()
  for n, data in G.nodes(data=True):
    M.add_node(n, **data)
  for u, v, data in G.edges(data=True):
    M.add_edge(u, v, **data)

  odd_nodes = [n for n, d in M.degree() if d % 2 == 1]
  print(f"[欧拉化] 检测到 {len(odd_nodes)} 个奇度顶点（三市交界T型路口），正在求解最小权完美匹配...")

  if odd_nodes:
    K_odd = nx.Graph()
    shortest_paths: Dict[Tuple[NodeKey, NodeKey], List[NodeKey]] = {}

    # 仅需从奇度顶点出发运行单源 Dijkstra（甘肃全图奇度点通常仅 ~30 个，耗时极短）
    for i, u in enumerate(odd_nodes):
      lengths, paths = nx.single_source_dijkstra(G, u, weight="weight")
      for j in range(i + 1, len(odd_nodes)):
        v = odd_nodes[j]
        K_odd.add_edge(u, v, weight=lengths[v])
        shortest_paths[(u, v)] = paths[v]
        shortest_paths[(v, u)] = list(reversed(paths[v]))

    # 求解奇度顶点的最小权完美匹配
    try:
      matching = nx.min_weight_matching(K_odd, weight="weight")
    except AttributeError:
      # 兼容旧版本 NetworkX
      max_w = max(d["weight"] for _, _, d in K_odd.edges(data=True)) + 1.0
      for _, _, d in K_odd.edges(data=True):
        d["neg_weight"] = max_w - d["weight"]
      matching = nx.max_weight_matching(
          K_odd, maxcardinality=True, weight="neg_weight"
      )

    # 将匹配出的最短路径添加为重复边（对应视频中黄点重走部分边界的过程）
    dup_edge_count = 0
    for u, v in matching:
      path = shortest_paths[(u, v)]
      for k in range(len(path) - 1):
        a, b = path[k], path[k + 1]
        w = G[a][b]["weight"]
        M.add_edge(a, b, weight=w, edge_type="duplicate")
        dup_edge_count += 1

    print(f"[欧拉化] 匹配完成！共添加了 {dup_edge_count} 条重复走线边，全图所有顶点度数现均为偶数。")

  assert nx.is_eulerian(M), "构建欧拉图失败：图中仍存在奇度顶点或非连通分量！"

  # 选择最北端（酒泉方向）的顶点作为一笔画起点
  pos_dict = nx.get_node_attributes(M, "pos")
  start_node = min(M.nodes(), key=lambda n: pos_dict[n][1])

  circuit_edges = list(nx.eulerian_circuit(M, source=start_node))
  euler_coords: List[Point2D] = [pos_dict[circuit_edges[0][0]]]
  for _, v in circuit_edges:
    euler_coords.append(pos_dict[v])

  print(f"[欧拉回路] 成功提取单一闭合回路！总计包含 {len(euler_coords)} 个连续轨迹点。")
  return euler_coords, M


# ==============================================================================
# 第三部分：导出单路径闭合 SVG 文件
# ==============================================================================


def export_single_path_svg(
    coords: List[Point2D],
    output_svg_path: str,
    canvas_width: float = 1000.0,
    canvas_height: float = 1000.0,
) -> None:
  """将欧拉回路坐标序列导出为仅含单一 <path> 闭合标签的标准 SVG 文件。"""
  if not coords:
    raise ValueError("坐标序列为空，无法生成 SVG。")

  path_cmds = [f"M {coords[0][0]:.3f},{coords[0][1]:.3f}"]
  for x, y in coords[1:]:
    path_cmds.append(f"L {x:.3f},{y:.3f}")
  path_cmds.append("Z")
  d_attr = " ".join(path_cmds)

  svg_content = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     viewBox="0 0 {int(canvas_width)} {int(canvas_height)}"
     width="{int(canvas_width)}" height="{int(canvas_height)}">
  <path id="gansu-euler-circuit"
        d="{d_attr}"
        fill="none"
        stroke="#000000"
        stroke-width="1.5"
        stroke-linecap="round"
        stroke-linejoin="round"/>
</svg>
"""
  with open(output_svg_path, "w", encoding="utf-8") as f:
    f.write(svg_content)
  print(f"[SVG 导出] 单一闭合路径 SVG 已成功生成: '{output_svg_path}'")


# ==============================================================================
# 主程序入口
# ==============================================================================

if __name__ == "__main__":
  # 自动匹配您上传的文件名 '甘肃省_市.geojson' 或 'gansu.geojson'
  INPUT_GEOJSON = (
      "./chinamap/陕西省_市.geojson"
      if os.path.exists("./chinamap/陕西省_市.geojson")
      else "./chinamap/陕西省_省.geojson"
  )
  OUTPUT_SVG = "euler_single_path.svg"

  if not os.path.exists(INPUT_GEOJSON):
    print("未找到输入文件！请将您的 GeoJSON 文件命名为 '陕西省_市.geojson' 并放在脚本同目录下。")
  else:
    print(f"====== 正在读取 '{INPUT_GEOJSON}' 并构建单笔闭合 SVG ======")
    raw_geojson = load_geojson(INPUT_GEOJSON)

    # 1. 提取多边形环并构建平面图（自动合并相邻市界）
    rings = extract_rings_from_geojson(raw_geojson)
    G_raw = build_planar_graph(rings, canvas_width=1000, canvas_height=1000, padding=50)
    print(f"[平面图] 合并公共边后共有 {G_raw.number_of_nodes()} 个顶点, {G_raw.number_of_edges()} 条边。")

    # 2. 拓扑保持抽稀（epsilon=0.5 可滤除冗余密集点；若想保留 100% 原始精度可设为 0）
    G_simp = simplify_planar_graph_topology(G_raw, epsilon=0.5)
    print(f"[拓扑抽稀] 抽稀后保留 {G_simp.number_of_nodes()} 个顶点, {G_simp.number_of_edges()} 条边。")

    # 3. 自动搭桥连接孤岛并求解中国邮递员问题，生成欧拉回路坐标序列
    euler_coords, _ = build_eulerian_circuit(G_simp)

    # 4. 导出最终的纯路径 SVG 文件
    export_single_path_svg(euler_coords, OUTPUT_SVG)