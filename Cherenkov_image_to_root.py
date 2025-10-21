#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import ROOT

# ===== 設定 =====
in_csv  = "pi_20250610-4_data.csv"   # 入力CSV（カレントに置くかフルパス指定）
outroot = "out_hist2d.root"          # 出力ROOT
start_col_1based = 6                 # 6列目から使用（1始まり）
nx, ny = 8, 8                        # 配置（横8 × 縦8）

# ===== 前処理 =====
start_col = start_col_1based - 1     # Pythonは0始まり
need_vals = nx * ny                  # 64

# ===== CSV読み込み =====
rows = []
with open(in_csv, newline="") as f:
    rdr = csv.reader(f)
    for r in rdr:
        if not r:
            continue
        # 数値化（空欄対策をするなら try/except を追加）
        rows.append([float(x) for x in r])

# ===== ROOTファイル作成 =====
fout = ROOT.TFile(outroot, "RECREATE")

for i, r in enumerate(rows):
    if len(r) < start_col + need_vals:
        print(f"[warn] 行{i}は列数不足のためスキップ（{len(r)}列）")
        continue

    # 8×8 のブロックを取り出し
    block = r[start_col : start_col + need_vals]

    # TH2D作成（x:1..8, y:1..8）
    hname = f"h2_event_{i:06d}"
    h = ROOT.TH2D(hname, f"Event {i};X bin;Y bin", nx, 0.5, nx + 0.5, ny, 0.5, ny + 0.5)

    # 充填：6列目～13列目 → y=1（最下段）, x=1..8、次の8要素で y=2 … の順
    # 画面上で「上を y=8」にしたい場合は SetBinContent(ix+1, ny-iy, val) のように y を反転してください。
    for iy in range(ny):        # 0..7
        for ix in range(nx):    # 0..7
            val = block[iy * nx + ix]
            h.SetBinContent(ix + 1, iy + 1, val)

    h.Write()

fout.Close()
print(f"[done] wrote {outroot}")
