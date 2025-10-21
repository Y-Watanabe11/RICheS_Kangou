#!/usr/bin/env bash
# csv2root_embed.sh
# CSVを選んで -> ch_eventtree作成 -> 全イベントの8x8ヒストをhist/に埋め込み
# さらに -> 各chのスペクトル(1D)をhist_spectrum/に作成（任意）
# さらに -> チェレンコフ・リング（3x3画像モーメント法）解析（任意）
# さらに -> DSSD ↔ チェレンコフ時刻同期（任意）

set -euo pipefail

# --- 準備 ---
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
MACRO_MAKE="${SCRIPT_DIR}/make_ch_tree_fixed.C"
MACRO_EMBED="${SCRIPT_DIR}/embed_all_hists.C"
MACRO_SPEC="${SCRIPT_DIR}/embed_channel_spectra.C"
MACRO_RING="${SCRIPT_DIR}/analyze_cherenkov_moments.C"   # ← モーメント法に変更
MACRO_SYNC="${SCRIPT_DIR}/sync_dssd_cherenkov.C"

# 依存チェック
command -v root >/dev/null 2>&1 || { echo "[ERR] ROOTが見つかりません。'root' コマンドを使えるようにしてください。"; exit 1; }
[[ -f "$MACRO_MAKE" ]] || { echo "[ERR] マクロが見つかりません: $MACRO_MAKE"; exit 1; }
[[ -f "$MACRO_EMBED" ]] || { echo "[ERR] マクロが見つかりません: $MACRO_EMBED"; exit 1; }
if [[ ! -f "$MACRO_SPEC" ]]; then
  echo "[WARN] スペクトル用マクロが見つかりません: $MACRO_SPEC"
  echo "       スペクトル作成はスキップします。"
  HAS_SPEC_MACRO=0
else
  HAS_SPEC_MACRO=1
fi

# --- CSV を選択（ドラッグ＆ドロップ想定） ---
read -r -p "CSVファイルをターミナルへドラッグ＆ドロップして Enter: " CSV_INPUT

# 前後空白と両端の引用符を除去
trim() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  [[ ${s:0:1} == "\"" || ${s:0:1} == "'" ]] && s="${s:1}"
  [[ ${s: -1} == "\"" || ${s: -1} == "'" ]] && s="${s:0:${#s}-1}"
  printf '%s' "$s"
}
CSV_PATH="$(trim "$CSV_INPUT")"

if [[ -z "$CSV_PATH" ]]; then
  echo "[ERR] CSVパスが空です。やり直してください。"; exit 1
fi
if [[ ! -f "$CSV_PATH" ]]; then
  echo "[ERR] ファイルが存在しません: $CSV_PATH"; exit 1
fi

# --- 出力ROOT名を自動決定（同ディレクトリに .root） ---
CSV_DIR="$(cd -- "$(dirname -- "$CSV_PATH")" && pwd)"
CSV_BASE="$(basename -- "$CSV_PATH")"
OUT_BASE="${CSV_BASE%.*}.root"
OUT_ROOT="${CSV_DIR}/${OUT_BASE}"

echo "[INFO] 入力CSV : $CSV_PATH"
echo "[INFO] 出力ROOT: $OUT_ROOT"

# --- 1) CSV -> ROOT（ch_eventtree 作成） ---
# make_ch_tree_fixed.C(csv, outroot, flipY)
CMD1=$(printf '%s("%s","%s",true)' "$MACRO_MAKE" "$CSV_PATH" "$OUT_ROOT")
echo "[RUN ] root -l -q \"$CMD1\""
root -l -q "$CMD1"

# --- 2) 全イベントのヒストを hist/ に埋め込み（上書き、flipY=true） ---
# embed_all_hists.C(rootfile, hist_dir, overwrite, flipY, first, last, every)
CMD2=$(printf '%s("%s","%s",true,true,0,-1,1)' "$MACRO_EMBED" "$OUT_ROOT" "hist")
echo "[RUN ] root -l -q \"$CMD2\""
root -l -q "$CMD2"

# --- 3) 各チャンネルのスペクトルを hist_spectrum/ に作成（任意） ---
if [[ $HAS_SPEC_MACRO -eq 1 ]]; then
  read -r -p "各chのスペクトル(hist_spectrum/)も作成しますか？ [Y/n]: " MAKE_SPEC
  MAKE_SPEC="${MAKE_SPEC:-Y}"
  if [[ "$MAKE_SPEC" =~ ^[Yy]$ ]]; then
    # 既定: nbins=512, 範囲=0..4095（1パス高速）
    read -r -p "ヒストのビン数（既定 512）: " NBINS
    NBINS="${NBINS:-512}"

    read -r -p "範囲を自動推定しますか？（2パス・遅い） [y/N]: " AUTORANGE
    AUTORANGE="${AUTORANGE:-N}"

    if [[ "$AUTORANGE" =~ ^[Yy]$ ]]; then
      # autoscan=true, xmin/xmaxはダミー
      CMD3=$(printf '%s("%s","%s",true,1,64,%d,0,1,true,-1)' \
            "$MACRO_SPEC" "$OUT_ROOT" "hist_spectrum" "$NBINS")
      echo "[RUN ] root -l -q \"$CMD3\""
      root -l -q "$CMD3"
    else
      read -r -p "x-min（既定 0）: " XMIN; XMIN="${XMIN:-0}"
      read -r -p "x-max（既定 4095）: " XMAX; XMAX="${XMAX:-4095}"
      # autoscan=false, 指定範囲で1パス高速
      CMD3=$(printf '%s("%s","%s",true,1,64,%d,%s,%s,false,-1)' \
            "$MACRO_SPEC" "$OUT_ROOT" "hist_spectrum" "$NBINS" "$XMIN" "$XMAX")
      echo "[RUN ] root -l -q \"$CMD3\""
      root -l -q "$CMD3"
    fi
  else
    echo "[SKIP] スペクトル作成はスキップしました。"
  fi
fi

echo "[INFO] ここまで: ch_eventtree / hist / (hist_spectrum) を作成済み"

# --- 4) チェレンコフ・リング解析（3x3画像モーメント法、総光子しきい値で暗電流を除外） ---
if [[ -f "$MACRO_RING" ]]; then
  echo
  read -r -p "チェレンコフ・リング解析を実行しますか？ [Y/n]: " DO_RING
  DO_RING="${DO_RING:-Y}"
  if [[ "$DO_RING" =~ ^[Yy]$ ]]; then
    read -r -p "総光子しきい値（sumQ >= ? 既定 50）: " MINPH
    MINPH="${MINPH:-50}"
    read -r -p "リング診断Canvasを保存しますか？ [y/N]: " SAVE_C
    SAVE_C="${SAVE_C:-N}"
    SAVE_FLAG=false; [[ "$SAVE_C" =~ ^[Yy]$ ]] && SAVE_FLAG=true

    PIXPITCH="6.0"   # mm（必要に応じて変更）
    # analyze_cherenkov_moments(root, min_photons, flipY, save_canv, canv_dir, pix_pitch_mm)
    CMD4=$(printf '%s("%s",%s,true,%s,"%s",%s)' \
          "$MACRO_RING" "$OUT_ROOT" "$MINPH" "$SAVE_FLAG" "ring_canv" "$PIXPITCH")
    echo "[RUN ] root -l -q \"$CMD4\""
    root -l -q "$CMD4"
  fi
else
  echo "[WARN] analyze_cherenkov_moments.C が見つからないためリング解析をスキップしました。"
fi

# --- 5) DSSD ↔ チェレンコフ 時刻同期（相互最近傍 + 閾値、未対応は別ツリーへ） ---
if [[ -f "$MACRO_SYNC" ]]; then
  echo
  read -r -p "DSSD とチェレンコフの時刻同期を実行しますか？ [Y/n]: " DO_SYNC
  DO_SYNC="${DO_SYNC:-Y}"
  if [[ "$DO_SYNC" =~ ^[Yy]$ ]]; then
    # DSSD ROOT を選択
    read -r -p "DSSDのROOTファイルをドラッグ＆ドロップして Enter: " DSSD_INPUT
    DSSD_ROOT="$(trim "$DSSD_INPUT")"
    if [[ -z "$DSSD_ROOT" || ! -f "$DSSD_ROOT" ]]; then
      echo "[ERR] DSSD ROOTファイルが無効です: $DSSD_ROOT"; exit 1
    fi

    # ツリー/枝名の既定
    read -r -p "CH側ツリー名（既定 ch_eventtree）: " CH_TREE;  CH_TREE="${CH_TREE:-ch_eventtree}"
    read -r -p "CH側時刻ブランチ名（既定 time）: " CH_TIME;    CH_TIME="${CH_TIME:-time}"
    read -r -p "DSSD側ツリー名（既定 eventtree）: " DSSD_TREE; DSSD_TREE="${DSSD_TREE:-eventtree}"
    read -r -p "DSSD側時刻ブランチ名（既定 ti）: " DSSD_TIME; DSSD_TIME="${DSSD_TIME:-ti}"

    # 既定パラメータ（あなたの仕様）
    CH_WRAP_SEC="356.0"          # CH: 1usカウンタが356秒で周回
    DSSD_TICK_SEC="1e-8"         # DSSD: 10 ns / tick
    DSSD_MODULUS="4294967296.0"  # 2^32
    read -r -p "許容時間差 [秒]（既定 0.001 = 1 ms）: " MAX_DT_SEC
    MAX_DT_SEC="${MAX_DT_SEC:-0.001}"
    read -r -p "一対一対応にしますか？（同じDSSDを再使用しない） [Y/n]: " OTO
    OTO="${OTO:-Y}"; ONE_TO_ONE=true; [[ "$OTO" =~ ^[Nn]$ ]] && ONE_TO_ONE=false

    # 出力ファイルの選択：空なら CH ファイル(OUT_ROOT)へ追記
    read -r -p "同期結果の出力ROOT（空=チェレンコフ側に追記）: " OUT_SYNC_INPUT
    OUT_SYNC_INPUT="$(trim "$OUT_SYNC_INPUT")"
    if [[ -z "$OUT_SYNC_INPUT" ]]; then
      OUT_SYNC_ARG='""'
    else
      OUT_SYNC_ARG=$(printf '"%s"' "$OUT_SYNC_INPUT")
    fi

    # sync_dssd_cherenkov.C(
    #   cher_root, cher_tree, cher_time, ch_wrap_sec,
    #   dssd_root, dssd_tree, dssd_time, dssd_tick_sec, dssd_modulus,
    #   max_dt_sec, one_to_one, out_root)
    CMD5=$(printf '%s("%s","%s","%s",%s,"%s","%s","%s",%s,%s,%s,%s,%s)' \
          "$MACRO_SYNC" \
          "$OUT_ROOT" "$CH_TREE" "$CH_TIME" "$CH_WRAP_SEC" \
          "$DSSD_ROOT" "$DSSD_TREE" "$DSSD_TIME" "$DSSD_TICK_SEC" "$DSSD_MODULUS" \
          "$MAX_DT_SEC" "$ONE_TO_ONE" "$OUT_SYNC_ARG")

    echo "[RUN ] root -l -q \"$CMD5\""
    root -l -q "$CMD5"
  fi
else
  echo "[WARN] sync_dssd_cherenkov.C が見つからないため時刻同期をスキップしました。"
fi

echo "[DONE] すべて完了しました。TBrowserで確認できます:"
echo "       root -l \"$OUT_ROOT\""
echo "       root [0] TBrowser b;"
echo "  → hist/ → histtree_******"
echo "  → hist_spectrum/ → spec_ch**（作成した場合）"
echo "  → ring_canv/（保存した場合）"
echo "  → sync/ → sync_tree, unmatched_ch, unmatched_dssd（同期した場合）"
