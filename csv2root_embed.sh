#!/usr/bin/env bash
# csv2root_embed.sh
# CSVを選んで -> ch_eventtree作成 -> 全イベントの8x8ヒストをhist/に埋め込み
# さらに -> 各chのスペクトル(1D)をhist_spectrum/に作成（任意）

set -euo pipefail

# --- 準備 ---
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
MACRO_MAKE="${SCRIPT_DIR}/make_ch_tree_fixed.C"
MACRO_EMBED="${SCRIPT_DIR}/embed_all_hists.C"
MACRO_SPEC="${SCRIPT_DIR}/embed_channel_spectra.C"   # ← 追加

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
  # 前後空白
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  # 先頭/末尾のシングル/ダブルクォートを剥がす
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

echo "[DONE] すべて完了しました。TBrowserで確認できます:"
echo "       root -l \"$OUT_ROOT\""
echo "       root [0] TBrowser b;"
echo "  → 左ペインの hist/ -> histtree_000000 などをダブルクリック"
echo "  → 左ペインの hist_spectrum/ -> spec_ch01, spec_ch02 ... をダブルクリック（作成した場合）"
