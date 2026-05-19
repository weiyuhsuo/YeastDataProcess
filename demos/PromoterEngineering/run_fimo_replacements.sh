#!/usr/bin/env bash
set -euo pipefail

# ================= 配置（集中管理路径与参数） =================
MOTIF_MEME="$(dirname "$0")/yeast_jaspar_motifs.meme"
INPUT_FASTA="$(dirname "$0")/replacements/replaced_sequences.fasta"
OUT_DIR="$(dirname "$0")/fimo_out_replacements"

# ================= 运行 =================
mkdir -p "$OUT_DIR"

echo "运行 FIMO..."
# 不设置任何默认阈值或verbosity，将用户传入参数原样透传
# 用法示例（与你平时一致）：
#   fimo [你的参数...] motif.meme replaced_sequences.fasta
# 现在等效为：
#   bash run_fimo_replacements.sh [你的参数...]
fimo \
  --oc "$OUT_DIR" \
  "$@" \
  "$MOTIF_MEME" \
  "$INPUT_FASTA"

echo "完成。FIMO 输出在：$OUT_DIR/ (包含 fimo.tsv)"


