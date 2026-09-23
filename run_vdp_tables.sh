#!/usr/bin/env bash
# =====================================================================
# run_vdp_tables.sh
#   improve_verifyVDP.cpp（提案手法）と verifyVDP_bdab.cpp（漸近対角優越法）を
#   コンパイルし，指定した m で順に実行して，ログと要約表を保存する．
#
#   置き場所: improve_verifyVDP.cpp と verifyVDP_bdab.cpp があるフォルダ
#             （-I.. でインクルードするので，vcp/ と kv/ はその 1 つ上にある前提）
#   実行方法: bash run_vdp_tables.sh
#   出力先  : ./results/<日付>_<hostname>/
#             ├─ build.log                    コンパイルのログ
#             ├─ env.txt                      実行環境・ソースの md5・パラメータ
#             ├─ improve_m<m>_<hostname>.log  提案手法の各 m のログ
#             ├─ bdab_m<m>_<hostname>.log     既存手法の各 m のログ
#             ├─ summary_improve_<hostname>.tsv  提案手法の要約（区間の上端）
#             └─ summary_bdab_<hostname>.tsv     既存手法の要約（区間の上端）
# =====================================================================
set -u

# ---------------- 設定 ----------------
M_LIST="55 60 70 80 90 100 110 200"

# 引数の順: [m_app] [alpha] [beta] [tomega] [tau] [gamma] [mu] [k] [n]
ALPHA=1
BETA=5
TOMEGA=1      # (1.2) の外力の角周波数 tilde{omega}（(1.6) の omega = tomega/n）
TAU=0.1
GAMMA=1
MU=0.1
K=0.1
N=2

RUN_IMPROVE=1   # 0 にすると提案手法を実行しない
RUN_BDAB=1      # 0 にすると既存手法を実行しない

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64}"
LIBS="${LIBS:--L${MKLROOT:-}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp}"

BIN_IMPROVE=./verifyVDP_improve.out
BIN_BDAB=./verifyVDP_bdab.out
# --------------------------------------

HOST=$(hostname)
STAMP=$(date +%Y%m%d-%H%M%S)
OUTDIR="./results/${STAMP}_${HOST}"
mkdir -p "$OUTDIR"

echo "Output directory: $OUTDIR"

# ---- 実行環境の記録 ----
{
  echo "host      : $HOST"
  echo "date      : $(date '+%Y-%m-%d %H:%M:%S %z')"
  echo "compiler  : $($CXX --version | head -1)"
  echo "CXXFLAGS  : $CXXFLAGS"
  echo "LIBS      : $LIBS"
  echo "OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}  MKL_NUM_THREADS=${MKL_NUM_THREADS:-unset}"
  echo "params    : alpha=$ALPHA beta=$BETA tomega=$TOMEGA tau=$TAU gamma=$GAMMA mu=$MU k=$K n=$N"
  echo "m list    : $M_LIST"
  echo "--- md5 of sources ---"
  md5sum improve_verifyVDP.cpp improve_VanDerPol.hpp verifyVDP_bdab.cpp VanDerPol.hpp 2>&1
} > "$OUTDIR/env.txt"

# ---- コンパイル ----
: > "$OUTDIR/build.log"
if [ "$RUN_IMPROVE" = 1 ]; then
  echo "Compiling improve_verifyVDP.cpp -> $BIN_IMPROVE"
  if ! $CXX $CXXFLAGS improve_verifyVDP.cpp -o "$BIN_IMPROVE" $LIBS >> "$OUTDIR/build.log" 2>&1; then
    echo "ERROR: compile failed (improve). See $OUTDIR/build.log"; exit 1
  fi
fi
if [ "$RUN_BDAB" = 1 ]; then
  echo "Compiling verifyVDP_bdab.cpp -> $BIN_BDAB"
  if ! $CXX $CXXFLAGS verifyVDP_bdab.cpp -o "$BIN_BDAB" $LIBS >> "$OUTDIR/build.log" 2>&1; then
    echo "ERROR: compile failed (bdab). See $OUTDIR/build.log"; exit 1
  fi
fi

# ---- 区間 [a,b] の上端を取り出す（見つからなければ "-"）----
# 使い方: upper "<行頭の文字列(正規表現)>" <logfile>
upper() {
  local v
  v=$(grep -E "$1" "$2" | tail -1 | sed -n 's/.*\[[^,]*,\([^]]*\)\].*/\1/p')
  [ -n "$v" ] && echo "$v" || echo "-"
}
result_of() {
  if grep -q "Verification Success" "$1"; then echo "success"
  elif grep -q "kk >= 1" "$1"; then echo "fail(kk>=1)"
  elif grep -q "Verification Failure" "$1"; then echo "fail"
  else echo "abnormal_end"; fi
}

# ---- 提案手法 ----
if [ "$RUN_IMPROVE" = 1 ]; then
  S="$OUTDIR/summary_improve_${HOST}.tsv"
  printf "m\tK0\tMm\tKmN\tsigmaK0\tkappa\tOldkappa\tCFinv\tFinvB\tFzh_Z\tetaZ\tetaD\tFinv_BZ\tFinv_BZD\tb_r0\terror_D\tresult\tminimal_period\tsec\n" > "$S"
  for m in $M_LIST; do
    LOG="$OUTDIR/improve_m${m}_${HOST}.log"
    echo "[improve] m=$m ..."
    t0=$(date +%s)
    "$BIN_IMPROVE" "$m" "$ALPHA" "$BETA" "$TOMEGA" "$TAU" "$GAMMA" "$MU" "$K" "$N" > "$LOG" 2>&1
    t1=$(date +%s)
    mp="-"
    if grep -q "Minimal period of z\* is 2" "$LOG"; then mp="2pi"
    elif grep -q "Minimal period is NOT verified" "$LOG"; then mp="not_verified"; fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$m" \
      "$(upper '^K0 =' "$LOG")" \
      "$(upper '^Mm =' "$LOG")" \
      "$(upper '^KmN =' "$LOG")" \
      "$(upper '^sigmam\*K0 =' "$LOG")" \
      "$(upper '^kappa =' "$LOG")" \
      "$(upper '^\(Old\) kappa =' "$LOG")" \
      "$(upper '^\|\| C Finv \|\| =' "$LOG")" \
      "$(upper '^\|\| Finv B \|\| =' "$LOG")" \
      "$(upper '^\|\| F\(zh\) \|\|_Z <=' "$LOG")" \
      "$(upper "^\|\| F'\[zh\]\^-1 F\(zh\) \|\|_Z <=" "$LOG")" \
      "$(upper "^\|\| F'\[zh\]\^-1 F\(zh\) \|\|_D <=" "$LOG")" \
      "$(upper '^\|\| DFinv \|\|_\{B\(Z\)\} <=' "$LOG")" \
      "$(upper '^\|\| DFinv \|\|_\{B\(Z, D\)\} <=' "$LOG")" \
      "$(upper '^b\(r0\) =' "$LOG")" \
      "$(upper '^\|\| u\* - u\^ \|\| <=' "$LOG")" \
      "$(result_of "$LOG")" \
      "$mp" \
      "$((t1-t0))" >> "$S"
    echo "          -> $(result_of "$LOG")  ($((t1-t0)) s)"
  done
fi

# ---- 既存手法（漸近対角優越）----
if [ "$RUN_BDAB" = 1 ]; then
  S="$OUTDIR/summary_bdab_${HOST}.tsv"
  printf "m\tK0\tkk_inf\tkk_one\tMn\tM\tFzh_Z\teta\tb_r0\tresult\tsec\n" > "$S"
  for m in $M_LIST; do
    LOG="$OUTDIR/bdab_m${m}_${HOST}.log"
    echo "[bdab]    m=$m ..."
    t0=$(date +%s)
    "$BIN_BDAB" "$m" "$ALPHA" "$BETA" "$TOMEGA" "$TAU" "$GAMMA" "$MU" "$K" "$N" > "$LOG" 2>&1
    t1=$(date +%s)
    # (kk) の行は infinity ノルム版 → 1 ノルム版の順に出力される
    kk1=$(grep -E '^\(kk\):' "$LOG" | sed -n '1s/.*\[[^,]*,\([^]]*\)\].*/\1/p'); [ -z "$kk1" ] && kk1="-"
    kk2=$(grep -E '^\(kk\):' "$LOG" | sed -n '2s/.*\[[^,]*,\([^]]*\)\].*/\1/p'); [ -z "$kk2" ] && kk2="-"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$m" \
      "$(upper '^K0 =' "$LOG")" \
      "$kk1" "$kk2" \
      "$(upper '^Mn =' "$LOG")" \
      "$(upper '^M =' "$LOG")" \
      "$(upper '^\|\| F\(zh\) \|\|_Z <=' "$LOG")" \
      "$(upper '^eta =' "$LOG")" \
      "$(upper '^b\(r0\) =' "$LOG")" \
      "$(result_of "$LOG")" \
      "$((t1-t0))" >> "$S"
    echo "          -> $(result_of "$LOG")  ($((t1-t0)) s)"
  done
fi

echo
echo "Done. Results are in: $OUTDIR"
ls -1 "$OUTDIR"
