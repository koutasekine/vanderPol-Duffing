# 遅延 van der Pol–Duffing 方程式の分数調波解に対する精度保証付き数値計算

論文「遅延van der Pol-Duffing方程式の分数調波解に対する無限次元ガウスの消去法を用いた精度保証付き数値計算法」の数値実験に用いたプログラムと実行結果です．

## 内容

| ファイル | 内容 |
|---|---|
| `improve_verifyVDP.cpp`, `improve_VanDerPol.hpp` | 提案手法（Table 1, 3） |
| `verifyVDP_bdab.cpp`, `VanDerPol.hpp` | 漸近優対角理論（Table 2, 4） |
| `run_vdp_tables.sh` | コンパイルと実行を行い，ログと要約表を保存するスクリプト |
| `results/` | 論文の表の元になった実行結果 |

## 必要なもの

- kv ライブラリ: http://verifiedby.me/kv/
- VCP ライブラリ: https://github.com/koutasekine/vcp
- Intel MKL，MPFR

kv と VCP は，このディレクトリの 1 つ上に置きます（`-I..` でインクルード）．

## コンパイルと実行

```sh
bash run_vdp_tables.sh
```

結果は `./results/<日付時刻>_<hostname>/` に保存されます．

個別に実行する場合:

```sh
./verifyVDP_improve.out [m] [alpha] [beta] [tomega] [tau] [gamma] [mu] [k] [n]
./verifyVDP_bdab.out    [m] [alpha] [beta] [tomega] [tau] [gamma] [mu] [k] [n]
```

`tomega` は (1.2) の $\tilde\omega$ です（(1.6) の $\omega=\tilde\omega/n$）．論文の例:

```sh
./verifyVDP_improve.out 70 1 5 1 0.1 1 0.1 0.1 2
```

## 出力

要約表（`summary_*.tsv`）の値は区間の上端で，論文の表はこれを有効数字 4 桁に上向きに丸めたものです．
