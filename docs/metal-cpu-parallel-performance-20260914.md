# Metal: CPU 並列化と CPU 最適化前からの累積高速化

2026-09-14、Apple M4 Max、macOS 26.6.2 で測定。
**Hessian・RT 回転の並列化で全体 1.600 倍、CPU 最適化前からの累積で 8.887 倍**。
前者は同じ Metal 設定でのコード変更比較、後者は旧 CPU 版と現在の Metal 版の比較。

## 条件と対象

- 入力：Hi-net tilt 水平成分、2004-01-01〜03、各試行20イベント・433採用窓。
- コンパイラ：Homebrew LLVM/Clang 23.1.1（C/C++）、Apple Clang（Objective-C++）。
  全版 Release `-O3 -DNDEBUG`、native arch 無効、16 OpenMP スレッド、dynamic teams 無効。
- GPU 設定：`AUTOFOCUSING_BACKEND=metal`、`AUTOFOCUSING_GPU_POWER=bootstrap`。
  初期グリッドは CPU。CUDA の既存比較と同じ power 設定を採用した。
- 同一マシン上で別プロセスを交互に実行。版ごとに最初のウォームアップ1回を除外。
  OS キャッシュは消去せず、計測中に他のエージェント起動ビルド・テストは実行しなかった。
- 検証用ドライバだけ FFTW を `ESTIMATE | UNALIGNED`、Bootstrap seed を
  `1837 + 104729 * replicate` に固定。production source は変更していない。
- wall time は起動、入力、GPU 転送・同期、出力を含む。

| 対象 | コミット | 用途 |
| --- | --- | --- |
| CPU 最適化前 | `393a93e030d74ace231dbd46ea85175a0795f59a` | 累積比較の基準。CPU のみ |
| Hessian・回転並列化前 | `6a7a4bd4d40a20f1383f683fe233776aed674efa` | 既存 Metal に共通 profiling を追加済み。数値カーネル変更なし |
| 現在の gpu | `dabb03727454421f4fbb0b42a61604539fe99a6d` | 並列化後。両比較の対象 |

`393a93e` は最初の CPU slant-stack 最適化より前であり、CUDA の累積比較と同じ基準。
`6a7a4bd` は CUDA の CPU 並列化比較と同じ基準。

## (i) Hessian・RT 回転並列化の効果

両版 Metal/Bootstrap、profiling 有効、各5回の中央値。

| 処理 | 並列化前 | 並列化後 | 高速化率 |
| --- | ---: | ---: | ---: |
| 全体 wall time | 21.039 s | 13.150 s | **1.600×** |
| Hessian | 8.945 s | 1.417 s | **6.311×** |
| RT 回転 | 0.342 s | 0.103 s | **3.312×** |
| Fitting | 10.083 s | 3.592 s | 2.807× |
| Bootstrap 段階 | 1.702 s | 0.660 s | 2.577× |
| 初期グリッド段階 | 1.144 s | 1.060 s | 1.080× |

全体の経過時間は **37.50% 短縮**。測定範囲は並列化前20.980〜21.105 s、
並列化後13.099〜13.430 s。全体と初期グリッドの回帰判定に合格した。
Hessian・回転はいずれも対象段階10%以上短縮を満たす。
Hessian、回転、Fitting は `serial_caller` の inclusive time であり、
段階間に重複があるため合算しない。
個々の変更を単独で有効化した実験ではなく、採用済み2変更の組合せでの内訳。

ウォームアップを含む12試行すべてで、イベント出力がバイト単位で一致。
採用窓・初期候補も一致した。比較では Bootstrap の GPU 誤差許容を緩めず、
同一バックエンドの変更として照合している。

[全試行・ハッシュ・段階時間](benchmarks/metal-cpu-parallel-20260914/parallel/report.json)
と、その隣に各試行のログ・イベントファイルを保存した。

## (ii) CPU 最適化前からの累積効果

旧版 CPU/off と現在版 Metal/Bootstrap、両方 profiling 無効、各3回の中央値。

| 対象 | 全体中央値 | 測定範囲 |
| --- | ---: | ---: |
| CPU 最適化前 `393a93e` | 121.214 s | 120.622〜121.527 s |
| 現在版 Metal＋CPU 最適化 | 13.639 s | 13.443〜13.660 s |

**8.887倍、経過時間88.75%短縮**。CPU slant stack、共通 I/O、Metal 化、
今回の Hessian・回転並列化を含む、同一マシンでの直接比較。
異なる測定の倍率を掛け合わせた値ではない。
(i) とは試行・profiling 条件が異なるため、現在版の中央値も別に報告する。

ウォームアップを含む8試行で20イベント、433採用窓、初期候補が一致。
各版のイベント出力は反復間でバイト単位に一致した。
CPU/GPU 間は既存の許容誤差を適用し、全38列を検証した。

- Beam 最大値・MAD（12〜13列）：`rtol=1e-4`、最大相対差 `1.5593e-5`。
- Bootstrap（21〜33列）：`rtol=1e-3`、最大相対差 `2.8718e-5`。
- 収束残差（15列）：`atol=1e-12`。
- その他：出力文字列一致。水平モードの既存 U 欠測 NaN は保持し、
  数値比較列の非有限値は許可しない。

[全試行・ハッシュ・精度・採用窓](benchmarks/metal-cpu-parallel-20260914/cumulative/report.json)
と、その隣に各試行のログ・イベントファイルを保存した。

## ビルド・再現手順

各コミットを別の完全なチェックアウトに用意する。ビルド設定は
[先行する Metal 動作検証](gpu-metal-verification-20260914.md)と同じ。
Objective-C++ を Apple Clang に指定する。旧 CPU 版には Metal のビルド設定は不要。
`build_io_probe.py` は現在版のスクリプトを旧版の `tests/` にもコピーして使用できる。
基準版の他のソースや数値カーネルを現在版に置き換えない。

```bash
# 各チェックアウトでビルド後に実行
python3 tests/build_io_probe.py build

# gpu チェックアウトのルートから、存在しない出力ディレクトリを指定
python3 tests/run_cpu_parallel_events.py \
  "$BEFORE_PROBE" "$CURRENT_PROBE" "$HINET_ROOT" "$CMT_CATALOG" \
  parallel-events --backend metal --power bootstrap --threads 16 --days 3 --repeats 5

# この検証で追加した再利用可能な比較スクリプト
python3 tests/run_legacy_gpu_events.py \
  "$LEGACY_PROBE" "$CURRENT_PROBE" "$HINET_ROOT" "$CMT_CATALOG" \
  cumulative-events --backend metal --power bootstrap --threads 16 --repeats 3
```

3版の CTest は旧版2/2、並列化前7/7、現在版7/7に合格。
現在版の Metal slant-stack 32条件・power 216条件とバッチ境界テストも合格。
最初の旧版ランチャーテストは archive コピーに Git 情報がないため失敗したが、
対象コミットの Git 情報を付けた再実行で合格した。ソース修正は不要だった。
Metal と `/usr/bin/time -l` はこのホストのエージェント sandbox 外で実行した。
[ビルド・CTest ログ](benchmarks/metal-cpu-parallel-20260914/build-checks/after/ctest.log)
も保存している。比較スクリプトは同一出力・Bootstrap 境界・イベント識別不一致・
非有限値拒否の人工例で確認した。

既存 CUDA 記録は `gpu@dabb037:docs/cpu-fitting-parallel.md` にある。
CUDA は3990X＋RTX PRO 2000上で並列化1.445倍、CPU最適化前から5.683倍。
今回の1.600倍・8.887倍は M4 Max 内の比較であり、CUDA/Metal の直接性能比較ではない。
実3成分データ、5日間の並列化前後比較、全アーカイブ、1スレッド性能の追加測定は
今回の範囲外で、すべての main マージ条件を完了したという意味ではない。

作業ブランチ `cuda` は維持し、製品コードの編集、インストール、コミット、push はしていない。
