# 共通 CPU I/O 最適化

## 実装（2026-09-12）

`perf/cpu-io` は CPU 最適化済み `main` (`51f60a7`) を基点とする。
フィルタ並列化は `ab175dc`、HDF5 情報の共有は `46dd8e0`。

- HDF5 は従来どおり単一スレッドで読み込む。読み込み完了後、観測点・成分ごとに独立した `hp_filt` を OpenMP で処理する。
- `load_h5` が日単位の初期化・読み込みをまとめ、ファイルハンドル、観測点列挙、成分情報を共有する。
- 観測点名の検索を索引化し、観測点配列を事前確保して直接構築する。
- 既存の `init_station` / `read_h5` API を維持する。単位変換、フィルタ係数・演算精度、観測点順序、選択半径、成分 QC は維持する。
- 座標属性の欠落・非有限値は、従来の未定義な座標使用に代わり明示的な例外とする。ファイルは例外時にも解放する。

日をまたぐ先読み、HDF5 の並列呼び出し、GPU 処理、キャッシュ容量変更は含まない。
`OMP_NUM_THREADS` がフィルタにも適用される。`AUTOFOCUSING_PROFILE=1` で
`#LOAD_PROFILE init_s=... read_decode_copy_s=... filter_s=... total_s=...` を出力する。
`read_decode_copy_s` は HDF5 読み込み・伸長・変換・コピー・残りのメタデータ処理の合計で、純粋な SSD 待ち時間ではない。
従来の `#Read data: init_station/read_h5` 行は `#Read data: load_h5` にまとまる。

## 検証結果

CPU ビルド、Metal 有効ビルドとも CTest 4 件が成功。
`tests/io_reference.h` は `51f60a7` のローダーを固定した比較用実装で、
水平2成分・3成分の合成データについて 1・4・16 スレッドで観測点順序・座標・
採用数・波形メタデータ・保存波形全サンプルが一致（波形はビット一致）。
従来 API、繰り返し読み込み、欠落成分、不正単位・サンプリング、非直交軸、
存在しないファイル、半径内に観測点がない場合、欠落座標と例外時のハンドル解放を確認。

Metal 検証ブランチ `test/metal-cpu-io` は共通 I/O の上に旧 Metal コミット
`27fccd3` を再適用したもの。Apple M4 Max の実 GPU テストも 32 ケース成功し、
ピーク不一致は 0、最大 scaled error は約 `5.05e-6`。
これは既存 Metal 計算の回帰確認で、I/O 高速化率の測定ではない。

SSD 接続後、2004-01-01～03、2014-01-01、2024-01-01 の全波形比較も成功。
採用観測点数は順に 653・653・652・725・725。各日 1・4・16 スレッド、
繰り返し読み込みと従来 API を含めて確認した。
[実波形検証記録](benchmarks/cpu-io-equivalence-20260912.json)。

2004-01-01 のキャッシュ済み入力を各 5 回、変更前後で交互に測定した中央値：

| スレッド数 | 変更前 | 変更後 | 高速化率 |
| --- | ---: | ---: | ---: |
| 1 | 2.119 s | 2.080 s | 1.02 倍 |
| 4 | 2.124 s | 1.298 s | 1.64 倍 |
| 8 | 2.117 s | 1.144 s | 1.85 倍 |
| 16 | 2.121 s | 1.080 s | 1.96 倍 |

16 スレッドでは入力時間を約 49% 削減。全試行の実ディスク読み込みカウンタは
0 byte で、USB 転送帯域や未キャッシュ時の測定ではない。
ピーク RSS の最大値は変更前約 1.035 GB、変更後約 1.011 GB。
中央値・範囲・CPU 時間・RSS の生データは [計測 JSON](benchmarks/cpu-io-2004001.json) を参照。

### 全体計算の比較

FFTW 計画と bootstrap seed を検証用ビルドで固定し、2004-01-01～03 を
変更前後それぞれ 2 回、実行順を逆にして測定した。バックエンドごとに、
20 イベントの全出力バイト（bootstrap を含む）と 433 の採用窓がすべて一致。

| バックエンド | 変更前の全体時間（中央値） | 変更後 | 短縮率 |
| --- | ---: | ---: | ---: |
| CPU、16 スレッド | 71.225 s | 68.325 s | 4.1% |
| Metal、CPU 16 スレッド | 45.210 s | 42.325 s | 6.4% |

複数スレッドの入力時間 10% 以上改善、全体時間の悪化 5% 以内という統合条件を満たした。
これは I/O 最適化前との比較であり、CPU slant stack 最適化以前との比較ではない。
通常ビルドは従来の FFTW_MEASURE と時刻由来の bootstrap seed を維持する。

- [CPU 計測・照合記録](benchmarks/cpu-io-events-cpu-20260912.json)、[CPU イベント出力](benchmarks/cpu-io-events-cpu-20260912.dat)
- [Metal 計測・照合記録](benchmarks/cpu-io-events-metal-20260912.json)、[Metal イベント出力](benchmarks/cpu-io-events-metal-20260912.dat)


## 再現方法

Homebrew LLVM を C/C++、Apple Clang を Metal の Objective-C++ に使用する。
後者を明示すると、この環境の LLVM 23 で生じた Objective-C メソッド呼び出しのリンクエラーを回避できる。

```bash
cmake -S . -B build-io -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_C_COMPILER="$(brew --prefix llvm)/bin/clang" \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++"
cmake --build build-io -j4
OMP_NUM_THREADS=4 ctest --test-dir build-io --output-on-failure
```

Metal 検証ブランチでは configure に以下を追加し、別ビルドディレクトリを使う。

```bash
-DDELTAP_ENABLE_METAL=ON -DCMAKE_OBJCXX_COMPILER="$(xcrun -f clang++)"
# ビルド後、GPU アクセス可能な環境で実行
OMP_NUM_THREADS=4 ./build-metal-io/test_metal_slant_stack
```

SSD 接続後の計測例（`test_io` は両ブランチにある）：

```bash
python3 tests/benchmark_io.py build-io/test_io \
  /Volumes/Seismic_Data/hdf5/Hi-net_tilt/2004/0101/20040010000.h5 \
  docs/benchmarks/cpu-io-2004001.json
```

スクリプトは先に全波形比較を実行し、キャッシュが温まった条件で
1・4・8・16 スレッド、各実装 5 回を交互に測定する。
wall time の中央値・範囲・比率、CPU 時間、実ディスク読み込みバイト数、
プロセスのピーク RSS を JSON に保存する。macOS の RSS は byte 単位。
各試行は別プロセスで、FFTW 初期化条件も両実装で同一。
これはコールドキャッシュや USB の帯域測定ではない。

## イベント出力の再現手順

`CMAKE_EXPORT_COMPILE_COMMANDS=ON` で configure し、通常ビルドした後、
テスト用ドライバだけ FFTW 計画と bootstrap seed を固定する。
CPU の変更前コミットは `51f60a7`、Metal では `27fccd3` を指定する。
同じバックエンドのビルドディレクトリを使用すること。

```bash
python3 tests/build_io_probe.py build-io --reference 51f60a7
python3 tests/build_io_probe.py build-io
python3 tests/run_io_events.py \
  build-io/io-fixed-reference/cal_ccf_io_probe \
  build-io/io-fixed-current/cal_ccf_io_probe \
  /Volumes/Seismic_Data/hdf5/Hi-net_tilt \
  ../moment_loc_76_24 build-io/io-events-cpu
```

Metal ではビルドパスを置き換え、実行に `--backend metal` を追加する。
出力先は未作成のディレクトリを指定する。3 日分だけの symlink 入力を作り、
前→後、後→前の順で各 2 回実行する。イベント数 20、セグメント数 12、
採用窓数と bootstrap を含む出力全バイトの一致、全体時間の悪化 5% 以内を確認する。
実行ログ・イベントファイル・`report.json` は指定出力先に残る。
重要な集計は `docs/benchmarks/` に保存する。

## TODO：翌日分の先読み（保留、2026-09-12）

現状の I/O 性能で十分なため、ユーザー判断で追加の先読み最適化は棚上げする。
現状はキャッシュ済み入力・前処理が約 1.08 秒/日（16 スレッド）。
NAS 利用、データ量の増加、計算部分のさらなる高速化によって入力待ちが
全体時間の無視できない割合になった時点で、再計測して優先度を判断する。

現行の `prefetchFile()` は Ubuntu/Linux 向けの
`posix_fadvise(..., POSIX_FADV_WILLNEED)` による OS への先読み依頼。
現在の macOS ビルドでは条件コンパイルでこの呼び出しが除外され、
翌日のファイルを開閉するだけになっている。CPU/I/O 最適化で削除したものではない。
macOS にも `fcntl(..., F_RDADVISE, ...)` はあるが、現行コードには未実装。
どちらの OS 向け依頼も、翌日の波形の展開・フィルタ完了を管理する仕組みとは異なる。

- [ ] 再開時は、当日の解析と翌日の読み込み・前処理を重ねる、CPU・Metal 共通の仕組みを検討する。
- [ ] RAM が十分な現環境では、翌日分 1 件を RAM に保持する方式を第一候補とする。
  容量は 1 GB 固定にせず実データに応じて確保し、保持件数・メモリ上限を管理する。
  現データは圧縮ファイル約 0.5 GB/日、展開後の水平波形約 0.9 GB/日で、別途作業領域が必要。
- [ ] 観測点中心座標などの日ごとの共有状態を分離し、HDF5 の非スレッドセーフな呼び出しを並行させない。
  先読み側の CPU 使用量も制限・計測し、当日の解析との競合を確認する。
- [ ] 連続する複数日で、波形・イベント出力の一致、初日と定常時の時間、ピーク RAM を比較する。
  NAS ではキャッシュ条件を区別して測定する。

ローカル SSD への一時保存は、RAM 制約や後日の再解析で再転送を避けたい場合の代案。
RAM に余裕がある現在の条件で、RAM 先読みより高速だと判断したわけではない。
この TODO は検討記録のみで、先読みのコード変更は行っていない。

## ブランチ統合

実データの検証条件を満たした共通 I/O コミットを `main` に fast-forward し、
その先頭の上に既存 Metal 開発コミットを再適用する。
`perf/cpu-io` は `main` と同じ先頭、`test/metal-cpu-io` は `metal` と同じ先頭に揃える。
旧 `metal` (`27fccd3`) は `backup/metal-before-cpu-io-20260912` に保存する。
push・バイナリのインストールは行わない。

```text
51f60a7  CPU slant stack 最適化
  └─ 共通 CPU I/O 最適化・検証  (main, perf/cpu-io)
       └─ Metal backend       (metal, test/metal-cpu-io)
```

調査資料は元の `investigate/metal-float-candidates` の `dd8bc8c` に保存済み。
GPU/FP32 調査専用コードは CPU ブランチへ取り込んでいない。
