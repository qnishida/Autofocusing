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

## 検証結果と未完了項目

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

外付け SSD `/Volumes/Seismic_Data` が未接続のため、実データの波形比較・
性能測定・イベント出力比較は未実施。小規模合成データで計測スクリプトの動作のみ確認した。
実データに対する高速化率はまだ確定していない。**`main` と `metal` は未更新**。

## 再現方法

Homebrew LLVM を C/C++、Apple Clang を Metal の Objective-C++ に使用する。
後者を明示すると、この環境の LLVM 23 で生じた Objective-C メソッド呼び出しのリンクエラーを回避できる。

```bash
cmake -S . -B build-io \
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

## 統合前の残作業

1. 2004-01-01～03 と 2014/2024-01-01 の実波形を `test_io --real FILE` で比較。
2. 上記の反復測定で、複数スレッドの入力時間が 10% 以上改善することを確認。
3. FFTW 計画と bootstrap の乱数を検証用ビルドで固定し、2004-01-01～03 の 20 イベントを各バックエンドの変更前後で比較。全体時間の悪化が 5% 以内であることを確認。生産コードの乱数処理は変更しない。
4. 成功後 `main` を `perf/cpu-io` へ fast-forward。旧 `metal` をバックアップ参照で保存し、共通 I/O 後の `main` を親とする検証済み Metal 先頭へ更新。push・インストールは行わない。

調査資料は元の `investigate/metal-float-candidates` の `dd8bc8c` に保存済み。
GPU/FP32 調査専用コードは CPU ブランチへ取り込んでいない。
