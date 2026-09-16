# Metal パワー評価の最適化

詳細なbenchmarkログと実行結果はローカルに保存しており、公開リポジトリには
含めない。本文中の詳細記録への言及は、そのローカル資料を指す。

[English](metal-power-optimization_en.md)

## 対象と使い方

共通 I/O 改善済み `metal` (`f745119`) を基点に、`perf/metal-power` で
Bootstrap と初期 fitting グリッドの目的関数を試験実装した。
slant stack の `(px, py)` グリッドとは別の処理である。

```bash
export AUTOFOCUSING_BACKEND=metal
export AUTOFOCUSING_METAL_POWER=bootstrap  # off / bootstrap / grid / all
export OMP_NUM_THREADS=16
```

既定値は `off`。`bootstrap` は100標本と最後の補正済みパワー、
`grid` は距離35候補・曲率40候補、`all` は両方を Metal で評価する。
CPU バックエンドでは常に既存の倍精度経路を使う。
実3成分データによる追加経路の検証は未実施なので、新経路は水平モードに限定し、
3成分モードは既存の CPU 目的関数を使う。既存の3成分 Metal slant stack は維持する。

最終 fitting の目的関数・勾配・Hessian、乱数生成方式、標本数、周辺の統計計算、
FFT、I/O は変更しない。通常実行の bootstrap は従来どおり時刻由来の seed を使うため、
異なる実行の bootstrap 値の完全一致は要求しない。比較は検証用の固定 seed で行う。

## 実装

`metal_power_batch` はパラメータ群、スペクトル、観測点座標、重みを受け取る内部 API。
GPU 内部は FP32、CPU への返却型は double。返却型が double でも GPU の誤差は残る。

- スペクトルは1バッチ呼び出しにつき1回準備し、GPU バッファを再利用する。
- 最大32候補ずつ評価し、拡張グリッドの一時 GPU メモリを制限する。
- Bootstrap の位相は標本間で共有する（各32候補の塊で生成）。
- 各周波数・時間窓で観測点を加算し、パワーと自己項補正を別々に階層的に集計する。
  Metal の fast math は無効にする。
- グリッドは GPU 最大値との差が `2e-3 * abs(max) + 1e-30` 以内の候補を
  CPU 倍精度で並列再評価し、最後の選択を元の順序・厳密な大小比較で行う。
  この幅は合成試験の最大誤差より余裕を持たせた工学的な基準で、一般的な誤差保証ではない。
- `AUTOFOCUSING_PROFILE=1` で `#POWER_PROFILE stage=grid|bootstrap total_s=...` を出力する。
  準備・転送・同期・再評価を含み、Bootstrap 時間には既存の Hessian も含む。

## 精度検証

- Metal 有効・無効の両ビルドで既存 CTest 4件が成功。
- 既存 Metal slant stack の実 GPU 回帰テストが成功。
- 新カーネルの合成216条件（2/31/650観測点、1/3窓、5/90/175度、
  ランダム・コヒーレント波形、振幅 `1e-12/1/1e6`、マスク・非一様重み、補正あり/なし）が成功。
  最大誤差は未補正の倍精度パワーに対して `1.73e-4`。判定基準は `5e-4`。
- 65個の異なるパラメータ・重みで32候補の境界をまたぐ添字を照合。
  ゼロ信号・同点候補、CPU/3成分での無効化も検証した。
- 2004-01-01～05 の既存 Metal 版29イベントを基準に、off、Bootstrap、grid、all を照合。
  後半2日は候補再評価幅の調整に使わない独立検証日とした。
  最初の3日は20イベント、後半2日は9イベント。
  以前の CPU 調査の30イベントと混同しない。
- すべての初期候補ログとイベント識別・件数、非 Bootstrap 列が一致。
  収束残差の許容絶対誤差は `1e-12`。Bootstrap 関連列21～33は
  `rtol=1e-3, atol=1e-30` で合格し、最大相対差は約 `2.87e-5`（0.0029%）。
  水平モードに元々ある欠測 U の NaN は維持し、新しい非有限値は認めない。
- 既存の `compare_event_results.py` は Bootstrap 列を除外するため、この判定には使わない。
  `run_metal_power.py` が該当列も比較する。

## 性能・採用結果（2026-09-12）

Apple M4 Max、CPU 16スレッド。2004-01-01～03 の20イベントを各経路5回、
実行順を交互にして測定した。各試行は別プロセスで初期化も含み、
OS・ドライバのキャッシュは消去していない。中央値は以下のとおり。

| 経路 | 全体時間 | Bootstrap 段階 | グリッド段階 |
| --- | ---: | ---: | ---: |
| off | 41.895 s | 22.912 s | 1.076 s |
| bootstrap | 20.641 s | 1.670 s | 1.068 s |
| grid | 41.484 s | 22.743 s | 0.859 s |
| all | 20.417 s | 1.668 s | 0.845 s |

両方有効で全体は **2.05倍**。個別の段階は Bootstrap
13.72倍、グリッド 1.25倍。
対象段階10%以上短縮・全体時間の悪化5%以内の採用条件を両候補とも満たした。
実イベントのグリッド選択・非 Bootstrap 出力は元の Metal 版と一致し、
Bootstrap 関連値の最大相対差は全反復で約0.0029%以内だった。

- off の最大プロセス RSS：2.970 GB。
- bootstrap の最大プロセス RSS：3.158 GB。
- grid の最大プロセス RSS：3.194 GB。
- all の最大プロセス RSS：3.193 GB。

これはプロセス全体のピークで、GPU メモリ単独の値ではない。
速度・範囲・CPU 時間・メモリの詳細は 反復測定 JSON。

合成データのグリッド単体では、標準75候補は0.89～0.97倍で GPU が少し遅い。
750候補では1.61～1.91倍、7500候補では2.38～2.45倍だった。
全条件で CPU と同じ候補を選択した。GPU がすべての小さい入力で速いとは限らない。
この拡張評価は合成データであり、実データの探索条件変更時には別途検証する。
拡張グリッド JSON。

両候補を採用対象とし、`metal` へ統合する。既定値 `off` は維持する。
`main`・先読み TODO・インストール済みバイナリは変更せず、push は行わない。

## 検証記録

- カーネル試験
- 5日分の初期照合
- 倍精度再評価の並列化後の照合
- 元の Metal 版29イベント

初期照合の時間は併走ビルドの影響を含み得るため、採用の性能判定には使わない。
数値的な変更を伴わない倍精度再評価の並列化後に、独立検証日も再照合済み。

## 再現

```bash
cmake -S . -B build-power -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_C_COMPILER="$(brew --prefix llvm)/bin/clang" \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++" \
  -DCMAKE_OBJCXX_COMPILER="$(xcrun -f clang++)" \
  -DDELTAP_ENABLE_METAL=ON
cmake --build build-power -j4
OMP_NUM_THREADS=4 ctest --test-dir build-power --output-on-failure
OMP_NUM_THREADS=4 build-power/test_metal_power
OMP_NUM_THREADS=4 build-power/test_metal_slant_stack
build-power/test_metal_power --benchmark-grid
python3 tests/build_io_probe.py build-power --reference f745119
python3 tests/build_io_probe.py build-power
python3 tests/run_metal_power.py \
  build-power/io-fixed-current/cal_ccf_io_probe "$HINET_ROOT" "$CMT_CATALOG" \
  build-power/qualification-new --days 5 --repeats 1 \
  --legacy build-power/io-fixed-reference/cal_ccf_io_probe
python3 tests/run_metal_power.py \
  build-power/io-fixed-current/cal_ccf_io_probe "$HINET_ROOT" "$CMT_CATALOG" \
  build-power/benchmark-new --days 3 --repeats 5
```

実 GPU テストには GPU アクセスが必要。出力先は新しいディレクトリを指定する。
固定 FFTW 計画・seed はテスト用ドライバだけに適用する。
計測は他のビルド・性能テストと同時実行しない。
`--benchmark-grid` は650観測点・3窓で、候補数1/10/100倍を各5回測定する。
元の範囲の細分化と範囲拡大（距離1～179度、曲率 -2e-5～4e-5 の上端を除く）を分け、
2段階探索の候補選択を CPU/GPU で照合する。通常運用の探索範囲・刻みは変更しない。
