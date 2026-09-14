# Metal: 5日間の厳密比較と1スレッド性能検証

2026-09-14、`6a7a4bd` と `gpu@7de4c22` を比較し、依頼された2項目はともに合格した。
**5日分の出力はバイト単位で一致し、1スレッドでは全体・Hessian・RT回転の
中央値および各測定ペアのすべてが5%以内の増加だった。**
`main` への統合は実施していない。

## 共通条件

| 項目 | 条件 |
| --- | --- |
| 並列化前 | `6a7a4bd4d40a20f1383f683fe233776aed674efa` |
| 現在版 | `7de4c2261d509a65f037f481270ffd579471451c` |
| ハードウェア | Apple M4 Max |
| OS | macOS 26.6.2、build 25G83 |
| C/C++ | Homebrew LLVM/Clang 23.1.1 |
| Objective-C++ | Apple Clang 17.0.0 (`clang-1700.6.3.2`) |
| ビルド | Release `-O3 -DNDEBUG`、native arch 無効、Metal 有効、CUDA 無効 |
| Backend / power | `metal` / `bootstrap`、初期グリッドは CPU |
| 入力 | Hi-net tilt、水平成分 |
| FFTW | 検証用ドライバで `ESTIMATE \| UNALIGNED` |
| Bootstrap seed | 検証用ドライバで `1837 + 104729 * replicate` |
| 計測 | `AUTOFOCUSING_PROFILE=1`、`OMP_DYNAMIC=FALSE`、別プロセス |

各版を独立したチェックアウトから同じ設定で新規ビルドした。
対象2版の Metal カーネルには差分がない。固定 seed と FFTW 条件は既存の
`tests/build_io_probe.py` で生成した検証用ドライバだけに適用し、製品ソースは変更していない。
両版の CTest は7/7に合格。
[基準版ログ](benchmarks/metal-main-qualification-20260914/build-checks/before/ctest.log)・
[現在版ログ](benchmarks/metal-main-qualification-20260914/build-checks/after/ctest.log)を保存した。

## 1. 5日分の並列化前後比較

2004-01-01〜05、16スレッド、各版1回。
精度確認用なのでウォームアップは省略し、この試行の時間を性能採用の測定値には使わない。

- 両版で29イベント、20セグメント、740採用窓が一致。
- セグメントごとの採用窓数と、58行の初期候補ログ（29イベント×2）が一致。
- 全38列を、Bootstrap列を含めて出力精度で一致させる既存の同一バックエンド判定を使用。
  収束残差（15列）だけは既存の `atol=1e-12` を維持。
  **GPU用の `rtol=1e-3` / `rtol=1e-4` には緩めていない。**
- 実際には残差も含めてイベントファイル全体がバイト単位で一致した。
- 共通イベント SHA-256：
  `92893ed8efa1d1a8ee428dd7dda0211f452ce9cdc16d7e797d6a2e87de4b0ec7`。

[機械可読レポート](benchmarks/metal-main-qualification-20260914/five-day/report.json)に
全試行、バイナリ・イベントハッシュ、採用窓、初期候補を保存。
同じディレクトリにログとイベントファイルも保存している。

## 2. 1スレッド性能の比較

既存の性能測定と同じ2004-01-01〜03を使用。
`OMP_NUM_THREADS=1`、各版1回のウォームアップを除外し、5組を交互の順序で測定した。
順序は基準→現在、現在→基準、基準→現在、現在→基準、基準→現在。
測定中に他のエージェント起動ビルド・テスト・ベンチマークは実行していない。
OSキャッシュは消去せず、全体時間は起動・I/O・GPU転送・同期・出力を含む。

| 項目 | 基準版中央値 | 現在版中央値 | 中央値の増減 | 各ペアの増減範囲 | 判定 |
| --- | ---: | ---: | ---: | ---: | --- |
| 全体 wall time | 36.30185 s | 36.36624 s | +0.177% | +0.022〜+0.623% | 合格 |
| Hessian | 8.92417 s | 8.97119 s | +0.527% | −0.359〜+0.694% | 合格 |
| RT回転 | 0.352154 s | 0.362434 s | +2.919% | +1.649〜+3.332% | 合格 |

中央値の比と5組のペアごとの比を別々に検査した。
**3項目すべてで5%超のペアは0/5**。この測定では、再現性のある5%超の性能低下は見られない。
回転には約3%の小さな増加があり、ゼロオーバーヘッドと主張する結果ではない。
Hessianと回転は `serial_caller` の inclusive time で、全体 wall time とは別の指標。

全体時間の範囲は基準版36.26295〜36.47868 s、現在版36.32404〜36.56101 s。
初期グリッドの中央値も10.27492 → 10.29204 sで、既存の5%回帰判定に合格した。
ウォームアップを含む12試行すべてで20イベント・433採用窓・初期候補が一致し、
Bootstrapを含む全イベント出力もバイト単位で一致した。

[機械可読レポート](benchmarks/metal-main-qualification-20260914/single-thread/report.json)の
`regression` に中央値、各ペアの増減、範囲、5%判定を保存。
同じディレクトリに全試行のログとイベントファイルも保存している。
保存後に14イベントファイルのハッシュ、件数、採用窓、初期候補を再確認した。

## 再現手順

各チェックアウトで、既存の [Apple Clang 併用設定](gpu-metal-verification-20260914.md)により
ビルドし、`python3 tests/build_io_probe.py build` で固定条件のドライバを作る。
`BEFORE_PROBE` は `6a7a4bd`、`CURRENT_PROBE` は `7de4c22` のドライバを指す。
以下は現在版チェックアウトのルートから実行する。出力ディレクトリは未作成のパスを使う。

```bash
python3 tests/run_cpu_parallel_events.py \
  "$BEFORE_PROBE" "$CURRENT_PROBE" "$HINET_ROOT" "$CMT_CATALOG" \
  five-day --backend metal --power bootstrap --threads 16 \
  --days 5 --repeats 1 --warmups 0

python3 tests/run_cpu_parallel_events.py \
  "$BEFORE_PROBE" "$CURRENT_PROBE" "$HINET_ROOT" "$CMT_CATALOG" \
  single-thread --backend metal --power bootstrap --threads 1 \
  --days 3 --repeats 5 --warmups 1
```

Metal GPU と `/usr/bin/time -l` にアクセスするため、エージェント sandbox 外で実行した。
使用した既存ドライバは出力の同一性を厳密に検査し、保存レポートにその結果を記録する。
Hessian・回転の5%判定は各試行の `cpu` 値から集計した。
比較式は `100 * (current / reference - 1)`、合格条件は `current <= 1.05 * reference`。

この記録は依頼された水平データでの2項目を完了するもの。
実3成分・全アーカイブ・別の機種への一般化は行っていない。
以前の [16スレッド高速化・CPU最適化前との比較](metal-cpu-parallel-performance-20260914.md)
と合わせてレビューできる。`main` への統合はユーザーによる結果確認後の別作業とする。
