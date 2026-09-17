# パラメータ共分散の計算精度（2026-09-17）

[ピーク探索の検証](peak-search-validation-20260917_jp.md)で見つかったNaNへの対応。
変更は `485d761` を基点とする `test/peak-search-boundaries` の変更に含める。
main、リモート、インストール済みバイナリは更新していない。

本資料は共分散のみの精度変更を記録する。後続の
[Newton法の倍精度化](newton-precision-20260917_jp.md)は別に検証する。
以下の「Newton法は単精度」という記述は、この先行する比較時点を指す。

## 原因と修正

追加されたT候補の1行で、パラメータ共分散の10列すべてがNaNになった。
対象は2005-03-16 11:51:28開始の窓、初期値 (0.063, -0.054) s/km。
乱数を固定した再実行で元の行を再現した。
H/sigmaの各要素は有限値だが、行列式は約5.205e39で、floatの最大値約3.403e38を超える。
単精度の4×4逆行列は非有限値になり、倍精度では有限値になった。

`est_dist_boot()` から共分散への変換を `set_parameter_covariance()` に切り出し、
H/sigmaの作成と逆行列計算を倍精度にした。出力規約は従来どおり。

```text
C_output(i,j) = -inverse(H/sigma)(i,j) * W(i) * W(j)
W = (0.06, pi/2, pi/2, 0.04/(30*111))
```

Wのfloatとしての表現値も維持した。位置誤差への変換時にJacobianでWを打ち消す
既存の規約を保つため、公開済みの位置誤差を0.06倍する変更ではない。
この精度修正では、bootstrapの抽出方法・sigma・推定の最適化・収束条件・
ピーク選択・スペクトル行列の計算は変更していない。

## 回帰テスト

`tests/parameter_covariance.cpp` は本番の共分散変換関数を呼ぶ。
非対角項を持つ通常のHessianに2通りのsigmaを与えた場合、実測のH/sigma、
実際のsigmaを用いて戻した生のHessianを検査する。
同じ入力に対して旧float逆行列が非有限値になることも確認する。
独立した、スケールを調整したlong doubleのLU解法と、規格化・非対角項を含め全要素を比較した。
このApple arm64環境ではlong doubleとdoubleの仮数は共に53ビット。
追加の有効桁ではなく、スケール調整と別の解法による独立検証である。
最大相対差は実測例で3.81e-12、通常のHessianで7e-16未満。
関連ターゲットの再ビルド後、選択したCTest 12件がすべて通った。

## 実データ検証

2005年3月15〜17日をCPUで再計算した。0.1〜0.25 Hz、±0.165 s/km、dp=0.005 s/km、
8スレッド。ピーク探索のみを修正した実行と、bootstrapの固定乱数・FFTW設定・
品質判定の履歴をそろえた。入力波形・参照カタログのハッシュも一致する。

- 初期候補147件、出力143行で同一。最適化の全記録、反復回数、終了理由も同一。
- 全143行が有限値になり、既存の検証条件を緩めずにすべて通過した。
- 共分散21〜30列を除く28列は、出力文字列として完全一致。
  位置、slowness、max/MAD、bootstrapのパワー・sigma、SRR/STT/SUU・交差項を含む。
  エネルギー比と後段の選別は変わらない。
- 選別された75行も同じで、R=31、T=10、U=34。
- 以前から有限値だった共分散要素の相対変化 `abs(new-old)/abs(old)` は、
  中央値0、95パーセンタイル5.90e-6、最大5.55e-5。
  選別された行に限ると最大3.80e-5（0.00380%）。これは有効数字6桁の出力値による
  比較で、内部精度での差ではない。以前から有限だった要素の符号反転はない。

NaNが解消したT行には負の対角要素があり、有限値でも誤差としては使えない。
max/MADは約-2.709で、従来の>7条件を満たさない。
別の選別対象外のR行にも修正前から負の対角要素があり、その状態は変わらない。
選別された75行には負の対角要素はないが、それだけで誤差の被覆率や正定値性が
検証されたという意味ではない。

## 結果の適用範囲

今回修正するのは確認された数値オーバーフローであり、誤差の統計モデルではない。
特異・悪条件のHessianやsigmaがゼロの場合も含め、倍精度なら常に有効な共分散になる
という保証ではない。正則化、イベントの除外、非有限値の置換は追加していない。

## 再現方法

```bash
cmake -S . -B build-clang -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-clang --target test_parameter_covariance --parallel 4
ctest --test-dir build-clang -R '^parameter_covariance_range$' --output-on-failure
```

`tests/build_covariance_probe.py` は保存済みの乱数固定ピーク探索検証用ソースから、
bootstrap関数と共分散変換関数を現在の本番コードに置き換える。乱数の固定は維持し、
元の関数が基準コードと一致することと、ピーク探索が同じであることを確認する。
`tests/compare_covariance_runs.py` は入力ハッシュ、品質判定後の窓、初期候補、最適化の記録、
選別結果、共分散以外の全出力列を比較する。修正後には既存の検証条件の全通過を要求する。

保存済みのピーク探索検証環境を準備した後のローカル実行例：

```bash
python3 -B tests/build_covariance_probe.py build-clang
python3 -B tests/run_matrix_audit_events.py \
  build-clang/covariance-double-probe/cal_ccf_covariance \
  /Volumes/Seismic_Data/hdf5/Hi-net ../moment_loc_76_24 \
  build-clang/covariance-double-mar15-17 --start 2005-03-15 --days 3 --threads 8
python3 -B tests/compare_covariance_runs.py \
  build-clang/peak-search-mar15-17-fixed build-clang/covariance-double-mar15-17 \
  build-clang/covariance-double-comparison.json
```

詳細カタログ、ログ、検証用生成ソース、ビルドコマンド、ハッシュは
`docs/benchmarks/parameter-covariance-precision-20260917/` にローカル保存し、Git管理外とする。
元のfloat計算の失敗記録は、別のピーク探索検証ディレクトリに残している。

## 補足：行列の区別・符号・CPUの計算精度

3×3のスペクトル行列 S_RTU と、4×4のパラメータ共分散は別の出力である。
S_RTUから観測点の自己項を引く操作は半正定値性を保証せず、補正後のパワー推定値が
負になることはある。一方、パラメータ共分散は -sigma W H^-1 W から求め、S_RTUを
逆行列化しているわけではない。1回のHessian評価で固定されるスペクトル・重みに
対して、スカラーの自己項は4つの位相パラメータに依存せず、2階微分はゼロ。
減算はbootstrapのsigmaに影響し得るが、正のスカラーsigmaは共分散の定値性を反転させない。

NaNだったT例の記録済みH/sigmaをDecimalの60桁でLDL分解すると、対角ピボットは
約 (-2.01770756e8, +4.83715177e6, -6.83665288e11, +7.80057330e12)。
合同変換による符号の保存から、正・負の固有値が2つずつあり、単精度で逆行列化する
前のHessian自体が不定である。有限値に戻ることは、極大点の有効な共分散になることを
意味しない。最適化中の曲率判定を通った理由は未切り分け。コード上は最終更新前に
曲率を判定し、最適化後にR/Tの回転を更新してからbootstrap用Hessianを評価するため、
評価点・入力の違いは今後の確認対象となる。

CPU処理も混合精度である。目的関数とHessianの集計は倍精度だが、Newton法の
固有値分解・更新ベクトルは現在も単精度。共分散の逆行列も修正前はCPU上の単精度で、
今回そのCPU処理を倍精度にした。確認したNaNはCPUのみの実行で発生している。
今回の比較はCPU対GPUの比較ではなく、残る単精度Newton処理を全倍精度版と比較して
検証したものでもない。
