# Mac と hiyodori の入力・計算速度を比較する計画

作成日：2026-10-22。測定前の計画であり、以下のMac側の測定は未実施。
測定結果のチェックポイントは `153b7ad`。
[hiyodoriの測定記録](hiyodori-performance-20261022_jp.md)と
[入力ハッシュ・測定要約JSON](hiyodori-performance-20261022.json)を比較の基準にする。

## 目的と判断方法

MacのUSB4 SSDで速かった理由を、次の三つに分けて確認する。

1. ファイルへのアクセス：同じMacでUSB4 SSDと同じNASを比較する。
2. HDF5処理：ファイル転送以外の展開・変換・コピー・前処理の差を調べる。
3. 実運用：同じ入力・計算条件でCPU/Metalとhiyodoriの結果を比較する。

MacのSSDとhiyodoriのNASだけを比べるとCPU・OS・ライブラリの差も入るため、
ストレージが原因とは断定しない。
最初にMac内でSSD/NASを比較し、その後、同じ手順をhiyodoriのNAS/ローカルSSDで行う。
Macの既存の約1.080 s/日は別入力のキャッシュ済み測定であり、今回の基準値にはしない。

## 固定する入力と条件

最初の比較は1日分に限定する。Macの過去のtilt水平2成分へ置き換えず、
hiyodoriで使用した **Denoiseの3成分** を用いる。

| 項目 | 固定値 |
| --- | --- |
| 波形の相対パス | `2004/0401/20040920000.h5` |
| ファイルサイズ | 897,926,923 byte |
| SHA-256 | `3431f694f21a40e3661f52cc0ace3ed41429c51d6b8b38abe792a5b5f2003407` |
| 成分 | `3c` |
| 日付 | 2004-04-01のみ |
| 周波数 | 要求0.1–0.25 Hz、実効下限0.099609375 Hz |
| Slant stack | 67×67、刻み0.005 s/km、各軸最大0.165 s/km |
| イベント選択 | `selected`、R/T/Uのmax/MAD閾値3/3/15 |
| 検証用乱数 | `1837u + 104729u * i` |
| 検証用FFTW | `ESTIMATE | UNALIGNED` |
| 主比較のOMP数 | 16、`OMP_DYNAMIC=FALSE` |
| アフィニティ | OS既定。`OMP_PROC_BIND`・`OMP_PLACES`は指定しない |

カタログと速度モデルも要約JSONのSHA-256と一致させる。
同一日の1日実行では既存コードの翌日先読みは呼ばれず、翌日ファイルへのアクセスが混ざらない。

Mac側で次を実際の保存先に設定する。NASはhiyodoriと同じデータを持つ共有を用いる。

```bash
export SSD_HINET_ROOT="/Volumes/<USB4-SSD>/<Denoiseのルート>"
export NAS_HINET_ROOT="/Volumes/<NASマウント>/<Denoiseのルート>"
export CMT_CATALOG="<moment_loc_76_24のパス>"
export COMPARE_OUT="<ワークスペース>/results/mac-storage-comparison/<UTC実行ID>"
```

ルート配下は両方とも `2004/0401/` に対象ファイルを1個だけ置く。
コピー元・先のサイズとSHA-256を検証し、パスだけが異なる同一ファイルであることを確認する。
コピーとハッシュ計算はキャッシュを温めるため、直後の実行を「コールド」と呼ばない。

## 準備：環境記録とMac向けビルド

1. 未コミットのユーザー変更を退避・破棄せず、`153b7ad` の独立したクリーンなチェックアウトを使う。
   このコミットはhiyodoriの測定時から文書・測定要約だけを追加しており、計算ソースは同じ。
2. Macの機種・チップ、物理/論理コア数、P/Eコア構成、RAM、macOS、電源モードを記録する。
   SSDの機種・ファイルシステム・実際のUSB接続状態、NASのNFS/SMB種別・マウントオプション、
   ネットワークの実際の接続速度も保存する。Macが10GbE接続とは決めつけない。
3. C/C++コンパイラ、OpenMP、HDF5、圧縮フィルタと関連ライブラリ、FFTWの版とリンク先を保存する。
   hiyodoriと異なる版を使用した場合は明記する。比較のために依存環境を移行・更新しない。
4. まずCPUのみのReleaseビルドを別ディレクトリに作る。既存環境のHomebrew LLVMを使う例：

```bash
cmake -S . -B build-mac-io-compare \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DDELTAP_ENABLE_NATIVE_ARCH=OFF \
  -DDELTAP_ENABLE_CUDA=OFF -DDELTAP_ENABLE_METAL=OFF \
  -DCMAKE_C_COMPILER="$(brew --prefix llvm)/bin/clang" \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++"
cmake --build build-mac-io-compare -j4
python3 tests/build_io_probe.py build-mac-io-compare
```

既存の依存関係の準備方法は[ビルド手順](../README.md)を優先する。
Linuxの実行ファイルをMacへコピーして使わず、同じソースをMac向けにコンパイルする。
Mac内のSSD/NAS比較には、同一のMac実行ファイルを使用する。
ビルド・既存のI/O/周波数/イベント選択の関連チェックを終えてから計測を開始する。
計測中はビルド、ハッシュ計算、他のベンチマークを並行させない。

## 測定A：単純なファイル読み取り

- Python標準ライブラリ等の小さな測定ドライバで、対象ファイルを先頭から末尾まで読む。
  固定8 MiBバッファの `readinto()` を使い、open→EOFまでの読み取り→closeを計時する。
  計時区間でハッシュ計算・HDF5処理・解析をしない。
  全ファイルをRAM上の巨大な戻り値として保持しない。
- 各ストレージについて計測対象外の全読み取りを1回行い、その後、新規プロセスで各5回測定する。
  順序はSSD→NAS、NAS→SSDを交互にする。
- wall time、プロセスCPU時間、実際に読んだbyte数、byte/wall timeを記録する。
  読み取りbyte数は常に897,926,923 byteであることを確認する。
- NASのネットワーク通信量はNIC/NASクライアント統計の前後差も保存する。
  これは他通信やメタデータ通信を含み得る参考値であり、プロセス専用の転送量とは呼ばない。
  ローカルSSDのdisk I/OカウンタとNASのネットワーク転送量を混同しない。

この手順の主結果は **キャッシュを温めた読み取り性能**。
実ストレージ/ネットワークの転送が観測されない試行から、USB4や10GbEの帯域を算出しない。
初回アクセスを追加記録する場合も、キャッシュ状態を保証できなければ
「初回・キャッシュ状態未確認」と分類する。OS全体のキャッシュ消去や再起動は本計画では行わない。
未キャッシュの持続転送性能は、この結果だけでは未測定として残す。

## 測定B：実際の3成分HDF5入力処理

専用のC++ローダーを新設せず、既存の固定乱数ドライバを1日分実行し、
`#LOAD_PROFILE` を集計する。比較の主バックエンドは両ストレージともCPUに固定する。
既存 `test_io --benchmark` は水平成分固定、`benchmark_io.py` は旧/新ローダー比較なので、
その出力を今回の3成分ストレージ比較の代わりには使わない。

1回の実行例（チェックアウトルートから、毎回別の `TAG` を使う）：

```bash
export BINARY="$PWD/build-mac-io-compare/io-fixed-current/cal_ccf_io_probe"
export INPUT_ROOT="$SSD_HINET_ROOT"  # NAS試行ではNAS_HINET_ROOT
export TAG="ssd-t16-warmup"          # 試行ごとに一意にする
OMP_NUM_THREADS=16 OMP_DYNAMIC=FALSE \
AUTOFOCUSING_BACKEND=cpu AUTOFOCUSING_GPU_POWER=off \
AUTOFOCUSING_EVENT_SELECTION=selected \
AUTOFOCUSING_MIN_MAX_MAD_R=3 AUTOFOCUSING_MIN_MAX_MAD_T=3 AUTOFOCUSING_MIN_MAX_MAD_U=15 \
AUTOFOCUSING_FREQ_MIN=0.1 AUTOFOCUSING_FREQ_MAX=0.25 \
AUTOFOCUSING_SLOWNESS_STEP=0.005 AUTOFOCUSING_SLOWNESS_MAX=0.165 \
AUTOFOCUSING_PROFILE=1 \
/usr/bin/time -l "$BINARY" 2004 mac-storage "$TAG" "$INPUT_ROOT" "$CMT_CATALOG" \
  "$COMPARE_OUT/results" 3c 2004-04-01 2004-04-01 \
  > "$COMPARE_OUT/$TAG.log" 2>&1
```

実行前に出力ディレクトリを新規作成し、継承されたverify系・CUDA profiling・
Metal旧設定・OMPアフィニティ変数を除去する。OS既定のアフィニティで統一する。
profile変数は上記の通常の段階計測だけを有効にする。
ワークスペースのランチャーはmachine-local設定で環境変数を上書きし得るため、
比較では上記のように検証用実行ファイルを直接起動する。

- OMP16でSSD/NASそれぞれ慣らし1回＋新規プロセス5回。順序は交互にする。
- `init_s`、`read_decode_copy_s`、`filter_s`、`total_s`、全体wall time、
  プロセスCPU時間、ピークRSS、受理観測点・成分・窓・初期候補・イベント数を保存する。
  MacのRSSはbyte、Linuxの`ru_maxrss`はKiBという単位差を正規化する。
- 入力以外の計算も実行されるが、ストレージ要因の主判定は `#LOAD_PROFILE` で行う。
  `read_decode_copy_s` を純粋なI/O待ちとは呼ばず、HDF5展開・変換・コピー等を含むと明記する。
- 現在のCUDA設定の短い検証は5イベントだったが、CPUバックエンドへ5を無条件に流用しない。
  最初のCPU基準試行で件数・受理窓・初期候補を確立し、同じMacのSSD/NAS全試行で一致を確認する。
  固定乱数の全38列は既存 `compare_power(..., bootstrap=False)` の規則で比較し、
  バイト単位一致の有無も別に記録する。許容誤差を新設・緩和しない。

同じ測定A/BをhiyodoriのNAS/ローカルSSDでも、同じソース・同じ入力・CPU/offで実施する。
既存のhiyodori CUDA/Bootstrapの全体10.44 sを、Mac CPU/offの全体時間の基準には流用しない。
現在のhiyodoriのOMP32での `read_decode_copy_s=6.208530 s` は参考として併記し、
  正式な機種比較にはOMP16の同条件の新規CPU測定を使う。

## 測定C：実運用のCPU/Metal比較

入力処理の結果を確定してから実施する。Macでは同じソースの別ディレクトリにMetalを有効化し、
既存手順どおりApple ClangをObjective-C++コンパイラとして指定する。
詳細は[Metal測定のビルド条件](metal-cpu-parallel-performance-20260914_jp.md)と
[I/O検証のビルド手順](cpu-io-optimization_jp.md)を参照。
Metal回帰チェック後、同じ検証用FFTW・乱数設定のドライバを生成する。

- 同一SSD入力、同じMetal対応ドライバ、OMP16でCPU/offとMetal/offを交互に各5回、各モードの慣らし1回を除外する。
  Slant stackだけをGPU化する効果を測る。入力・全体・Slant stack・Bootstrap・初期探索を分ける。
- 3成分のMetal Bootstrapは現行ではCPU OpenMP経路であり、CUDAのGPU FP64 Bootstrapと同等ではない。
  Metal/bootstrapを「3成分BootstrapもGPU」と解釈しない。
- Mac内のCPU/GPU出力比較は既存のGPU比較規則を使う。機種間の丸め差が既存基準を超えた場合は
  不一致を記録して調査し、速度比較のために基準を緩和しない。
- 同じ1日分で入力差と計算差を説明できれば終了する。
  12年分や14日分の再計算は、このMac側の切り分けには含めない。

## 判定・保存する成果物

| 観測 | 解釈と次の対象 |
| --- | --- |
| 同一MacでNASのraw読み取りもHDF5入力も遅い | データアクセスの寄与あり。キャッシュと実転送量も併記する |
| raw読み取りは同程度、HDF5入力だけNASが遅い | HDF5のアクセス粒度・メタデータ・NFS/SMB往復などを調査する |
| Mac内のSSD/NAS差は小さいがMac/hiyodoriのHDF5入力差は大きい | CPU、HDF5/圧縮ライブラリ、変換・コピーなどの候補を調査する |
| raw/HDF5入力差は小さいが全体時間が異なる | GPUバックエンド、Bootstrap経路、探索・Hessianなど計算側を調べる |

上記は切り分けであり、単一のタイマーだけで原因を確定しない。
全5回の最小–最大と、試行順を反転した場合の傾向も示す。

成果物は日時付きの新規 `results/mac-storage-comparison/<UTC実行ID>/` に保存する。
環境・ソース/実行ファイル/入力/カタログ/速度モデルのハッシュ、実行設定、
各回のログ・イベント出力・検証結果、段階別CSV/JSON、中央値・最小・最大・時間比を残す。
大きなログや外部波形はGitへ追加せず、確定後に要約だけを `docs/` に記録する。
Mac実測値がそろうまではUSB4 SSDやNASを主因と断定せず、未実施の条件も明記する。
