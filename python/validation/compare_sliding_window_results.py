#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python-Juliaスライディングウィンドウ結果比較スクリプト

両言語で実行したスライディングウィンドウ計算の結果を読み込み、
性能・精度を比較してMarkdownレポートを生成する。

実行方法:
  python python/validation/compare_sliding_window_results.py --phase 1
  python python/validation/compare_sliding_window_results.py --phase 2 --python-result shared/results/python_sliding_window_cgm20000.npz

出力:
  shared/results/sliding_window_comparison_cgm{N}.md
"""

import numpy as np
import argparse
from pathlib import Path
from datetime import datetime
import sys

# ---------------------------------------------------------------------------
# データ読み込み
# ---------------------------------------------------------------------------

def load_result_file(npz_path):
  """
  npzファイルを読み込み、結果データを返す。

  Parameters
  ----------
  npz_path : Path
      npzファイルパス

  Returns
  -------
  dict
      読み込んだデータの辞書
  """
  if not npz_path.exists():
    raise FileNotFoundError(f"結果ファイルが見つかりません: {npz_path}")

  data = np.load(npz_path, allow_pickle=True)

  # windows_infoは配列形式で保存されているので取得
  windows_info = data['windows_info'].item() if data['windows_info'].ndim == 0 else data['windows_info']

  result = {
    'q_global': data['q_global'],
    'T_final': data['T_final'],
    'elapsed_total': float(data['elapsed_total']),
    'windows_info': windows_info,
    'cgm_max_iter': int(data['cgm_max_iter']),
    'nt': int(data['nt']),
    'window_size': int(data['window_size']),
    'overlap': int(data['overlap']),
  }

  return result

# ---------------------------------------------------------------------------
# 性能比較
# ---------------------------------------------------------------------------

def compare_performance(python_data, julia_data):
  """
  性能指標を比較する。

  Parameters
  ----------
  python_data, julia_data : dict
      読み込んだ結果データ

  Returns
  -------
  dict
      性能比較結果
  """
  py_total = python_data['elapsed_total']
  jl_total = julia_data['elapsed_total']

  py_windows = python_data['windows_info']
  jl_windows = julia_data['windows_info']

  # ウィンドウ数
  n_py_windows = len(py_windows)
  n_jl_windows = len(jl_windows)

  # ウィンドウ平均時間
  py_avg_window = np.mean([w['elapsed_window'] for w in py_windows])
  jl_avg_window = np.mean([w['elapsed_window'] for w in jl_windows])

  # 総CGM反復数
  py_total_cgm = sum(w['cgm_iterations'] for w in py_windows)
  jl_total_cgm = sum(w['cgm_iterations'] for w in jl_windows)

  perf = {
    'python_total_time': py_total,
    'julia_total_time': jl_total,
    'julia_speedup': py_total / jl_total,
    'julia_speedup_percent': (1.0 - jl_total / py_total) * 100,
    'python_n_windows': n_py_windows,
    'julia_n_windows': n_jl_windows,
    'python_avg_window_time': py_avg_window,
    'julia_avg_window_time': jl_avg_window,
    'python_total_cgm_iters': py_total_cgm,
    'julia_total_cgm_iters': jl_total_cgm,
  }

  return perf

# ---------------------------------------------------------------------------
# 精度比較
# ---------------------------------------------------------------------------

def compare_accuracy(python_data, julia_data):
  """
  精度を比較する。

  Parameters
  ----------
  python_data, julia_data : dict
      読み込んだ結果データ

  Returns
  -------
  dict
      精度比較結果
  """
  q_py = python_data['q_global']  # (ni, nj, nt)
  q_jl = julia_data['q_global']

  T_py = python_data['T_final']  # (ni, nj, nk)
  T_jl = julia_data['T_final']

  # 配列形状確認
  if q_py.shape != q_jl.shape:
    print(f"警告: q_global形状不一致 - Python: {q_py.shape}, Julia: {q_jl.shape}")
  if T_py.shape != T_jl.shape:
    print(f"警告: T_final形状不一致 - Python: {T_py.shape}, Julia: {T_jl.shape}")

  # 熱流束誤差
  q_diff = np.abs(q_jl - q_py)
  q_rel_diff = q_diff / (np.abs(q_py) + 1e-10)  # ゼロ除算回避

  # 温度場誤差
  T_diff = np.abs(T_jl - T_py)
  T_rel_diff = T_diff / (np.abs(T_py) + 1e-10)

  acc = {
    # 熱流束
    'q_abs_max': np.max(q_diff),
    'q_abs_mean': np.mean(q_diff),
    'q_abs_rms': np.sqrt(np.mean(q_diff**2)),
    'q_rel_max': np.max(q_rel_diff) * 100,  # %
    'q_rel_mean': np.mean(q_rel_diff) * 100,
    'q_rel_rms': np.sqrt(np.mean(q_rel_diff**2)) * 100,
    # 温度場
    'T_abs_max': np.max(T_diff),
    'T_abs_mean': np.mean(T_diff),
    'T_abs_rms': np.sqrt(np.mean(T_diff**2)),
    'T_rel_max': np.max(T_rel_diff) * 100,
    'T_rel_mean': np.mean(T_rel_diff) * 100,
    'T_rel_rms': np.sqrt(np.mean(T_rel_diff**2)) * 100,
  }

  return acc

# ---------------------------------------------------------------------------
# ウィンドウ別比較
# ---------------------------------------------------------------------------

def compare_windows(python_data, julia_data):
  """
  ウィンドウ別の詳細比較を行う。

  Parameters
  ----------
  python_data, julia_data : dict
      読み込んだ結果データ

  Returns
  -------
  list of dict
      ウィンドウごとの比較結果
  """
  py_windows = python_data['windows_info']
  jl_windows = julia_data['windows_info']

  n_windows = min(len(py_windows), len(jl_windows))

  window_comparisons = []
  for i in range(n_windows):
    py_w = py_windows[i]
    jl_w = jl_windows[i]

    comp = {
      'window_id': py_w['window_id'],
      'start_idx': py_w['start_idx'],
      'end_idx': py_w['end_idx'],
      # CGM反復数
      'py_cgm_iters': py_w['cgm_iterations'],
      'jl_cgm_iters': jl_w['cgm_iterations'],
      'cgm_iters_diff': jl_w['cgm_iterations'] - py_w['cgm_iterations'],
      # 実行時間
      'py_time': py_w['elapsed_window'],
      'jl_time': jl_w['elapsed_window'],
      'time_speedup': py_w['elapsed_window'] / jl_w['elapsed_window'],
      # 熱流束統計
      'py_q_mean': py_w['q_mean'],
      'jl_q_mean': jl_w['q_mean'],
      'q_mean_diff': abs(jl_w['q_mean'] - py_w['q_mean']),
    }

    window_comparisons.append(comp)

  return window_comparisons

# ---------------------------------------------------------------------------
# レポート生成
# ---------------------------------------------------------------------------

def generate_markdown_report(python_data, julia_data, perf, acc, window_comps, output_path):
  """
  比較結果をMarkdownレポートとして出力する。

  Parameters
  ----------
  python_data, julia_data : dict
      読み込んだ結果データ
  perf : dict
      性能比較結果
  acc : dict
      精度比較結果
  window_comps : list of dict
      ウィンドウ別比較結果
  output_path : Path
      出力ファイルパス
  """
  lines = []

  # ヘッダー
  lines.append("# Python-Julia スライディングウィンドウ比較レポート")
  lines.append("")
  lines.append(f"**生成日時**: {datetime.now().strftime('%Y年%m月%d日 %H:%M:%S')}")
  lines.append("")
  lines.append("---")
  lines.append("")

  # 1. 実行環境
  lines.append("## 1. 実行環境")
  lines.append("")
  lines.append(f"- **CGM反復数**: {python_data['cgm_max_iter']}")
  lines.append(f"- **時間ステップ数**: {python_data['nt']}")
  lines.append(f"- **ウィンドウサイズ**: {python_data['window_size']}")
  lines.append(f"- **オーバーラップ**: {python_data['overlap']}")
  lines.append(f"- **実行ウィンドウ数**: Python={perf['python_n_windows']}, Julia={perf['julia_n_windows']}")
  lines.append("")

  # 2. 性能比較
  lines.append("## 2. 性能比較")
  lines.append("")
  lines.append("### 2.1 総実行時間")
  lines.append("")
  lines.append("| 項目 | Python | Julia | Julia/Python比 |")
  lines.append("|------|--------|-------|----------------|")
  lines.append(f"| 総実行時間 | {perf['python_total_time']:.2f}s | {perf['julia_total_time']:.2f}s | {perf['julia_speedup']:.2f}x ({perf['julia_speedup_percent']:+.1f}%) |")
  lines.append(f"| 平均ウィンドウ時間 | {perf['python_avg_window_time']:.2f}s | {perf['julia_avg_window_time']:.2f}s | {perf['python_avg_window_time']/perf['julia_avg_window_time']:.2f}x |")
  lines.append(f"| 総CGM反復数 | {perf['python_total_cgm_iters']} | {perf['julia_total_cgm_iters']} | {perf['julia_total_cgm_iters']-perf['python_total_cgm_iters']:+d} |")
  lines.append("")

  # 3. 精度比較
  lines.append("## 3. 精度比較")
  lines.append("")
  lines.append("### 3.1 熱流束誤差 (q)")
  lines.append("")
  lines.append("| 誤差指標 | 絶対誤差 [W/m²] | 相対誤差 [%] |")
  lines.append("|----------|----------------|--------------|")
  lines.append(f"| 最大値 | {acc['q_abs_max']:.3e} | {acc['q_rel_max']:.3f} |")
  lines.append(f"| 平均値 | {acc['q_abs_mean']:.3e} | {acc['q_rel_mean']:.3f} |")
  lines.append(f"| RMS | {acc['q_abs_rms']:.3e} | {acc['q_rel_rms']:.3f} |")
  lines.append("")

  lines.append("### 3.2 温度場誤差 (T_final)")
  lines.append("")
  lines.append("| 誤差指標 | 絶対誤差 [K] | 相対誤差 [%] |")
  lines.append("|----------|--------------|--------------|")
  lines.append(f"| 最大値 | {acc['T_abs_max']:.3e} | {acc['T_rel_max']:.3f} |")
  lines.append(f"| 平均値 | {acc['T_abs_mean']:.3e} | {acc['T_rel_mean']:.3f} |")
  lines.append(f"| RMS | {acc['T_abs_rms']:.3e} | {acc['T_rel_rms']:.3f} |")
  lines.append("")

  # 4. ウィンドウ別比較
  lines.append("## 4. ウィンドウ別詳細比較")
  lines.append("")
  lines.append("| ID | 範囲 | CGM反復数 (Py) | CGM反復数 (Jl) | 差分 | 時間 (Py) [s] | 時間 (Jl) [s] | 高速化率 |")
  lines.append("|----|------|----------------|----------------|------|---------------|---------------|----------|")

  for comp in window_comps:
    lines.append(
      f"| {comp['window_id']} "
      f"| [{comp['start_idx']}, {comp['end_idx']}] "
      f"| {comp['py_cgm_iters']} "
      f"| {comp['jl_cgm_iters']} "
      f"| {comp['cgm_iters_diff']:+d} "
      f"| {comp['py_time']:.2f} "
      f"| {comp['jl_time']:.2f} "
      f"| {comp['time_speedup']:.2f}x |"
    )

  lines.append("")

  # 5. 結論
  lines.append("## 5. 結論")
  lines.append("")

  # 成功基準判定
  phase = 1 if python_data['cgm_max_iter'] <= 10 else 2

  if phase == 1:
    lines.append("### Phase 1 成功基準判定")
    lines.append("")
    check_q_rel = acc['q_rel_max'] < 5.0
    check_windows = perf['python_n_windows'] == perf['julia_n_windows']
    lines.append(f"- [ ] 熱流束相対誤差 < 5%: **{acc['q_rel_max']:.2f}%** {'✅' if check_q_rel else '❌'}")
    lines.append(f"- [ ] ウィンドウ数の一致: Python={perf['python_n_windows']}, Julia={perf['julia_n_windows']} {'✅' if check_windows else '❌'}")
    lines.append("")
  else:
    lines.append("### Phase 2 成功基準判定")
    lines.append("")
    check_T_rel = acc['T_rel_max'] < 2.0
    check_q_rel = acc['q_rel_max'] < 1.0
    check_speedup = perf['julia_speedup'] >= 2.0
    lines.append(f"- [ ] 温度場相対誤差 < 2%: **{acc['T_rel_max']:.2f}%** {'✅' if check_T_rel else '❌'}")
    lines.append(f"- [ ] 熱流束相対誤差 < 1%: **{acc['q_rel_max']:.2f}%** {'✅' if check_q_rel else '❌'}")
    lines.append(f"- [ ] Julia高速化率 >= 2.0x: **{perf['julia_speedup']:.2f}x** {'✅' if check_speedup else '❌'}")
    lines.append("")

  lines.append("### サマリー")
  lines.append("")
  lines.append(f"- Julia版は **{perf['julia_speedup']:.2f}倍** の高速化を達成")
  lines.append(f"- 熱流束の相対誤差は最大 **{acc['q_rel_max']:.2f}%**")
  lines.append(f"- 温度場の相対誤差は最大 **{acc['T_rel_max']:.2f}%**")
  lines.append("")

  lines.append("---")
  lines.append("")
  lines.append("**生成スクリプト**: `python/validation/compare_sliding_window_results.py`")
  lines.append("")

  # ファイル書き込み
  output_path.parent.mkdir(parents=True, exist_ok=True)
  with open(output_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

  print(f"✅ レポート生成完了: {output_path}")

# ---------------------------------------------------------------------------
# メイン処理
# ---------------------------------------------------------------------------

def main():
  parser = argparse.ArgumentParser(
    description='Python-Juliaスライディングウィンドウ結果比較',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
使用例:
  python python/validation/compare_sliding_window_results.py --phase 1
  python python/validation/compare_sliding_window_results.py --phase 2 --python-result custom_py.npz
    """
  )

  parser.add_argument(
    '--phase',
    type=int,
    choices=[1, 2],
    required=True,
    help='Phase番号 (1: CGM 3回, 2: CGM 20000回)'
  )
  parser.add_argument(
    '--python-result',
    type=str,
    default=None,
    help='Pythonの結果ファイル (デフォルト: shared/results/python_sliding_window_cgm{N}.npz)'
  )
  parser.add_argument(
    '--julia-result',
    type=str,
    default=None,
    help='Juliaの結果ファイル (デフォルト: shared/results/julia_sliding_window_cgm{N}.npz)'
  )
  parser.add_argument(
    '--output-report',
    type=str,
    default=None,
    help='出力レポート (デフォルト: shared/results/sliding_window_comparison_cgm{N}.md)'
  )

  args = parser.parse_args()

  # ファイル名決定
  cgm_iter_str = "3" if args.phase == 1 else "20000"

  if args.python_result is None:
    python_result_path = Path(f"shared/results/python_sliding_window_cgm{cgm_iter_str}.npz")
  else:
    python_result_path = Path(args.python_result)

  if args.julia_result is None:
    julia_result_path = Path(f"shared/results/julia_sliding_window_cgm{cgm_iter_str}.npz")
  else:
    julia_result_path = Path(args.julia_result)

  if args.output_report is None:
    output_report_path = Path(f"shared/results/sliding_window_comparison_cgm{cgm_iter_str}.md")
  else:
    output_report_path = Path(args.output_report)

  print("=" * 80)
  print("Python-Julia スライディングウィンドウ結果比較")
  print("=" * 80)
  print(f"Phase: {args.phase}")
  print(f"Python結果: {python_result_path}")
  print(f"Julia結果: {julia_result_path}")
  print(f"出力レポート: {output_report_path}")
  print()

  # データ読み込み
  print("データ読み込み中...")
  try:
    python_data = load_result_file(python_result_path)
    print(f"  ✅ Python結果読み込み完了")
  except Exception as e:
    print(f"  ❌ Python結果読み込み失敗: {e}")
    return 1

  try:
    julia_data = load_result_file(julia_result_path)
    print(f"  ✅ Julia結果読み込み完了")
  except Exception as e:
    print(f"  ❌ Julia結果読み込み失敗: {e}")
    return 1

  print()

  # 性能比較
  print("性能比較中...")
  perf = compare_performance(python_data, julia_data)
  print(f"  ✅ 完了 (Julia高速化率: {perf['julia_speedup']:.2f}x)")
  print()

  # 精度比較
  print("精度比較中...")
  acc = compare_accuracy(python_data, julia_data)
  print(f"  ✅ 完了 (q相対誤差最大: {acc['q_rel_max']:.2f}%, T相対誤差最大: {acc['T_rel_max']:.2f}%)")
  print()

  # ウィンドウ別比較
  print("ウィンドウ別比較中...")
  window_comps = compare_windows(python_data, julia_data)
  print(f"  ✅ 完了 ({len(window_comps)}ウィンドウ)")
  print()

  # レポート生成
  print("レポート生成中...")
  generate_markdown_report(python_data, julia_data, perf, acc, window_comps, output_report_path)
  print()

  print("=" * 80)
  print("比較完了")
  print("=" * 80)

  return 0

if __name__ == '__main__':
  sys.exit(main())
