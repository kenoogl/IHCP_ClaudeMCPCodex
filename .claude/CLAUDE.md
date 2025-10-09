# CLAUDE.md

逆熱伝導問題（IHCP）をCGMで解くJulia実装プロジェクト。

## 状態
✅ Python→Julia移植完了（Phase 1-6、505テスト全合格）
✅ マトリクスフリーPBICGSTAB!実装（89倍高速化）
**ブランチ**: tuning2 | **更新**: 2025-10-08

## 実行

### Julia版
```bash
# 全テスト（505項目）
julia --project=. -e 'using Pkg; Pkg.test()'

# メイン実行
julia --project=. src/main.jl --config config/test.toml
```

### Phase構成（全6完了）
| Phase | 内容 | テスト |
|-------|------|-------|
| 1 | 熱物性値 | 25 |
| 2 | DHCP（マトリクスフリー版） | 298 |
| 3 | Adjoint | 13 |
| 4 | CGM | 18 |
| 5 | スライディングウィンドウ | 29 |
| 6 | 検証・実データ読込 | 122 |

### 依存
- Julia 1.10+
- コア: LinearAlgebra, FLoops
- I/O: NPZ, JLD2, MAT, CSV
- 詳細: `julia/Project.toml`

## アーキテクチャ

### 主要モジュール
```
julia/src/
├── solvers/
│   ├── CommonSolver.jl      # マトリクスフリーPBICGSTAB!
│   ├── DHCPSolver.jl        # 直接解法（311行、89倍高速）
│   ├── SensitivitySolver.jl # 感度問題
│   ├── AdjointSolver.jl     # 随伴解法
│   ├── CGMSolver.jl         # 共役勾配法
│   └── SlidingWindowSolver.jl
└── utils/
    ├── boundary_conditions.jl
    └── commons.jl           # WorkBuffers
```

### コアアルゴリズム
1. **DHCP**: `solve_dhcp!()` - 前進時間差分、マトリクスフリーPBICGSTAB!
2. **Sensitivity**: `solve_sensitivity!()` - 熱流束微小変化に対する温度応答
3. **Adjoint**: `solve_adjoint!()` - 後退時間差分、勾配計算
4. **CGM**: `solve_cgm!()` - 最適化メインループ

## 規約

- **言語**: 日本語のみ、コメント詳細記述
- **インデント**: 2スペース
- **配列**: `(ni, nj, nk)`順序、1-indexed
- **TDD**: テスト先行開発

## 検証結果

### 完全一致達成 ✅
- 熱流束相対誤差: **最大0.0041%**
- 温度場相対誤差: **最大0.0047%**
- Python完全一致（356ステップ検証）

### 性能
| 版 | 計算時間 | 高速化 |
|----|---------|--------|
| cg!版（旧） | 53.2秒 | 1x |
| **マトリクスフリー版** | **0.6秒** | **89x** |

## 重要事項
- **データ**: `T_measure_700um_1ms.npy`（1.1GB、gitignore済み）必須
- **メモリ**: 8GB以上推奨
- **数値精度**: 性能改善時は必ず検証

## ドキュメント
- `docs/INDEX.md`: 全ドキュメント索引
- `docs/logs/`: Phase 1-6実装ログ
- `shared/results/validation/`: 検証結果