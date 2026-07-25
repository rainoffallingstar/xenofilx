# xenofilx 实现状态

## 已完成的工作

### 1. 项目架构 ✅
- 完整的 Go 项目结构
- 使用 cobra 的 CLI 框架
- 模块化的内部包设计

### 2. 核心算法实现 ✅
- **编辑距离计算**: NM tag + CIGAR 插入(I) + 软剪裁(S)
- **单端分类**: 比较两个参考基因组上的编辑距离
- **双端分类**: 成对 reads 的平均分数比较
- **配置管理**: 参数验证和默认值

### 3. BAM I/O 纯 Go 实现 ✅
- `internal/bamnative/` - 纯 Go 实现的 BAM 读写
- 支持 BAM 压缩/解压缩
- 支持 BAM 索引构建
- 支持 BAM 排序检查

### 4. NM Tag 计算模块 ✅
- `internal/bamnative/nmtag.go` - NM 计算
- 支持与参考基因组比对计算 NM
- 支持 bisulfite 模式 (C→T 不计入 mismatch)

### 5. FASTA 索引功能 ✅
- `internal/bamnative/faidx.go` - FASTA 索引
- 自动检测 gzip 压缩
- 自动生成 .fai 索引文件
- LRU 缓存优化内存使用

### 6. CLI 接口 ✅
```
xenofilx run \
  --graft <files> \
  --host <files> \
  --output <dir> \
  [--mm-threshold 4] \
  [--unmapped-penalty 8] \
  [--threads 1] \
  [--graft-ref <fa>] \
  [--host-ref <fa>] \
  [--recalculate-nm] \
  [--bisulfite]
```

### 7. 结果统计表格 ✅
- Total: 总 reads 数
- GraftOnly: 仅比对到 graft
- HostOnly: 仅比对到 host
- Both: 比对到两者
- Graft(%): 分类为 graft 比例
- Host(%): 分类为 host 比例
- Discard(%): 超出阈值比例
- TotalGraft(%): 总 graft / 总 reads
- Thresh: 使用的阈值

## 项目文件结构

```
xenofilx/
├── cmd/xenofilx/main.go        # CLI 入口
├── internal/
│   ├── bamnative/              # BAM I/O 纯 Go 实现
│   │   ├── bamnative.go        # BAM 读写核心
│   │   ├── nmtag.go           # NM tag 计算
│   │   ├── faidx.go           # FASTA 索引
│   │   ├── sort.go            # BAM 排序
│   │   └── index.go          # BAM 索引
│   ├── classifier/             # 分类算法
│   │   ├── classifier.go
│   │   ├── edit_distance.go
│   │   ├── single_end.go
│   │   └── paired_end.go
│   ├── filter/                  # 过滤逻辑
│   │   └── filter.go
│   └── config/                 # 配置管理
│       └── config.go
├── pkg/cli/                     # CLI 命令
│   └── run.go
├── go.mod                      # 依赖管理
├── README.md                   # 文档
└── STATUS.md                   # 状态
```

## 功能验证

- ✅ BAM 读取 (纯 Go)
- ✅ BAM 写入 (BGZF 压缩)
- ✅ BAM 排序检查
- ✅ BAM 自动排序
- ✅ BAM 索引构建
- ✅ 过滤管道
- ✅ NM tag 计算
- ✅ FASTA 索引生成
- ✅ gzip 压缩 FASTA 支持
- ✅ 结果统计表格

## 性能考虑

- **纯 Go 实现**: 无外部依赖
- **内存使用**: 比 R 版本低
- **流式处理**: FASTA 按需加载 + LRU 缓存
- **并发**: 使用 goroutines 实现多样本并行处理

## 下一步计划

1. 添加更多测试用例
2. 性能优化
3. 跨平台二进制发布
