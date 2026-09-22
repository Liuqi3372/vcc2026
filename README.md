# VCC2026：最初版 Atlas Transfer 模型

用于预测单细胞 CRISPRi 扰动响应的 Atlas Transfer 模型。此代码包来自 `casia1:/ssd2/liuqi/xiexiu/code` 及同级最早的数据准备脚本，未纳入后续 `v1`、`v2`、`v3` 或网页评测适配版本。

本次只整理目录、移除 Python 注释与文档字符串、修正迁移后的文件路径，并补充说明。模型计算逻辑、超参数、随机种子和源数据处理公式均未调整。原有命令行帮助文本已转为普通字符串保留。

## 版本与评分

源代码 Git 提交为 `ffd0c66bad1c7a57153bda2da91a842f3ba58ef2`，取回时工作区干净。6 个核心文件的 SHA256 与当时保存的运行记录一致。

最初目录保存的官方评分 `score_avg` 为 **0.15592419253448622**，提交名为 `Xiexiu-X1-LocalAtlas-ABC-20260915`，提交 ID 为 `4HYR6o6tILZhdXbCoJru`。记录见 [official_latest.json](results/official_latest.json)。用户口述的 1.59 / 0.59 尚未匹配到评分证据，因此不作为本项目已验证的分数。

上游 README 另记载了 0.1545618019；它对应上游数据流程，不能与这里的本地数据适配结果混用。这里的“最初版”指 xiexiu 根目录中、后续 v1/v2/v3 之前保存的版本，不代表未经本地数据适配的上游实验。

## 目录

```text
vcc2026-original/
├── README.md
├── LICENSE
├── requirements.txt
├── requirements-pack.txt
├── src/
│   ├── model.py                  响应迁移、融合及整数计数生成
│   ├── predict.py                A/B/C 上下文推理入口
│   ├── prepare.py                上游公开数据准备流程
│   ├── sources.json              上游数据地址、大小及 SHA256
│   ├── pack.py                   VCC 提交包生成
│   └── compact.py                可选稀疏存储重排
├── scripts/
│   ├── prepare_original_local.py 本地 CD4 与 promoter 统计准备
│   ├── prepare_local_atlas.py     本地 K562/HCT116/HEK293T 统计准备
│   ├── run_local_sources.py       顺序运行本地数据准备
│   ├── audit_inputs.py            输入与源覆盖率校验
│   └── run_prediction_pipeline.py 原服务器历史流水线
├── docs/                         来源、数据与历史运行说明
├── results/                      历史评分和本次整理验证记录
├── tools/verify_delivery.py       交付文件完整性检查
├── data/                         预处理统计、官方控制细胞与基因面板
├── raw/                          原始数据
├── output/                       预测 H5AD 和 VCC 文件
└── logs/                         运行日志与状态
```

`data/`、`raw/`、`output/` 和 `logs/` 仅提供目录占位，不含大型数据、预测产物或运行环境。此模型通过源统计进行迁移，随包不提供另行训练的 checkpoint。自动提交和账号监控脚本未纳入模型交付包。

## 方法与固定参数

模型将 K562、HCT116、HEK293T 和 H1 四个来源的响应按 `[2, 1, 1, 2]` 融合，再引入 CD4 响应及 promoter 邻近基因校正。在目标上下文的控制细胞基础上，同时拟合 mean CPM 和 pseudobulk 目标，生成非负整数计数。

保留的参数包括随机种子 `20260910`、每个扰动 400 个细胞、`prior_counts=100000`、`log2fc` / `bulk_delta` 振幅 0.6 / 0.3、CD4 权重 0.5 和 promoter 校正系数 0.15。具体执行逻辑以 `src/model.py` 和 `src/predict.py` 为准。

## 环境

本地数据适配版历史推理环境为 Python 3.10；`requirements.txt` 使用当时记录的依赖版本。上游依赖文件单独保存在 `docs/requirements.upstream.txt`，不要将两个环境混装。以下命令面向 Linux / WSL，在项目根目录执行：

```bash
python3.10 -m venv env
env/bin/python -m pip install -r requirements.txt
python3.11 -m venv vcc-env
vcc-env/bin/python -m pip install -r requirements-pack.txt
```

打包使用独立 Python 3.11 环境和 `vcc-cli==0.2.0`。本次未重新安装这两个环境或执行完整推理。桌面副本可用于阅读和保存代码；`prepare_local_atlas.py` 使用 Linux 的 `fcntl`，不能直接在原生 Windows 环境执行。

## 输入数据

将最初版本的已准备文件放入 `data/`：

```text
gene_names.csv
pert_counts.csv
context_A.h5ad
context_B.h5ad
context_C.h5ad
K562_GWPS_CPM_full_statistics.npz
HCT116_full_statistics.npz
HEK293T_full_statistics.npz
H1_2025_full_statistics.npz
CD4_DE_statistics.npz
official_pairs.csv
```

其中 CSV 确定官方基因顺序和扰动面板，H5AD 提供各上下文的控制细胞，NPZ 提供源响应统计。目标是 300 个扰动、3 个上下文、每组 400 个细胞、18,533 个基因。数据来源及适配差异见 [DATA.md](docs/DATA.md)。

如果要复现最初的本地适配结果，应取得服务器原有的已准备统计文件。`src/prepare.py` 是上游公开数据流程，不能假设重新下载公开数据后得到的统计与本地适配版本相同。

## 推理与打包

准备好上述输入后执行：

```bash
env/bin/python scripts/audit_inputs.py
OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4 env/bin/python src/predict.py --data-dir data --output output/prediction_ABC_local.h5ad
vcc-env/bin/python src/pack.py output/prediction_ABC_local.h5ad --data-dir data --output output/prediction_ABC_local.vcc --scratch-dir scratch
```

最终 H5AD 应有 360,000 个细胞。打包器校验输入并确认打包前后预测文件的 SHA256 一致。预留充足磁盘空间；历史 H5AD 为约 4.76 GB，VCC 文件约 4.00 GB，打包另需临时存储。不要覆盖未完成的预测或已有提交包。

`scripts/run_prediction_pipeline.py` 保留历史行为，依赖已有的 `logs/local_sources_status.json`、`vcc-env/bin/python` 和 Linux `/dev/shm`。建议新目录使用上面的分步命令；历史流水线不是从空目录一键运行的入口。

## 本次验证

```bash
python tools/verify_delivery.py
```

本次检查了历史 SHA256、10 个 Python 文件的语法、注释清理前后忽略文档字符串的 AST 一致性，以及交付文件完整性。路径替换和命令行帮助保留属于显式整理修改，详见 `results/cleanup_verification.json`。未重新运行模型或重新评测，因此不将历史运行结果视为本次复现结果。

## 来源与许可

保留上游 MIT [LICENSE](LICENSE) 及 [原始 README](docs/README.upstream.md)。数据仍适用各自许可。历史结果与来源说明见 [PROVENANCE.md](docs/PROVENANCE.md)。本地交付包不含访问 token、账号凭据或自动提交流程。
