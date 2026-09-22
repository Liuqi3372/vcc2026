# 数据来源与准备

本次只下载代码、源数据清单与小型验证记录；没有下载数据或预测文件。

| 输入 | 原服务器来源 / 处理方式 |
| --- | --- |
| 官方控制细胞、基因和扰动面板 | `/ssd1/PubData/vcc2026-val-1` |
| H1 | `/ssd1/PubData/VCC2025`，使用原始准备逻辑 |
| CD4 | `raw/GWCD4i.DE_stats.h5ad`，原发布者 DE 统计及质量筛选 |
| Promoter | GENCODE v47 注释，5 kb 邻近关系 |
| K562 | `/ssd1/PubData/arc_cviplus_panel/by_cell_line/Replogle_K562_gwps/01_preprocess/adata_processed.h5ad`，`layers/counts` |
| HCT116 / HEK293T | 同一根目录下 `XAtlas_HCT116` / `XAtlas_HEK293T` 的预处理 H5AD |

本地 K562 数据已做 QC，重复基因符号求和，以 barcode 后缀恢复批次；HCT116/HEK293T 使用每个样本前 250 个保留控制细胞。所有本地源统计使用保留基因上的行和作为 CPM 分母，保留未测量基因标记。它们与上游公开原始数据流程存在差异，详见 `RUN_NOTES.original.md`。

若使用服务器已准备的数据，将 `/ssd2/liuqi/xiexiu/data` 下 README 列出的文件复制到项目 `data/`。控制 H5AD 在原服务器上是符号链接，传输时需要复制实际内容。不要将原服务器的悬空符号链接直接放入桌面交付包。

本地重新准备统计还需要 `raw/pert_counts_Training.csv` 和 `raw/gene_names.csv` 等 H1 面板。CD4/promoter 适配脚本只负责这两项，不能独立完成所有数据准备；准备 H1 时需参考原始 `prepare.py` 与 `sources.json`。

`prepare_local_atlas.py` 中的 `PUB` 保留原服务器数据路径；如果数据迁移，应先调整到相同数据集的新位置。源文件版本、顺序、QC 规则或统计文件变化都可能改变预测结果。
