# 来源与整理记录

- 下载日期：2026-09-23。
- SSH：`casia1`，迁移后地址 `211.137.21.4`，通过旧地址已信任的主机密钥验证连接。
- 原始模型目录：`/ssd2/liuqi/xiexiu/code`。
- 原始同级脚本目录：`/ssd2/liuqi/xiexiu`。
- 上游工作区提交：`ffd0c66bad1c7a57153bda2da91a842f3ba58ef2`；读取时 `git status --short` 无输出。
- 历史核心文件校验：`results/reproducibility.json`；6 个文件均匹配。
- 后续 `v1-in-barchmark-20260920`、`v2`、`v3`、`vccweb-20260919` 均未取入。

## 整理范围

`code/` 移至 `src/`；根目录的准备、审计和历史流水线脚本移至 `scripts/`。脚本以自身文件位置推导项目根目录，并适配新的相对位置。核心模型代码仅移除注释和文档字符串；CLI 的模块文档字符串改为显式帮助字符串，保留原帮助内容。

没有修改模型参数、统计公式、控制细胞选择方法或整数计数生成逻辑。此次 AST 比较是在显式路径和 CLI 文本调整之后、去注释之前建立基线；它证明注释清理没有额外改变代码结构，不等同于实际推理结果验证。

## 历史结果与当前检查

`official_latest.json` 和 `wyh_latest.json` 都记录 `score_avg=0.15592419253448622`，属于不同提交 ID 的历史结果。用户提供的分数线索没有在这些记录中出现；不能据此断言找到了 0.59 或 1.59 分的版本。

`final_validation.json` 和 `input_validation.json` 是当时服务器运行的检查记录，其绝对路径仍指向原服务器。`cleanup_verification.json` 和 `delivery_manifest.json` 才是本次交付整理生成的检查信息。

保留上游 README 和 requirements 作为历史来源说明；根 README 中的推理环境采用本地运行记录。历史记录不是新一次运行的成功证明。
