# Luolab Skills

Reusable Codex skills developed for Luolab workflows.

## Skills

### Point mutation design

Design and generate point-mutant SnapGene plasmids and mutagenesis primers for
site-directed mutagenesis, back-to-back PCR, Gibson assembly, and batch mutation
tables.

Location: [`point-mutation-design/`](point-mutation-design/)

### Transfection rate

Provide a simple, exploratory estimate of reporter transfection from
fluorescence microscopy images. The output is intended for rough relative
comparison under matched imaging conditions, not validated absolute
quantification.

Location: [`transfection-rate/`](transfection-rate/)

## Install

Copy the desired skill folder into your Codex skills directory:

```text
%USERPROFILE%\.codex\skills\
```

Restart Codex or open a new task after installation so the skill can be
discovered.

## 使用说明

本仓库收录 Luolab 工作流程中使用的 Codex 技能。安装时，将所需技能文件夹复制到
`%USERPROFILE%\.codex\skills\`，然后重启 Codex 或新建任务。

当前包含：

- `point-mutation-design`：批量设计点突变质粒、背靠背 PCR/Gibson 引物，并生成和验证 SnapGene `.dna` 文件。
- `transfection-rate`：基于 GFP 等报告基因荧光显微图像，对转染效率进行简化、粗略的探索性估算，适合相同成像条件下的组内相对比较，不应视为严谨的绝对定量结果。

## License

No license has been selected yet. All rights are reserved unless a license is
added later.
