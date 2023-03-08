# 哈尔滨工业大学（深圳）计算方法实验

该仓库内容主要为哈尔滨工业大学祖传的计算方法实验，如果实验指导书有更新，可能内容将不再适用

我们实验时编程语言不限，故我选择使用 Julia 语言完成该实验

如果这个项目对你有帮助，请在右上角点一下 star ~

> 这只是个想骗星星的仓库罢了 :(

## 项目结构

> 因实验首先用 Julia 在 Jupyter Notebook 上完成，而后因为提交验收问题，不得已改为 Julia 脚本。而提交的实验报告 PDF 由 Jupyter Notebook 导出的 Markdown 转换而来，故 Notebook 所撰写的报告部分内容可能没有更新，但实验结果应当是正确的。

- `docs` 目录下包含完整的作业和实验的 PDF 提交文档，代码和源文件在 `hws` 和 `labs` 目录下查看
- `hws` 目录下包含作业的源 Markdown 文件，各班作业内容有差别，仅供参考
- `labs` 目录为实验目录，包括实验代码和报告源 Markdown 文件
  - `handout` 为实验指导书等内容
  - `julia` 为 `*.jl` 脚本文件，运行需要 `julia` 环境和相关包的依赖
  - `jupyer` 为 `*.ipynb` 笔记本，执行代码需要在笔记本安装 IJulia 核
  - `matlab` 为 `.m` 格式的 MATLAB 脚本文件，运行需要安装 MATLAB （当前实现为 ChatGPT 翻译 Julia 代码，可能无法执行）
  - `reports` 目录下为实验报告的源 Markdown 文件，包括实验报告所用的图片等，如仅需参考报告，建议查看 `docs` 目录下排版完整的 PDF 文件
    - `lab1-lagrange` 为拉格朗日插值实验报告
    - `lab2-romberg` 为龙贝格积分实验报告
    - `lab3-runge-kutta` 为龙格-库塔方法实验报告
    - `lab4-newton` 为牛顿插值实验报告
    - `lab5-gauss` 为高斯消元法实验报告

完整目录结构如下所示：

```bash
.
|-- archive
|-- docs
|   |-- homework-pdfs
|   `-- lab-reports
|-- hws
|   `-- assets
`-- labs
    |-- handout
    |-- julia
    |-- jupyter
    |-- matlab
    `-- reports
        |-- lab1-lagrange
        |-- lab2-romberg
        |-- lab3-runge-kutta
        |   `-- assets
        |-- lab4-newton
        `-- lab5-gauss
```

## 安装使用

### Jupyter Notebook

如果你有 Jupyter Notebook 环境，添加了 IJulia 核之后，可以直接打开 `labs/jupyter` 目录下的 `*.ipynb` 文件，运行代码即可。运行报错时，根据提示安装需要的第三方包。

### Julia

Julia 脚本运行需要解释器，可以在 [Julia 官网](https://julialang.org/) 下载安装包，或者使用包管理器安装。国内安装有中文官网，可以参考 [Julia 中文官网](https://cn.julialang.org/)。

### MATLAB

MATLAB 脚本的运行需要安装 MATLAB

> 注意：当前的 MATLAB 实现为 Julia 代码通过 ChatGPT 生成，尚未验证可执行性，但因两者的语法相似（数组、向量等的下标起点为 1 ，列优先保存等，如果是参考算法二者实现的差异不大）

### Release Binary

如果使用 Release 版本的二进制文件，因执行时编译运行的耗时，其执行效率较低，且只适用于 Windows 平台

## 参考资料

> 具体的参考资料引用已在各个实验报告中给出，这里不再列出
