# <p align='center'>panX 基本操作</p>

## 安装
```
# install by conda
mamba search panx
mamba create -n py_2.7 python=2.7
mamba install panx=1.6.0

```


## 基本使用
### 命令及参数
```
# 基本
# -fn 输入文件目录
# -sl species_name 及输出目录
# -t --threads, number of threads (default:1)
panX.py -fn data/TestSet -sl TestSet -t 32 > TestSet.log 2> TestSet.err

```



