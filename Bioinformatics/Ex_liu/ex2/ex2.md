# 第二次上机操作

## 实验目的

1. 学习并掌握使用本地工具MEGA构建系统发育树
2. 学习并掌握使用在线工具CIPRES构建系统发育树
3. 学习并掌握使用Figtree编辑展示树文件

## 实验内容
### MEGA构建系统发育树

> bioedit多序列比对
> 步骤
> + 打开序列
> + Ctrl+A选中全部序列
> + Accessory Application中选中Clustalw Multiple alignment
> + 设置参数，勾选output
> + 导出比对结果

#### 序列导入
可以通过DATA导入数据也可以通过File导入数据：Open A File/Session。

如果打开的文件是比对结果，选择Analyze；如果打开的文件是序列文件，选择Align。

#### 多序列比对
序列导入后，选中全部序列(ctrl+A),点击软件上方Alignment，选择 Align by ClustalW，比对参数保持默认。
> 如果比对的是蛋白质序列，选择 Align by MUSCLE

比对完成后保存为fasta或meg格式。Data -> Export Alignment -> 导出的格式。

#### 构建进化树
导入比对结果，在MEGA主界面的PHYLOGENY选项中选择一个算法来构建进化树。

参数选择：进化树评估：Bootstrap method，Bootstrap value 一般选择1000~1500；第一次绘图时建议选择500，这样运行速度会比较快，结果合适再调至1000重新进行进化分析。其余参数保持默认

进化树构建好后，导出为nwk格式，方便后面用figtree进行美化



### EBI、CIPRES在线构建系统发育树

#### EBI在线工具MAFFT多序列比对
打开网站[MAFFT](https://www.ebi.ac.uk/Tools/msa/mafft/) -> 选择比对的序列类型 -> 上传或粘贴fasta文件 -> 参数保持默认,输出格式为fasta -> Be notified by email，设置邮箱提醒 -> submit
比对完后，下载比对结果（Download Alignment File）。查看构建的进化树，没有设置Bootstrap。

**[CIPRES](https://www.phylo.org/)需注册**
#### 上传数据
注册完后，进入Folders。点击data -> upload data(已比对的fasta文件) -> save

#### 新建task

+ 流程 点击tasks -> create new task -> 填写task描述 -> Select Data 选择数据 -> 选择工具Select Tool 点击RAxML-HPC2 on XSEDE (8.2.12) -> 设置参数 Set Parameters -> Save Parameters -> Save and Run task

+ 具体参数的设置
  + Maximum Hours to Run (click here for help setting this correctly) ：1
  + Set a name for output files ：给输出文件命名
  + Enter the number of patterns in your dataset ：比对结束后每条序列的长度
  + Please select the Data Type ：数据类型
  + Bootstrap iterations (-#|-N) ：Bootstrap Value 设置为1000-1500。


结果文件解释
结果文件 | 含义
:-- | :--
RAxML_bestTree.ex2.tree | Contains the best-scoring ML tree of a thorough ML analysis. 
RAxML_bipartitionsBranchLabels.ex2.tree | Contains the same information as the file above, but support values are correctly displayed as Newick branch labels and not node labels! Support values always refer to branches/splits of trees and never to nodes of the tree. Note that, some tree viewers have problems displaying branch labels, they are however part of the standard Newick format. 
RAxML_bipartitions.ex2.tree |  If you used the -f b option, this file will contain the input tree with confidence values from 0 to 100 drawn on its nodes! It is also printed when -f a -x have been specified, at the end of the analysis the program will draw the BS support values on the best tree found during the ML search.
RAxML_bootstrap.ex2.tree | If a multiple bootstrap is executed by -# and -b or -x all final bootstrapped trees will be written to this one, single file.
RAxML_info.ex2.tree |  contains information about the model and algorithm used and how RAxML was called. The final GAMMA-based likelihood(s) as well as the alpha shape parameter(s) are printed to this file. In addition, if the rearrangement setting was determined automatically (-i  has not been used) the rearrangement setting found by the program will be indicated. This is the most important output file because it tells you what RAxML did and is always written irrespective of the command line option. In addition it contains information about all other output files that were written by your run.

选择RAxML_bipartitionsBranchLabels.ex2.tree。如果figtree打不开这个文件，需要用Notepad++等文本编辑器打开文件：选择正则表达式，`:([\.\d]+)\[(\d+)\]`替换为`$2:$1`，保存后就可以用Figtree打开。原因是figtree识别newick里节点标号不支持[方括号]。

### 进化树美化figtree

<!--看懂进化树
![看懂进化树](https://pic2.zhimg.com/v2-f9f3043dd508d5abfcc7b0a9bf2658ad_r.jpg)-->
+ 菜单栏
 	+ File 一般是输入、输出，新建、保存，打印等功能
	+ Edit 复制、粘贴、选择、查找等功能
	+ Tree 和工具栏功能基本一样，但是有清除功能
	+ Help 大家都懂~

+ 工具栏
	+ Layout 布局标签，红框部分是不同的呈现格式，按要求选择；Zoom 标签是调整图的大小； Expansion 标签是调整图的高度；Fish eye 进行布局调整；Root length 是调节根的长度；Cuevature 是 调整角分支的线条；
	+ Current tree：切换进化树
	+ Appearance 标签 Line weight 可调整分支的粗细；
	+ Tip lables 标签，加上样本名；还可以改变样本名大小以及颜色；
    + Tip shapes 改变tip前的节点的形状（二叉树叶的形状），有圆形、方形和菱形三种形状，可以改变大小，颜色，outline等等
	+ Node labels 是添加节点标签；可以改变样本名大小以及颜色，sig digits 选项可以改变小数保留位数；
	+ Node shapes 是改变节点形状以及大小；
	+ Branch labels 是添加分支标签，可以改变字体颜色和大小， sig digits 选项可以改变小数保留位数
	+ Scale bar 是增加标尺
	+ Scale axis 是增加最下方的标尺；勾选掉 show grid 选项可以删除竖线部分；
	+ Time scale 标签是调节最下方坐标轴的，scale by factor 是默认选项，scale root to 选项可改变轴的最大值。
	+ Hilight 标签可以改变分支的颜色；



作业详见

[课后作业](https://github.com/xujunbi/linux-bioinformatics/tree/master/bioinformatics/Ex_liu/ex2/homework)




---
参考资料
+ [FigTree构建进化树（柱状）教程](http://www.360doc.com/content/19/0125/12/52645714_811199917.shtml)
+ [一文读懂进化树（图文详解](https://zhuanlan.zhihu.com/p/141835886)
+ [BioEdit 7.0.5.3安装、多序列比对及结果导出](http://blog.sciencenet.cn/blog-3375649-1104907.html)
