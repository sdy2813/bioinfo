# TCR

## 生物学

T细胞受体（TCR，T cell receptor）是T细胞表面的一种蛋白质复合物，其功能是识别抗原，并在这个过程中激活T细胞。每个T细胞受体都能特异性地识别某一特定的抗原。

T细胞受体（T-cell receptor, TCR）是所有T细胞表面的特征性标志。T细胞介导的抗原识别依赖于T细胞受体与抗原-主要组织相容性复合体（major histocompatibilitycomplex, MHC）的相互作用。当TCR、抗原肽、MHC结合时，T淋巴细胞通过信号转导被激活，进入后续的免疫应答过程。TCR是由α、β两条肽链组成异二聚体，α链的可变区由VJ编码，β链的可变区由VDJ编码。由于疾病进展中，TCR会发生很大变化，因此，目前TCR研究应用于癌症、自身免疫、炎症和传染病等疾病研究。95%的T细胞的受体由α亚基和β亚基构成，另外5%的受体由γ亚基和δ亚基构成。这个比例会因为个体发育或是疾病而变化。每个链都包含一个可变区（V）和一个常量区（C）。可变区负责识别抗原，而常量区则与细胞膜和信号转导相关的其他蛋白质相结合。

![](https://upload.wikimedia.org/wikipedia/commons/4/4d/Antigen_presentation.svg)

每个T细胞都表达一种独特的T细胞受体，可以特异性地识别一种特定的抗原。这种多样性是由基因重排过程产生的，即在T细胞发育过程中，T细胞受体的基因会进行随机的切割和重组，从而产生大量的不同的T细胞受体。

在T细胞激活过程中，T细胞受体首先识别并结合到由抗原呈递细胞（APC，Antigen Presenting Cell）上的主要组织相容性复合物（MHC，Major Histocompatibility Complex）上呈递的抗原肽段。然后，T细胞受体通过一系列信号转导过程，激活T细胞，使其进行增殖和分化，产生免疫效应。

辅助型T细胞表面的CD4分子，负责识别第二类主要组织相容性复合体（MHC II）

细胞毒性T细胞表面的CD8分子，负责识别第一类主要组织相容性复合体（MHC I）

总的来说，T细胞受体是免疫系统中非常重要的分子，对于T细胞识别抗原和进行免疫应答起着关键的作用。

## T cell receptor repertoire sequencing

T细胞受体库测序（T cell receptor repertoire sequencing，TCR-Seq）是以生物信息学（Bioinformatics）全面高速分析高通量测序技术（High-throughput sequencing，HTS）检测靶向扩增后的T细胞抗原识别决定性表面分子，即T细胞受体（T cell receptor，TCR）多样性的检测技术，用以揭示机体在生理和病理状态下T细胞介导的细胞免疫应答（T cell-mediated immune response）状态改变。

T细胞介导的细胞免疫应答过程中，抗原呈递细胞（antigen-presenting cell，APC）摄取抗原（Ag）、消化形成抗原-MHC分子复合物，并呈递给T细胞。T细胞通过自身T细胞受体β链中V-D-J基因重排后的CDR3β参与抗原识别。

![](http://seqhealth.cn/upload/20190107/5c32c55021503.png)

（TCR的基因由可变区（V）、多变区（D）、结合区（J）和恒定区（C）四部分基因片段组成，形成互补决定区（complementarities determining region, CDR）和间隔的4个骨架区（framework region, FR），基因结构如下图所示。在T细胞发育过程中CDR1、2和FR区域相对保守，CDR3区由V、D和J 进行重排而形成具有功能的TCR编码基因（T细胞克隆），由于V（65~100种）、D（2种）、J（13种）基因片段本身具有多样性，此外，由于在重排的过程中，在VD及D-J的连接区经常有非模板的核苷酸的随机插入或删除，进一步增加了CDR3区的多样性。这种基因片段连接的不准确性使TCR的表达呈多样性，以识别各种不同的抗原。每个TCR链的可变结构域有三个互补决定区(CDR): CDR1、CDR2和CDR3。CDR1和CDR2，由V基因片段编码，主要通过与MHC的保守a-螺旋接触促进TCR和MHC之间的相互作用。CDR3是由V和J或D和J基因片段的之间编码，导致了高变异性。该区域负责结合MHC呈递的抗原肽(8,10)。由于其与抗原的直接相互作用和固有的高变异性，CDR3区为TCR的特异性提供了丰富的多样性，因此是TCR测序常用的靶点区域。



![](http://seqhealth.cn/upload/20181229/5c26d1152dcef.png)
<p align="center">TCR基因结构，左：α链基因结构，右：β链基因结构</p>

TCR-seq常用于评价各种免疫相关疾病和遗传性突变引起的某个物种所有T细胞或特定T细胞激活介导的细胞免疫反应中TCR基因重排碱基序列，以及各序列的丰度，用于研究不同T细胞克隆的转录情况和相互间关系，从而揭示更深层次的T细胞功能特异性，继而解释免疫应答机制、免疫耐受原因，免疫调节形式等相关生命现象。

TCR测序的应用方法包括：（1）肿瘤免疫治疗后的监测和治疗指导，（2）抗感染免疫治疗效果和耐受性评估，（3）移植排斥反应中急性排斥反应的预测和干预指导，（4）自身免疫疾病的临床试验和科学研究，（5）感染性和神经可塑性疾病的生物标志物。





---
参考资料
1. [免疫组库测序：从入门到进阶](https://ming-lian.github.io/2019/04/28/Learning-ImmuSeq/)
2. [TCR测序在肿瘤免疫治疗中的应用](http://www.immuquad.com/tcr%E6%B5%8B%E5%BA%8F%E5%9C%A8%E8%82%BF%E7%98%A4%E5%85%8D%E7%96%AB%E6%B2%BB%E7%96%97%E4%B8%AD%E7%9A%84%E5%BA%94%E7%94%A8/)









