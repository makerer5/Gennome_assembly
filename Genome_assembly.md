
# 1)  数据准备
ls
test_1.fq.gz  test_2.fq.gz
# 通过ls命令查看，发现当前目录下有两个fastq的压缩文件，它们是illumina测序获得的原始数据，也是通常要提交到NCBI SRA数据库的文件。关于如何提交请参考生信入门：如何将测序原始数据上传NCBI。
2)  数据质控
# 使用fastqc软件对原始测序reads进行质控，生成网页统计报告，快速获得数据质量的好坏。
fastqc -t 10 test_1.fq.gz test_2.fq.gz
# 3） 查看运行fastqc获得的结果文件
ls -t  test*html #-t sort by modification time, newest first
test_1_fastqc.html  test_2_fastqc.html

# 使用fastp软件对原始测序数据进行过滤，去除低质量、adapter及N碱基等，得到cleandata。
fastp -i test_1.fq.gz -I test_2.fq.gz -o clean_1.fq.gz -O clean_2.fq.gz

# 查看发现运行fastp获得结果文件clean1.fq.gz和clean2.fq.gz
ls -t
clean_2.fq.gz  fastp.html  test_1_fastqc.html  test_1_fastqc.zip  test_1.fq.gz
clean_1.fq.gz  fastp.json  test_2_fastqc.html  test_2_fastqc.zip  test_2.fq.gz

# 4）基因组组装
# 使用SPAdes (版本：3.11) 短序列组装软件对Clean Data进行组装，经多次调整参数后获得最优组装结果；然后reads将比对回组装获得的Contig上，再根据reads的paired-end和overlap关系，对组装结果进行局部组装和优化
spades.py -h
# 运行spades程序进行组装，可以自定义-k 选项获得最佳组装结果
spades.py -k 21,33,55,77,99,127 --careful  --only-assembler -1 clean_1.fq.gz  -2 clean_2.fq.gz  -o assembly -t 20

# 5）组装结果统计

# QUAST 用于基因组和宏基因组的拼接评估

    #查看quast帮助文档
    quast.py -h

# 运行quast程序
quast.py assembly/scaffolds.fasta  

# 查看拼接结果统计
more quast_results/results_2019_02_23_12_05_14/report.tsv

# 6）基因组注释
# Prokka是一个快速注释原核生物（细菌、古细菌、病毒等）基因组的软件工具。它产生GFF3、GBK和SQN文件，能够在Sequin中编辑并最终上传到Genbank/DDJB/ENA数据库。不能注释真核生物。
# 使用的注释工具：
    Prodigal 编码序列（CDS）
    RNAmmer 核糖体RNA基因（rRNA）
    Aragorn 转运RNA基因（tRNA）
    SignalP 信号肽
    Infernal 非编码RNA

# 输出：
    .fna 原始输入contigs的FASTA文件（核苷酸）
    .faa 翻译的编码基因的FASTA文件（蛋白质）
    .ffn 所有基因组特征的FASTA文件（核苷酸）
    .fsa 用于提交的Contig序列（核苷酸）
    .tbl 用于提交的特征表（Feature table）
    .sqn 用于提交的Sequin可编辑文件
    .gbk 包含序列和注释的Genbank文件
    .gff 包含序列和注释的GFF v3文件
    .log 日志文件
    .txt 注释汇总统计

# 基因组注释命令prokka
prokka --outdir annotation/ --force --locustag test --prefix test assembly/scaffolds.fasta --force  --cpus 20  --centre test --compliant

# 在线分析网站：
基因组注释 RAST
http://rast.nmpdr.org/
基因预测
http://prodigal.ornl.gov/
http://ccb.jhu.edu/software/glimmer/index.shtml
RNAmer
http://www.cbs.dtu.dk/services/RNAmmer/
tRNAscan
http://lowelab.ucsc.edu/tRNAscan-SE/
trf预测
http://tandem.bu.edu/trf/trf407b.linux.download.html
操纵子预测
http://www.microbesonline.org/operons/
http://operondb.cbcb.umd.edu/cgi-bin/operondb/operons.cgi
基因岛预测 IslandViewer
http://www.pathogenomics.sfu.ca/islandviewer/browse/
预测噬菌体
http://phaster.ca
预测信号肽
http://www.cbs.dtu.dk/services/SignalP/
毒力基因数据 VFDB
http://www.mgc.ac.cn/VFs/
耐药基因数据 CARD
https://card.mcmaster.ca/

方案二
## 1.2 FastQC质量评估

# 启动质控软件环境
    conda activate kneaddata
# (可选)使用指定位置的(别人安装的)conda
# source /home/liuyongxin/miniconda2/bin/activate
# 第一次使用软件要记录软件版本，文章方法中必须写清楚
    fastqc --version # 
# time统计运行时间，fastqc质量评估
# *.gz为原始数据，-t指定多线程
    time fastqc seq/*.gz -t 6

质控报告见`seq`目录，详细解读见[《数据的质量控制软件——FastQC》](https://mp.weixin.qq.com/s/tDMih7ISLJcL4F4sWBq3Vw)。
# multiqc将fastqc的多个报告生成单个报告查看和比较
# 记录软件版本
    multiqc --version 
# 整理fastqc报告，输出multiqc_report.html至result/qc目录
    multiqc -d seq/ -o result/qc -f 

查看右侧result/qc目录中multiqc\_report.html，单击，选择`View in Web Browser`查看可交互式报告。

## 1.3 质量控制

    mkdir -p temp/qc

### Fastp质量控制环境样品

适用于无宿主污染的环境样品，质控速度快，自动识别接头和低质量，详见：[极速的FASTQ文件质控+过滤+校正fastp](http://mp.weixin.qq.com/s/u3U-AJW7oRYTx5h13c19UQ)
# 开头记录软件版本，0.23.4，结尾为iMeta引文
# Citation: Shifu Chen. 2023. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2: e107
    fastp

# 单样本质控
    i=C1
    fastp -i seq/${i}_R1.fq.gz  -I seq/${i}_R2.fq.gz \
      -j temp/qc/${i}_fastp.json -h temp/qc/${i}_fastp.html \
      -o temp/qc/${i}_R1.fastq -O temp/qc/${i}_R2.fastq

MEGAHIT组装

# 启动工作环境
    conda activate megahit
    
    megahit -v 
    metaspades.py -v 
    quast.py -v 
    prodigal -v 

# 删除旧文件夹，否则megahit无法运行
    /bin/rm -rf temp/megahit
# 组装，10~30m，TB级数据需几天至几周
    megahit -t 6 \
        -1 temp/qc/49642_HQ_R1.fastq \
        -2 temp/qc/49642_HQ_R2.fastq \
        -o temp/megahit 
# 统计大小通常300M~5G，如果contigs太多，可以按长度筛选，降低数据量，提高基因完整度，详见附录megahit
    seqkit stat temp/megahit/final.contigs.fa
# 预览重叠群最前6行，前60列字符
    head -n6 temp/megahit/final.contigs.fa | cut -c1-60

# 备份重要结果
    mkdir -p result/megahit/
    ln -f temp/megahit/final.contigs.fa result/megahit/
# 删除临时文件
    rm -rf temp/megahit/intermediate_contigs

### 方法2. metaSPAdes精细拼接

# 精细但使用内存和时间更多，15~65m
    metaspades.py -t 6 -1 temp/qc/49642_HQ_R1.fastq  -2 temp/qc/49642_HQ_R2.fastq  -o temp/metaspades -m 200
	
# 查看软件时间User time和内存Maximum resident set size
    cat metaspades.py.log
# 2.3M，contigs体积更大
    ls -sh temp/metaspades/contigs.fasta
    seqkit stat temp/metaspades/contigs.fasta

# 备份重要结果
    mkdir -p result/metaspades/
    ln -f temp/metaspades/contigs.fasta result/metaspades/
# 删除临时文件
    /bin/rm -rf temp/metaspades
注：metaSPAdes支持二、三代混合组装，见附录，此外还有OPERA-MS组装二、三代方案

### QUAST评估
 
# QUAST评估，生成report文本tsv/txt、网页html、PDF等格式报告
    quast.py result/megahit/final.contigs.fa \
      -o result/megahit/quast -t 6

# (可选) megahit和metaspades比较
    quast.py --label "megahit,metapasdes" \
        result/megahit/final.contigs.fa \
        result/metaspades/contigs.fasta \
        -o result/quast

 # (可选)metaquast评估，更全面，但需下载相关数据库，受网速影响可能时间很长(经常失败)
# metaquast based on silva, and top 50 species genome, 18min
    time metaquast.py result/megahit/final.contigs.fa \
      -o result/megahit/metaquast

## 3.2 基因预测、去冗余和定量Gene prediction, cluster & quantitfy

### metaProdigal基因预测Gene prediction

# 输入文件：组装的序列 result/megahit/final.contigs.fa
# 输出文件：prodigal预测的基因序列 temp/prodigal/gene.fa
# 基因大，可参考附录prodigal拆分基因文件，并行计算    
    mkdir -p temp/prodigal
# prodigal的meta模式预测基因，>和2>&1记录分析过程至gene.log
    prodigal -i result/metaspades/contigs.fasta \
        -d temp/prodigal/gene.fa \
        -o temp/prodigal/gene.gff \
        -p meta -f gff > temp/prodigal/gene.log 2>&1 
# 查看日志是否运行完成，有无错误
    tail temp/prodigal/gene.log
# 统计基因数量
    seqkit stat temp/prodigal/gene.fa 
# 统计完整基因数量，数据量大可只用完整基因部分
    grep -c 'partial=00' temp/prodigal/gene.fa 
# 提取完整基因(完整片段获得的基因全为完整，如成环的细菌基因组)
    grep 'partial=00' temp/prodigal/gene.fa | cut -f1 -d ' '| sed 's/>//' > temp/prodigal/full_length.id
    seqkit grep -f temp/prodigal/full_length.id temp/prodigal/gene.fa > temp/prodigal/full_length.fa
    seqkit stat temp/prodigal/full_length.fa

### cd-hit基因聚类/去冗余cluster & redundancy

# 输入文件：prodigal预测的基因序列 temp/prodigal/gene.fa
# 输出文件：去冗余后的基因和蛋白序列：result/NR/nucleotide.fa, result/NR/protein.fa

    mkdir -p result/NR
# aS覆盖度，c相似度，G局部比对，g最优解，T多线程，M内存0不限制
# 2万基因2m，2千万需要2000h，多线程可加速
    cd-hit-est -i temp/prodigal/gene.fa \
        -o result/NR/nucleotide.fa \
        -aS 0.9 -c 0.95 -G 0 -g 0 -T 0 -M 0
# 统计非冗余基因数量，单次拼接结果数量下降不大，多批拼接冗余度高
    grep -c '>' result/NR/nucleotide.fa
# 翻译核酸为对应蛋白序列, --trim去除结尾的*
    seqkit translate --trim result/NR/nucleotide.fa \
        > result/NR/protein.fa 

### salmon基因定量quantitfy
# 输入文件：去冗余后的基因序列：result/NR/nucleotide.fa
# 输出文件：Salmon定量：result/salmon/gene.count, gene.TPM

    mkdir -p temp/salmon
    salmon -v 

# 建索引, -t序列, -i 索引，10s
    salmon index -t result/NR/nucleotide.fa \
      -p 3 -i temp/salmon/index 

# 定量，l文库类型自动选择，p线程，--meta宏基因组模式, 2个任务并行2个样
    salmon quant -i temp/salmon/index -l A -p 3 --meta \
        -1 temp/qc/49642_HQ_R1.fastq -2 temp/qc/49642_HQ_R2.fastq \
        -o temp/salmon/49642_HQ.quant

# 合并
    mkdir -p result/salmon
    salmon quantmerge --quants temp/salmon/*.quant \
        -o result/salmon/gene.TPM
    salmon quantmerge --quants temp/salmon/*.quant \
        --column NumReads -o result/salmon/gene.count
    sed -i '1 s/.quant//g' result/salmon/gene.*

# 预览结果表格
    head -n3 result/salmon/gene.*

## 3.3 功能基因注释Functional gene annotation

# 输入数据：上一步预测的蛋白序列 result/NR/protein.fa
# 中间结果：temp/eggnog/protein.emapper.seed_orthologs
#           temp/eggnog/output.emapper.annotations
#           temp/eggnog/output

# COG定量表：result/eggnog/cogtab.count
#            result/eggnog/cogtab.count.spf (用于STAMP)

# KO定量表：result/eggnog/kotab.count
#           result/eggnog/kotab.count.spf  (用于STAMP)

# CAZy碳水化合物注释和定量：result/dbcan2/cazytab.count
#                           result/dbcan2/cazytab.count.spf (用于STAMP)

# 抗生素抗性：result/resfam/resfam.count
#             result/resfam/resfam.count.spf (用于STAMP)

# 这部分可以拓展到其它数据库

### eggNOG基因注释gene annotation(COG/KEGG/CAZy)
软件主页：https://github.com/eggnogdb/eggnog-mapper

# 运行并记录软件版本
    conda activate eggnog
    emapper.py --version 
# emapper-2.1.10 / Expected eggNOG DB version: 5.0.2 / 
# Installed eggNOG DB version: 5.0.2 / (eggnog6在线可以用，本地用不了，还没释放)
# Diamond version found: diamond version 2.0.15 / 
# MMseqs2 version found: 13.45111

# 运行emapper，18m，默认diamond 1e-3
    mkdir -p temp/eggnog
    time emapper.py --data_dir ${db}/eggnog \
      -i result/NR/protein.fa --cpu 6 -m diamond --override \
      -o temp/eggnog/output

# 格式化结果并显示表头
    grep -v '^##' temp/eggnog/output.emapper.annotations | sed '1 s/^#//' \
      > temp/eggnog/output
    csvtk -t headers -v temp/eggnog/output

 # 生成COG/KO/CAZy丰度汇总表
    mkdir -p result/eggnog
 # 显示帮助
    summarizeAbundance.py -h
# 汇总，7列COG_category按字母分隔，12列KEGG_ko和19列CAZy按逗号分隔，原始值累加
    summarizeAbundance.py \
      -i result/salmon/gene.TPM \
      -m temp/eggnog/output \
      -c '7,12,19' -s '*+,+,' -n raw \
      -o result/eggnog/eggnog
    sed -i 's#^ko:##' result/eggnog/eggnog.KEGG_ko.raw.txt
    sed -i '/^-/d' result/eggnog/eggnog*
 # eggnog.CAZy.raw.txt  eggnog.COG_category.raw.txt  eggnog.KEGG_ko.raw.txt

# 添加注释生成STAMP的spf格式
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
      ${db}/EasyMicrobiome/kegg/KO_description.txt \
      result/eggnog/eggnog.KEGG_ko.raw.txt | \
      sed 's/^\t/Unannotated\t/' \
      > result/eggnog/eggnog.KEGG_ko.TPM.spf
    head -n 5 result/eggnog/eggnog.KEGG_ko.TPM.spf
# KO to level 1/2/3
    summarizeAbundance.py \
      -i result/eggnog/eggnog.KEGG_ko.raw.txt \
      -m ${db}/EasyMicrobiome/kegg/KO1-4.txt \
      -c 2,3,4 -s ',+,+,' -n raw \
      -o result/eggnog/KEGG
     
# CAZy
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
       ${db}/EasyMicrobiome/dbcan2/CAZy_description.txt result/eggnog/eggnog.CAZy.raw.txt | \
      sed 's/^\t/Unannotated\t/' > result/eggnog/eggnog.CAZy.TPM.spf

# COG
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1],$0}' \
      /${db}/EasyMicrobiome/eggnog/COG.anno result/eggnog/eggnog.COG_category.raw.txt > \
      result/eggnog/eggnog.COG_category.TPM.spf

### CAZy碳水化合物酶

    conda activate eggnog
 # 比对CAZy数据库, 用时2~18m
    mkdir -p temp/dbcan2
# --sensitive慢10倍，dbCAN2推荐e值为1e-102，此处结果3条太少，以1e-3为例演示
    diamond blastp \
      --db ${db}/dbcan2/CAZyDB.08062022 \
      --query result/NR/protein.fa \
      --threads 6 -e 1e-3 --outfmt 6 --max-target-seqs 1 --quiet \
      --out temp/dbcan2/gene_diamond.f6
    wc -l temp/dbcan2/gene_diamond.f6
# 整理比对数据为表格 
    mkdir -p result/dbcan2
# 提取基因与dbcan分类对应表
    format_dbcan2list.pl \
      -i temp/dbcan2/gene_diamond.f6 \
      -o temp/dbcan2/gene.list 
# 按对应表累计丰度，依赖
    summarizeAbundance.py \
      -i result/salmon/gene.TPM \
      -m temp/dbcan2/gene.list \
      -c 2 -s ',' -n raw \
      -o result/dbcan2/TPM
 # 添加注释生成STAMP的spf格式，结合metadata.txt进行差异比较
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
       ${db}/EasyMicrobiome/dbcan2/CAZy_description.txt result/dbcan2/TPM.CAZy.raw.txt | \
      sed 's/^\t/Unannotated\t/' \
      > result/dbcan2/TPM.CAZy.raw.spf
      
    head result/dbcan2/TPM.CAZy.raw.spf
 # 检查未注释数量，有则需要检查原因
    grep 'Unannotated' result/dbcan2/TPM.CAZy.raw.spf|wc -l

### CARD耐药基因

CARD在线分析平台：https://card.mcmaster.ca/ 
本地软件使用教程: https://github.com/arpcard/rgi
参考文献：http://doi.org/10.1093/nar/gkz935

    mkdir -p result/card
# 启动rgi环境和记录版本
    conda activate rgi
    rgi main -v # 5.2.1
    
# 简化蛋白ID
    cut -f 1 -d ' ' result/NR/protein.fa > temp/protein.fa
 # 这个错误忽略即可，不是报错，没有任何影响  grep: 写错误: 断开的管道
    grep '>' result/NR/protein.fa | head -n 3
    grep '>' temp/protein.fa | head -n 3
 # 蛋白层面注释ARG
    rgi main -i temp/protein.fa -t protein \
      -n 6 -a DIAMOND --include_loose --clean \
      -o result/card/protein
    head -n3 result/card/protein.txt
    
# 基因层面注释ARG
    cut -f 1 -d ' ' result/NR/nucleotide.fa > temp/nucleotide.fa
    grep '>' temp/nucleotide.fa | head -n3
    rgi main -i temp/nucleotide.fa -t contig \
      -n 6 -a DIAMOND --include_loose --clean \
      -o result/card/nucleotide
    head -n3 result/card/nucleotide.txt
    
# 重叠群层面注释ARG
    cut -f 1 -d ' ' result/metaspades/contigs.fasta > temp/contigs.fa
    grep '>' temp/contigs.fa | head -n3
    rgi main -i temp/contigs.fa -t contig \
      -n 6 -a DIAMOND --include_loose --clean \
      -o result/card/contigs
    head result/card/contigs.txt

结果说明：
- protein.json，在线可视化
- protein.txt，注释基因列表
