# extractUpStreamSequence
用来提取给定基因列表的基因上游区域用法如下
可以根据你的文件中要提取的feature的行来调整--feature
`python extractUpStreamSequence_n.py --genome genome.fasta --gene lc_genelist.txt --gff merged_sigua.gtf --feature "transcript"`
