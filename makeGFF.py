# file: makeGFF.py
import pandas


class GFFmaker:
    def __init__(self, file, sep=','):
        if sep == ',':
            self.init_list = pandas.read_csv(file)
        elif sep == '\t':
            self.init_list = pandas.read_table(file)

    def make_list(self):
        '''
        提取输入文件内容
        '''
        must_need_header = ['type', 'start', 'end', 'geneID']
        self.gff_list = pandas.DataFrame()
        for item in must_need_header:
            if item not in self.init_list.columns:
                print('没有找到' + item + '。它是必需的。')
            else:
                self.gff_list[item] = self.init_list[item]
        if 'seqid' in self.init_list.columns:
            self.gff_list['seqid'] = self.init_list['seqid']
        else:
            self.gff_list['seqid'] = 'unknown'
        if 'source' in self.init_list.columns:
            self.gff_list['source'] = self.init_list['source']
        else:
            self.gff_list['source'] = 'GFFmaker'
        if 'score' in self.init_list:
            self.gff_list['score'] = self.init_list['score']
        else:
            self.gff_list['score'] = '.'
        if 'strand' in self.init_list:
            self.gff_list['strand'] = self.init_list['strand']
        else:
            self.gff_list['strand'] = '?'
        if 'phase' in self.init_list:
            self.gff_list['phase'] = self.init_list['phase']
        else:
            self.gff_list['phase'] = '0'
        if 'attribute' in self.init_list:
            self.gff_list['more'] = self.init_list['attribute']
        del self.init_list

    def make_gff_content(self):
        '''
        构造 gff 矩阵
        '''
        self.make_list()
        # tem_table 保存 self.gff_list 并列排序
        tem_table = pandas.DataFrame()
        must_need_header = [
            'seqid', 'source', 'type', 'start', 'end', 'score', 'strand',
            'phase'
        ]
        for item in must_need_header:
            tem_table[item] = self.gff_list[item]
        tem_table['geneID'] = self.gff_list['geneID']
        tem_table['attribute'] = None
        del self.gff_list

        def sort_feature(gene_feature_table: pandas.DataFrame):
            '''
            排序一个基因内部的所有 feature
            '''
            gene_feature_table_modified = pandas.DataFrame()
            for feature in gene_feature_table.iterrows():
                tem_feature = feature[1]
                tem_feature['attribute'] = 'ID=' + tem_feature['type'] + ':' + tem_feature['geneID'] + ';Parent=gene:' + tem_feature['geneID']
                gene_feature_table_modified = gene_feature_table_modified.append(tem_feature)
            gene_feature_table_modified.sort_values(by=['end', 'start'], ascending=[False, True])
            return gene_feature_table_modified

        # gene_table contain all genes
        gene_table = tem_table[tem_table['type'] == 'gene']
        gene_table = gene_table.sort_values(by='start')
        self.gff_table = pandas.DataFrame()
        for gene in gene_table.iterrows():
            tem_feature = tem_table[tem_table['geneID'] == gene[1]['geneID']]
            tem_gene = gene[1]
            tem_gene['attribute'] = 'ID=' + tem_gene['type'] + ':' + tem_gene['geneID'] + ';' + 'Name:' + tem_gene['geneID']
            tem_feature = sort_feature(tem_feature)
            self.gff_table = self.gff_table.append(tem_feature)

    def export_gff(self):
        self.make_gff_content()
        self.gff_table[['start', 'end']] = self.gff_table[['start', 'end']].astype('int')
        gff_last = pandas.DataFrame()
        header = [
            'seqid', 'source', 'type', 'start', 'end', 'score', 'strand',
            'phase', 'attribute'
        ]
        for item in header:
            gff_last[item] = self.gff_table[item]
        gff_last.to_csv('out.gff', sep='\t', header=False, index=False)
