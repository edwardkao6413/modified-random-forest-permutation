import pandas as pd
import random
from collections import Counter
class helps:
    def __init__(self):
        pass

    def __str__(self):
        return 'You need other packages to run it. Using .packages to see the package requirements.'


class packages:
    def __init__(self):
        pass

    def __str__(self):
        return 'pandas, random, collections\' Counter.'


class MRF:
    def __init__(self, df, file_safe_pos_name):
        self.df = df
        self.file_safe_pos_name = file_safe_pos_name

    def run(self, bootstrap_rounds, permuted_rounds, ranks):
        bootstrap_genes_record = []
        ct_genes_score = self.df
        ct_genes_score.columns = ['gene', 'score']
        new_gene_score = ct_genes_score.copy()
        new_gene_score.columns = ['gene', 'permuted score']

        all_genes = list(ct_genes_score['gene'])
        original_time = [0 for i in range(len(ct_genes_score['gene']))]
        boots_over_permute_times = dict(zip(all_genes, original_time))

        for bootstrap_times in range(bootstrap_rounds):
            bootstrap_genes = ct_genes_score.sample(n=len(ct_genes_score.index), replace=True)
            bootstrap_score = bootstrap_genes.sort_values(['score'], ascending=False)[0:ranks]

            for bootstrap_gene in bootstrap_score['gene']:
                bootstrap_genes_record.append(bootstrap_gene)

            for permutation_times in range(permuted_rounds):
                permuted_genes = random.sample(list(ct_genes_score['gene']), len(ct_genes_score['gene']))
                new_gene_score['gene'] = permuted_genes

                score = []
                permu_score_dict = dict(zip(list(new_gene_score['gene']), list(new_gene_score['permuted score'])))
                for i in range(len(bootstrap_score['gene'])):
                    score.append(permu_score_dict[list(bootstrap_score['gene'])[i]])

                final_genes_2score = bootstrap_score
                final_genes_2score['permuted'] = score

                for i in range(len(final_genes_2score['gene'])):
                    if final_genes_2score.iloc[i, :]['score'] > final_genes_2score.iloc[i, :]['permuted']:
                        boots_over_permute_times[final_genes_2score.iloc[i, :]['gene']] = boots_over_permute_times[
                                                                                              final_genes_2score.iloc[i,
                                                                                              :]['gene']] + 1
            print('finished round {} bootstrap'.format(bootstrap_times + 1))

        all_result = pd.DataFrame()
        all_result['gene'] = list(boots_over_permute_times.keys())
        all_result['boots over times'] = list(boots_over_permute_times.values())

        gene_frequency = sorted(Counter(bootstrap_genes_record).items(), key=lambda x: x[1], reverse=True)
        gene2 = []
        frequency = []
        for i in range(len(gene_frequency)):
            gene2.append(gene_frequency[i][0])
            frequency.append(gene_frequency[i][1])
        file2 = pd.DataFrame()
        file2['gene'] = gene2
        file2['frequency'] = frequency
        all_result.to_csv(self.file_safe_pos_name + 'boots.csv')
        file2.to_csv(self.file_safe_pos_name + 'freq.csv')

    def run_nofreq(self, bootstrap_rounds, permuted_rounds, ranks):
        ct_genes_score = self.df.iloc[:, [0, -3]]
        ct_genes_score.columns = ['gene', 'score']
        new_gene_score = ct_genes_score.copy()
        new_gene_score.columns = ['gene', 'permuted score']

        all_genes = list(ct_genes_score['gene'])
        original_time = [0 for i in range(len(ct_genes_score['gene']))]
        boots_over_permute_times = dict(zip(all_genes, original_time))

        for bootstrap_times in range(bootstrap_rounds):
            bootstrap_genes = ct_genes_score.sample(n=len(ct_genes_score.index), replace=True)
            bootstrap_score = bootstrap_genes.sort_values(['score'], ascending=False)[0:ranks]

            for permutation_times in range(permuted_rounds):
                permuted_genes = random.sample(list(ct_genes_score['gene']), len(ct_genes_score['gene']))
                new_gene_score['gene'] = permuted_genes

                score = []
                permu_score_dict = dict(zip(list(new_gene_score['gene']), list(new_gene_score['permuted score'])))
                for i in range(len(bootstrap_score['gene'])):
                    score.append(permu_score_dict[list(bootstrap_score['gene'])[i]])

                final_genes_2score = bootstrap_score
                final_genes_2score['permuted'] = score

                for i in range(len(final_genes_2score['gene'])):
                    if final_genes_2score.iloc[i, :]['score'] > final_genes_2score.iloc[i, :]['permuted']:
                        boots_over_permute_times[final_genes_2score.iloc[i, :]['gene']] = boots_over_permute_times[
                                                                                              final_genes_2score.iloc[i,
                                                                                              :]['gene']] + 1
            print('finished round {} bootstrap'.format(bootstrap_times + 1))

        all_result = pd.DataFrame()
        all_result['gene'] = list(boots_over_permute_times.keys())
        all_result['boots over times'] = list(boots_over_permute_times.values())
        all_result.to_csv(self.file_safe_pos_name + 'boots.csv')

