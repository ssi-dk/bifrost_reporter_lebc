#!/usr/bin/env python

# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import numpy as np

def concatenate_vector(x, sep=","):
    return ",".join([str(i) for i in x])



def parse_csvs(folder, analysis_type, rename_dict = None, sort_by = ''):
    pd_dfs = [pd.read_csv(i, index_col=0) for i in os.listdir(folder) if re.match(f'.*{analysis_type}.*csv', i)]
    ssi_df = pd_dfs[0].loc[pd_dfs[0].index.str.contains('SSI')]
    pd_dfs = [i.loc[~i.index.str.contains('SSI')] for i in pd_dfs] # remove ssi duplicates form each of them
    pd_dfs.append(ssi_df) # and readd the ssi df by itself
    altered_dfs = []
    for df in pd_dfs:
        list_of_joined_dfs = []
        if any(df.index.duplicated()): # case with multiple row values per sample, i.e. resistance genes
            unique_indexes = pd.unique(df.index)
            joined_vectors = []
            for idx in unique_indexes:
                pd_df_subset = df.loc[[idx]]
                #print(pd_df_subset)
                subset_as_frame = pd_df_subset.apply(concatenate_vector, axis=0).to_frame().T
                subset_as_frame.rename(index={0:idx}, inplace=True)
                joined_vectors.append(subset_as_frame)
                #print(subset_as_frame)
            joined_df = pd.concat(joined_vectors)
            list_of_joined_dfs.append(joined_df)
        else:
            list_of_joined_dfs.append(df)
        altered_dfs.append(pd.concat(list_of_joined_dfs))
    concatenated_df = pd.concat(altered_dfs)
    if rename_dict:
        concatenated_df.rename(columns= rename_dict, inplace=True)
    pattern = r"(.*_WGS)_(EQA_[0-9]+.*)"
    kma_pattern_matches = [re.match(pattern, i).groups() for i in concatenated_df.index]
    kmas = sorted(list(set(i[0] for i in kma_pattern_matches)))
    samples = sorted(list(set(i[1] for i in kma_pattern_matches)))
    # sort data frame
    if sort_by:
        concatenated_df = concatenated_df.sort_values(by=sort_by).sort_index()
    return [concatenated_df, kmas, samples]

def parse_dfs(df_list, analysis_type, rename_dict = None, sort_by = ''):
    pd_dfs = df_list # I copy pastaed this function from the other one because I was working with csvs locally originally
    ssi_df = pd_dfs[0].loc[pd_dfs[0].index.str.contains('SSI')]
    pd_dfs = [i.loc[~i.index.str.contains('SSI')] for i in pd_dfs] # remove ssi duplicates form each of them
    pd_dfs.append(ssi_df) # and readd the ssi df by itself
    altered_dfs = []
    for df in pd_dfs:
        list_of_joined_dfs = []
        if any(df.index.duplicated()): # case with multiple row values per sample, i.e. resistance genes
            unique_indexes = pd.unique(df.index)
            joined_vectors = []
            for idx in unique_indexes:
                pd_df_subset = df.loc[[idx]]
                subset_as_frame = pd_df_subset.apply(concatenate_vector, axis=0).to_frame().T
                subset_as_frame.rename(index={0:idx}, inplace=True)
                joined_vectors.append(subset_as_frame)
            joined_df = pd.concat(joined_vectors)
            list_of_joined_dfs.append(joined_df)
        else:
            list_of_joined_dfs.append(df)
        altered_dfs.append(pd.concat(list_of_joined_dfs))
    concatenated_df = pd.concat(altered_dfs)
    if rename_dict:
        concatenated_df.rename(columns= rename_dict, inplace=True)
    pattern = r"(.*_WGS)_(EQA_[0-9]+.*)"
    kma_pattern_matches = [re.match(pattern, i).groups() for i in concatenated_df.index]
    kmas = sorted(list(set(i[0] for i in kma_pattern_matches)))
    samples = sorted(list(set(i[1] for i in kma_pattern_matches)))
    # sort data frame
    if sort_by:
        concatenated_df = concatenated_df.sort_values(by=sort_by).sort_index()
    return [concatenated_df, kmas, samples]


def anonymize_kma_df(input_df, kmas, kma_name):
    kma_remapped_names = {}
    kma_count = list(range(1, len(kmas)))
    kma_count.reverse() # this metho is in place
    for i in kmas:
        if i != kma_name and i != 'SSI_BTP_WGS':
            kma_remapped_names[i] = 'other_kma_' + str(kma_count.pop())
        else:
            kma_remapped_names[i] = i
    anonymized_df = input_df.copy()
    for i in kma_remapped_names.keys():
        anonymized_df.index = anonymized_df.index.str.replace(i, kma_remapped_names[i])
    return [anonymized_df, kma_remapped_names]

def create_kma_sample_df(anonymized_kma_df, kma_name, sample_name, kma_name_map):
    kma_sample = '_'.join([kma_name, sample_name])
    #kma_match_regex = f"({kma_name}_EQA_[0-9])+.*)"
    #kma_sample_match = re.match(kma_match_regex, kma_sample)
    #kma_sample_base = kma_sample_match.groups[0]
    other_kma_matches = list(anonymized_kma_df.index[anonymized_kma_df.index.str.contains(f"{kma_name}.*{sample_name[0:6]}")]) # match duplicates for the same kma
    other_kma_matches = [i for i in other_kma_matches if i!= kma_sample] # remove original element
    ssi_match = '_'.join(['SSI_BTP_WGS', sample_name[0:6]])
    anon_kmas = [kma_name_map[i] for i in kma_name_map.keys() if re.match('other', kma_name_map[i])]
    other_samples = []
    for i in anon_kmas:
        matches = list(anonymized_kma_df.index[anonymized_kma_df.index.str.contains(f"{i}.*{sample_name[0:6]}")]) # catch duplicates
        kma_sample_base = '_'.join([i, sample_name[0:6]])
        if kma_sample_base not in matches:
            matches.append(kma_sample_base)
        other_samples += matches
    kma_samples = [kma_sample] + other_kma_matches + [ssi_match] + other_samples
    for i in kma_samples:
        if i not in anonymized_kma_df.index: # we want to catch missing data as well
            #if kma_name not in KMA_WITH_DUPLICATES and re.match('.*-[2-9]', i): # enumerate the sample but skip if it's not the ones we know have duplicates
                #continue
            #else:
            new_row = pd.Series(['no_data'] * anonymized_kma_df.shape[1], index=anonymized_kma_df.columns, name=i)
            anonymized_kma_df = pd.concat([anonymized_kma_df, new_row.to_frame().T])
    sub_df = anonymized_kma_df.loc[kma_samples]
    return sub_df

def generate_color_matrix(input_df) -> np.array:
    """
        by the design the input to this function has
        row 0 = KMA being compared
    """
    df_array = input_df.to_numpy()
    reference_row = df_array[0, :]
    #ssi_row = df_array[1, :]
    #other_kma_rows = df_array[2:, :]
    other_rows = df_array[1:, :]
    color_matrix = np.full(df_array.shape, 'w')
    color_matrix [0, :] = 'cyan' # row being compared
    for i in range(other_rows.shape[0]):
        bool_array = reference_row == other_rows[i, :]
        #print(bool_array.shape)
        string_array = np.where(bool_array, 'g', 'r')
        color_matrix[i+1, :] = string_array
    return color_matrix

def append_metrics_to_df_and_cmat(input_df:pd.DataFrame, cmat:np.array, comparison_col:str = 'all', column_subset:list = []) -> list:
    reference_row = input_df.iloc[0]
    if column_subset:
        col_idx = [input_df.columns.get_loc(i) for i in column_subset]
        input_df = input_df.iloc[:, col_idx].copy()
        cmat = cmat[:, np.array(col_idx)].copy()
    if comparison_col == 'all':
        matches = input_df.apply(lambda x: str(sum(x==reference_row)) + f'/{len(input_df.columns)}', axis = 1)
        #pct_match = input_df.apply(lambda x: str(sum(x==reference_row)/len(input_df.columns)), axis = 1)
        input_df['n_match'] = matches
        #input_df['pct_match'] = pct_match
    else:
        reference_vector = set(reference_row[comparison_col].split(','))
        jcd = input_df[comparison_col].apply(lambda x: str(len(reference_vector.intersection(set(x.split(',')))) / len(reference_vector.union(set(x.split(','))))))
        input_df['jaccard'] = jcd
    #print(input_df.shape, cmat_extension.shape
    cmat_extension = color_matrix = np.full((len(input_df.index), input_df.shape[1] - cmat.shape[1]), 'w')
    #print(cmat_extension)
    cmat = np.concatenate([cmat, cmat_extension], axis = 1)
    return [input_df, cmat]

def generate_table(df, cmat, path, figsize = (10,5), max_width=70, fontsize=10, row_height=0.2, shorten=False):
    columns = df.columns
    rows = df.index

    if shorten:
        cell_text = df.map(lambda x: textwrap.shorten(str(x), width=max_width)).to_numpy()
    else:
        cell_text = df.to_numpy()
    # Add a table at the bottom of the axes
    colors = cmat
    
    fig, ax = plt.subplots()
    fig.set_size_inches(*figsize)
    #plt.figure(figsize=(30, 19))
    ax.axis('tight')
    ax.axis('off')
    table = ax.table(cellText=cell_text,cellColours=colors,
                         colLabels=columns,loc='center', rowLabels=rows, )

    
    #fig.set_size_inches(20,10)
    #plt.show()
    plt.savefig(path, dpi=300)
    plt.close()

def generate_concatenated_dfs_for_a_kma(kma_list:list, kma_name:str, input_df:pd.DataFrame, sample_names, concatenate = True):
    anonymized_df, kma_name_map = anonymize_kma_df(input_df, kmas, kma_name)
    kma_sample_dfs = []
    color_matrices = []
    for sample in sample_names:
        #print(sample)
        if '_'.join([kma_name, sample]) not in anonymized_df.index:
            continue # if sample name is a duplicate but not duplicated for this kma, skip
        kma_sample_df = create_kma_sample_df(anonymized_df, kma_name, sample, kma_name_map)
        color_matrix = generate_color_matrix(kma_sample_df)
        kma_sample_df_extended, color_matrix_extended = append_metrics_to_df_and_cmat(kma_sample_df, color_matrix)
        kma_sample_dfs.append(kma_sample_df_extended)
    concatenated_df = pd.concat(kma_sample_dfs, axis = 0)
    concatenated_colors = np.concatenate(color_matrices, axis = 0)
    if concatenate:
        return dict(df = pd.concat(kma_sample_dfs, axis = 0), cmat = np.concatenate(color_matrices, axis = 0))
    else:
        return dict(dfs = kma_sample_dfs, cmats = color_matrices)

def spam_tables(parsed_df:pd.DataFrame, samples:list, kmas:list, kma:str, comparison_col:str, column_subset:list, output_folder:str, analysis_type:str, file_ending='.png'):
    anonymized_df, kma_name_map = anonymize_kma_df(parsed_df, kmas, kma)
    for sample in samples:
        if kma not in ['RH_BTP_WGS'] and bool(re.match('.*-[2-9]', sample)): # so far RH is the only one with duplicates
            #print(kma, sample)
            continue # we don't want to loop over kma_samples that don't have duplicates
        kma_sample_df = create_kma_sample_df(anonymized_df, kma, sample, kma_name_map)
        cmat = generate_color_matrix(kma_sample_df)
        extended_df, extended_cmat = append_metrics_to_df_and_cmat(kma_sample_df, cmat, comparison_col, column_subset)
        output_file_name = os.path.join(output_folder, '_'.join([kma, analysis_type, sample]) + file_ending) 
        generate_table(extended_df, extended_cmat, output_file_name, figsize = (15, 5))

def spam_tables_for_all_kmas(parsed_df:pd.DataFrame, samples:list, kmas:list, comparison_col:str, column_subset:list, output_folder:str, analysis_type:str, file_ending='.png'):
    for kma in kmas:
        spam_tables(parsed_df, samples, kmas, kma, comparison_col, column_subset, output_folder, analysis_type, file_ending)
