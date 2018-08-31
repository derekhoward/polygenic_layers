from pathlib import Path
import pandas as pd

RAW_DATA = Path(
    "/Users/derek_howard/projects/HBAsets/data/raw/allen_human_fetal_brain")
DATA_DIR = Path("/Users/derek_howard/projects/polygenic_layers/data")


def read_expression_file(file_name):
    expression_df = pd.read_csv(file_name, index_col=0, header=None)
    expression_df.index.rename('probe_id', inplace=True)
    return expression_df


def read_samples_file(samples_file):
    sample_df = pd.read_csv(samples_file)
    sample_df.set_index(sample_df.index + 1, inplace=True)
    sample_df.index.rename('sample_id', inplace=True)
    return sample_df


def read_probes_data(probes_file):
    probes_df = pd.read_csv(probes_file)

    # rename columns for consistency between adult and fetal brain datasets
    if 'probeset_name' in probes_df.columns:
        probes_df.rename(columns={'probeset_name': 'probe_name',
                                  'probeset_id': 'probe_id'}, inplace=True)
    cols = ['probe_id', 'probe_name', 'gene_symbol']
    probes_df = probes_df.loc[:, cols]
    probes_df.set_index('probe_id', inplace=True)
    return probes_df


def get_donor_data(donor_file_list):
    probe_file_strings = ['Probes', 'rows_meta']
    samples_file_strings = ['Sample', 'columns_meta']
    expression_file_strings = ['Expression', 'expression']
    
    for file_path in donor_file_list:
        if any(string in file_path.stem for string in probe_file_strings):
            probes_df = read_probes_data(file_path)
        if any(string in file_path.stem for string in samples_file_strings):
            samples_df = read_samples_file(file_path)
        if any(string in file_path.stem for string in expression_file_strings):
            exp_df = read_expression_file(file_path)
        else:
            continue
    return exp_df, samples_df, probes_df


def get_exp_by_genes(probes_df, exp_df):
    annotated_probes = probes_df[['gene_symbol']].merge(exp_df,
                                                        left_index=True,
                                                        right_index=True)

    return annotated_probes.groupby('gene_symbol').mean()


def merge_sampleinfo_gene_expression(exp_by_genes, samples_df):
    df = samples_df.merge(exp_by_genes.T, left_index=True, right_index=True)
    return df.drop(['structure_id', 'well_id', 'structure_acronym'], axis=1)


def get_layer_samples(annotated_samples_df):
    layer_samples = annotated_samples_df[(annotated_samples_df.structure_name.str.contains('SG|MZ|CP|SP|IZ|SZ|VZ')) & (~annotated_samples_df.structure_name.str.contains('hippoc'))]
    return layer_samples


#if __name__ == '__main__':
