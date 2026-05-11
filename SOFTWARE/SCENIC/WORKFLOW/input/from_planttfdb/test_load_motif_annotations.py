# input
fname = "/data/work/scenicplus/input/from_planttfdb/process_motif/Os_TF_binding_motifs_information.tbl"
column_names=('#motif_id', 'gene_name',
             'motif_similarity_qvalue', 'orthologous_identity', 'description')
motif_similarity_fdr = 0.001
orthologous_identity_threshold = 0.0

# test load_motif_annotations

import pandas as pd

if fname is None:
    if specie == 'mus_musculus':
        name='mgi'
    elif specie == 'homo_sapiens':
        name='hgnc'
    elif specie == 'drosophila_melanogaster':
        name='flybase'
    fname = 'https://resources.aertslab.org/cistarget/motif2tf/motifs-'+version+'-nr.'+name+'-m0.001-o0.0.tbl'
df = pd.read_csv(fname, sep='\t', usecols=column_names)
df.rename(columns={'#motif_id':"MotifID",
                   'gene_name':"TF",
                   'motif_similarity_qvalue': "MotifSimilarityQvalue",
                   'orthologous_identity': "OrthologousIdentity",
                   'description': "Annotation" }, inplace=True)
df = df[(df["MotifSimilarityQvalue"] <= motif_similarity_fdr) &
        (df["OrthologousIdentity"] >= orthologous_identity_threshold)]

# Direct annotation
df_direct_annot = df[df['Annotation'] == 'gene is directly annotated']
try:
    df_direct_annot = df_direct_annot.groupby(['MotifID'])['TF'].apply(lambda x: ', '.join(list(set(x)))).reset_index()
except:
    pass
df_direct_annot.index = df_direct_annot['MotifID']
df_direct_annot = pd.DataFrame(df_direct_annot['TF'])
df_direct_annot.columns = ['Direct_annot']
# Indirect annotation - by motif similarity
motif_similarity_annot = df[df['Annotation'].str.contains('similar') & ~df['Annotation'].str.contains('orthologous')]
try:
    motif_similarity_annot = motif_similarity_annot.groupby(['MotifID'])['TF'].apply(lambda x: ', '.join(list(set(x)))).reset_index()
except:
    pass
motif_similarity_annot.index =  motif_similarity_annot['MotifID']
motif_similarity_annot = pd.DataFrame(motif_similarity_annot['TF'])
motif_similarity_annot.columns = ['Motif_similarity_annot']
# Indirect annotation - by orthology
orthology_annot = df[~df['Annotation'].str.contains('similar') & df['Annotation'].str.contains('orthologous')]
try:
    orthology_annot = orthology_annot.groupby(['MotifID'])['TF'].apply(lambda x: ', '.join(list(set(x)))).reset_index()
except:
    pass
orthology_annot.index = orthology_annot['MotifID']
orthology_annot = pd.DataFrame(orthology_annot['TF'])
orthology_annot.columns = ['Orthology_annot']
# Indirect annotation - by orthology
motif_similarity_and_orthology_annot = df[df['Annotation'].str.contains('similar') & df['Annotation'].str.contains('orthologous')]
try:
    motif_similarity_and_orthology_annot = motif_similarity_and_orthology_annot.groupby(['MotifID'])['TF'].apply(lambda x: ', '.join(list(set(x)))).reset_index()
except:
    pass
motif_similarity_and_orthology_annot.index = motif_similarity_and_orthology_annot['MotifID']
motif_similarity_and_orthology_annot = pd.DataFrame(motif_similarity_and_orthology_annot['TF'])
motif_similarity_and_orthology_annot.columns = ['Motif_similarity_and_Orthology_annot']
# Combine
df = pd.concat([df_direct_annot, motif_similarity_annot, orthology_annot, motif_similarity_and_orthology_annot], axis=1, sort=False)

print(df.Direct_annot.unique())

print(df.Orthology_annot.unique())