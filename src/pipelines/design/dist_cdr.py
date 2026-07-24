'''
example: 
    - python app.py dist_cdr
functions:
    - calculate distance and analyze PPI
    - put meta data to table dist_cdr
'''
from src.helper import *
from bioomics import ProcessPickle

def retrieve(params, meta):
    source_table_name = params['source_table_name']
    query = f"""
        SELECT d.summary_id, d.relative_pkl, p.design_id, t.binder_id,
	        v.fragment, v.seq_from, v.seq_to, v.seq
        FROM {source_table_name}  d
        LEFT JOIN predict_summary s ON d.summary_id = s.summary_id
        LEFT JOIN predict         p ON s.predict_id=p.predict_id
        INNER JOIN design_dimer   t ON p.design_id = t.design_id
        RIGHT JOIN align_vfrag    v ON t.binder_id = v.chain_id
        WHERE s.ranking = 1
            AND d.relative_pkl IS NOT NULL
            AND v.fragment IN ('CDR1', 'CDR2', 'CDR3')
            AND d.summary_id NOT IN (
                SELECT summary_id  FROM {params['table_name']}
            )
    """
    rows = QueryComplex().list_data(query)
    meta['rows'] = len(rows)
    g = pd.DataFrame(rows).groupby(['summary_id', 'binder_id', 'relative_pkl'])
    return g


def build(g, params, meta):
    for (summary_id, binder_id, relative_pkl), sub in g:
        pfile = os.path.join(params['design_dir'], relative_pkl)
        if os.path.isfile(pfile):
            data = ProcessPickle(pfile).load_pickle([])[0]
            dist2 = data.get('dist2')
            # label CDRs
            cdr = pd.Series([None,] * len(dist2), index=dist2['res_no'])
            for _, row in sub.iterrows():
                frag = row['fragment']
                start = row['seq_from']
                end = row['seq_to']
                cdr.loc[start:end] = frag
            dist2['frag'] = cdr
            # build record
            _dist = dist2.dropna()
            for _, row in _dist.iterrows():
                yield (
                    summary_id,
                    binder_id,
                    row['frag'],
                    row['res_no'],
                    row['aa'],
                    row['value'],
                )

if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'chunk_size': 20,
        'dist_cutoff': 10,
        'span': 10,
    })

    # unrelaxed
    params['source_table_name'] = f'unrelaxed_dist'
    params['table_name'] = f'unrelaxed_dist_cdr'
    pdb_iter = retrieve(params, meta)
    record_iter = build(pdb_iter, params, meta)
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])


    # relaxed
    params['source_table_name'] = f'relaxed_dist'
    params['table_name'] = f'relaxed_dist_cdr'
    pdb_iter = retrieve(params, meta)
    record_iter = build(pdb_iter, params, meta)
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)