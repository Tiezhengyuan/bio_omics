'''
example:
    python src/collect/distance.py
'''
from src.helper import *
from bioomics import DistanceAnalyze, DistanceProcess

def get_meta(table_name):
    query = f"""
        select D.pdb_id, D.combo_id, C.chain_combo, D.relative_pkl
        from {table_name} D
        left join meta_combo2 C on D.combo_id = C.combo_id
        where D.pdb_id in (
            select distinct pdb_id from view_antibody
        )
        AND D.relative_pkl is not null
    ;"""
    rows = QueryComplex(True).list_data(query)
    df = pd.DataFrame(rows)
    info = {
        'table_name': table_name,
        'pdb': df['pdb_id'].nunique(),
        'combo': df['combo_id'].nunique(),
        'regions': len(rows),
    }
    print(info)
    return rows

if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': False,
        'data_dir': os.getenv('data_dir'),
    }
    outdir = Path(args['data_dir']) / 'collect' / 'distance'
    outdir.mkdir(parents=True, exist_ok=True)

    tables = ['combo2_distance', 'pdb_distance']
    spans = [10, 5, 1,]
    cutoffs = [10, 5,]
    combinations = list(itertools.product(tables, spans, cutoffs))
    for (table_name, span, cutoff) in combinations:
        dist_meta = get_meta(table_name)
        args['max_span'] = span
        args['dist_cutoff'] = cutoff
        outfile = outdir / f"{table_name}_cutoff{cutoff}_span{span}.csv"
        if args['overwrite'] or not outfile.is_file():
            df = DistanceProcess(args).retrieve_interface_seq(dist_meta)
            df.to_csv(outfile, index=None)
            info = {
                'table_name': table_name,
                'span': span,
                'cutoff': cutoff,
                'pdb': df['pdb_id'].nunique(),
                'combo': df['combo_id'].nunique(),
                'out_csv': outfile,
            }
            print(info)
        else:
            pass
        