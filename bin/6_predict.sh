#! /usr/bin/bash
#
# NOTE: af_params, af_msa, af_pred, af_summary, af_plddt, af_rmsd and
# af_tmalign moved to structure_profiling:
#   see ~/bio/structure_profiling/bin/6_predict.sh
#
# check pending combos with:
#   SELECT a.combo_id AS input_id, STRING_AGG(a.category, ',')
#   FROM label_complex a
#   LEFT JOIN meta_pdb p ON a.pdb_id = p.pdb_id
#   WHERE p.pdb_group = 0 AND a.combo = 5 AND NOT EXISTS (
#       SELECT 1 FROM af_params b WHERE a.combo_id = b.input_id)
#   GROUP BY a.combo_id;



