from pathlib import Path
import json
R=Path(__file__).resolve().parents[1]
c=json.loads((R/'results/certificate_d69.json').read_text('utf8'))
m=json.loads((R/'publication/zenodo_metadata_blueprint.json').read_text('utf8'))
assert c['status']=='PASS' and c['dimension']==m['scientific_dimension']==69
assert c['vertices']==m['norm_one_vertices']==74
assert c['pairs']==m['exact_labelled_squared_pairs']==2701
assert c['affine_rank']==c['centered_gram_rank']==69
assert c['eigenvalue_37_multiplicity']==m['model_eigenvalue_37_multiplicity']==63
assert c['eigenvalue_73_over_2_multiplicity']==m['model_eigenvalue_73_over_2_multiplicity']==4
assert c['parent_tag_commit']=='a431c711973df7ea1e14489ad0699b4b715d745f'
assert c['parent_doi']==m['scientific_parent_doi']=='10.5281/zenodo.24433099'
assert m['related_identifiers']==[{'identifier':c['parent_doi'],'relation':'isDerivedFrom','resource_type':'publication-preprint','scheme':'doi'}]
assert len(m['keywords'])==9 and len(set(m['keywords']))==9
assert 'affine Euclidean dimension sixty-nine' in m['keywords'] and 'sixty-threefold ideal eigenvalue' in m['keywords']
assert not set(m['keywords']).intersection({'affine Euclidean dimension sixty-eight','sixty-twofold ideal eigenvalue','affine Euclidean dimension seventy','sixty-fourfold ideal eigenvalue'})
assert len(m['references'])==2 and c['parent_doi'] in m['references'][0] and '10.5281/zenodo.24407952' in m['references'][1]
assert all(x in (R/'publication/POLYBRUTUS_69D_NOTE.md').read_text('utf8') for x in ('## References',c['parent_doi'],'10.5281/zenodo.24407952'))
print('D69 METADATA 9 KEYWORDS / 2 REFERENCES / DIRECT D68 DOI PASS')
