import json

import pandas as pd

from common import DATA, ROOT

db = pd.read_csv(DATA / "database.csv")
cols = ["db_id", "system", "dopant", "source_id", "lab_group", "year", "doi_or_arxiv", "implant_temp_C", "anneal_temp_C",
        "anneal_time_min", "implanted_conc_cm3", "dose_cm2", "measurement_technique", "activation_pct", "definition_class",
        "value_basis", "text_basis", "origin"]
recs = json.loads(db[cols].to_json(orient="records"))
html = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Acceptor activation database</title>
<style>
body{font-family:system-ui,sans-serif;margin:16px;max-width:1200px}
table{border-collapse:collapse;font-size:12px;width:100%}
th,td{border-bottom:1px solid #ccc;padding:3px 6px;text-align:left;white-space:nowrap}
th{position:sticky;top:0;background:#f4f4f4;cursor:pointer}
.controls{display:flex;gap:12px;flex-wrap:wrap;margin:10px 0}
.wrap{overflow:auto;max-height:70vh}
</style></head><body>
<h1>Acceptor activation database</h1>
<p>Literature-compiled electrical activation of ion-implanted acceptors. Each row has a DOI. Click a column header to sort.</p>
<div class="controls">
<label>System <select id="system"><option value="">all</option></select></label>
<label>Dopant <select id="dopant"><option value="">all</option></select></label>
<label>Value basis <select id="basis"><option value="">all</option></select></label>
<label>Search <input id="q" type="search" placeholder="source, DOI, text"></label>
<span id="count"></span>
</div>
<div class="wrap"><table id="t"><thead></thead><tbody></tbody></table></div>
<script>
const DATA = __DATA__;
const cols = Object.keys(DATA[0]);
let sortKey = null, asc = true;
function fill(id, key){const s=document.getElementById(id);[...new Set(DATA.map(r=>r[key]).filter(v=>v!==null))].sort().forEach(v=>{const o=document.createElement('option');o.value=v;o.textContent=v;s.appendChild(o)})}
fill('system','system');fill('dopant','dopant');fill('basis','value_basis');
const thead=document.querySelector('#t thead');
thead.innerHTML='<tr>'+cols.map(c=>'<th data-k="'+c+'">'+c+'</th>').join('')+'</tr>';
thead.addEventListener('click',e=>{const k=e.target.dataset.k;if(!k)return;asc=(sortKey===k)?!asc:true;sortKey=k;render()});
['system','dopant','basis','q'].forEach(id=>document.getElementById(id).addEventListener('input',render));
function render(){
  const sys=document.getElementById('system').value, dop=document.getElementById('dopant').value,
        bas=document.getElementById('basis').value, q=document.getElementById('q').value.toLowerCase();
  let rows=DATA.filter(r=>(!sys||r.system===sys)&&(!dop||r.dopant===dop)&&(!bas||r.value_basis===bas)&&
    (!q||JSON.stringify(r).toLowerCase().includes(q)));
  if(sortKey){rows=rows.slice().sort((a,b)=>{const x=a[sortKey],y=b[sortKey];if(x===y)return 0;if(x===null)return 1;if(y===null)return -1;
    return (x<y?-1:1)*(asc?1:-1)})}
  document.getElementById('count').textContent=rows.length+' of '+DATA.length+' rows';
  document.querySelector('#t tbody').innerHTML=rows.map(r=>'<tr>'+cols.map(c=>'<td>'+(r[c]===null?'':r[c])+'</td>').join('')+'</tr>').join('');
}
render();
</script></body></html>
""".replace("__DATA__", json.dumps(recs))
(ROOT / "explorer.html").write_text(html)
print(len(recs), "rows written to explorer.html")
