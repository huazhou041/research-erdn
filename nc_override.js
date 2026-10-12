function ncMadLookup(){
  var m={};
  LONG.filter(function(r){return r.model==='glm_mad';}).forEach(function(r){
    m[[r.set,r.cohort,r.datatype,r.featureset,r.subset,r.varcut].join('|')]=r.mean;
  });
  return m;
}
function ncRowsWithOwn(){
  var setF=document.getElementById('ncSet').value, taskF=document.getElementById('ncTask').value,
      fsF=document.getElementById('ncFs').value, sortF=document.getElementById('ncSort').value;
  var mad=ncMadLookup();
  var rows=LONG.filter(function(r){return r.model==='glm_natcommad'
     && (setF==='all'||r.set===setF) && (taskF==='all'||r.subset===taskF) && (fsF==='all'||r.featureset===fsF);})
    .map(function(r){
      var own=mad[[r.set,r.cohort,r.datatype,r.featureset,r.subset,r.varcut].join('|')];
      if(own===undefined) own=null;
      var delta=(own!=null&&r.mean!=null)?(r.mean-own):null;
      return Object.assign({},r,{own:own,delta:delta});
    });
  if(sortF==='mean') rows.sort(function(a,b){return (b.mean==null?-1:b.mean)-(a.mean==null?-1:a.mean);});
  else if(sortF==='delta') rows.sort(function(a,b){return (b.delta==null?-99:b.delta)-(a.delta==null?-99:a.delta);});
  else rows.sort(function(a,b){return String(a.cohort).localeCompare(String(b.cohort))||((b.mean==null?-1:b.mean)-(a.mean==null?-1:a.mean));});
  return rows;
}
function renderNatcommad(){
  if(!document.getElementById('ncTable')) return;
  var rows=ncRowsWithOwn();
  document.querySelector('#ncTable tbody').innerHTML=rows.map(function(r){
    var own=r.own==null?'':r.own.toFixed(3);
    var nc=r.mean==null?'':r.mean.toFixed(3);
    var dcol='';
    if(r.delta!=null){var c=r.delta>=0?'#1e7a3f':'#b03030';var s=r.delta>=0?'+':'';dcol='<span style="color:'+c+'">'+s+r.delta.toFixed(3)+'</span>';}
    return '<tr><td>'+r.set+'</td><td>'+r.cohort+'</td><td>'+r.datatype+'</td><td>'+r.featureset+'</td><td>'+r.subset+'</td><td>'+r.varcut+'</td><td>'+own+'</td><td>'+nc+'</td><td>'+dcol+'</td><td>'+(r.n==null?'':r.n)+'</td></tr>';
  }).join('');
  var best=rows.length?rows.reduce(function(m,r){return (r.mean==null?-1:r.mean)>(m.mean==null?-1:m.mean)?r:m;}):null;
  var wins=rows.filter(function(r){return r.delta!=null&&r.delta>0;}).length;
  var tot=rows.filter(function(r){return r.delta!=null;}).length;
  document.getElementById('ncSummary').textContent=rows.length+' configs ('+(useSmoking()?'ever-smoker':'full')+'); nature_com genes beat own-cohort MAD in '+wins+'/'+tot+' matched configs'+(best?('; best ncMAD '+best.mean.toFixed(3)+' @ '+best.cohort+'/'+best.featureset+'/top-'+best.varcut):'');
}
function downloadNatcommad(){
  var rows=ncRowsWithOwn();
  var hdr=['set','cohort','datatype','featureset','task','top_n_genes','auc_own_cohort_mad','auc_naturecom_mad','delta_nc_minus_own','sd_nc','n'];
  var lines=[hdr.join(',')];
  rows.forEach(function(r){
    var vals=[r.set,r.cohort,r.datatype,r.featureset,r.subset,r.varcut,
      r.own==null?'':r.own, r.mean==null?'':r.mean, r.delta==null?'':r.delta, r.sd==null?'':r.sd, r.n==null?'':r.n];
    lines.push(vals.map(function(x){return '"'+String(x).replace(/"/g,'""')+'"';}).join(','));
  });
  var blob=new Blob([lines.join('\n')],{type:'text/csv'});
  var a=document.createElement('a'); a.href=URL.createObjectURL(blob);
  a.download='natcommad_vs_owncohort_'+(useSmoking()?'eversmoker':'full')+'.csv'; a.click(); URL.revokeObjectURL(a.href);
}
