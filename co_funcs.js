// ---- Operating point & optimal cutoff (Youden's J) ----
// ROC curves were computed at thresholds c(Inf, seq(1,0,by=-0.01), -Inf): 103 points.
// JS index k (1..101) <-> probability cutoff (101-k)/100 ; predict positive if prob>=cutoff.
// sensitivity = tpr[k]; specificity = 1 - fpr[k]. Uses RAW (un-flipped) curve.
var CO_LIST = [];
function cutoffIndexFromP(p){ return Math.max(1, Math.min(101, 101 - Math.round(p*100))); }
function populateCutoff(){
  var sel=document.getElementById('coCurve'); if(!sel) return;
  CO_LIST = (typeof rocFiltered==='function'? rocFiltered() : []).slice().sort(function(a,b){return b.auc-a.auc;}).slice(0,300);
  var cur = sel.value;
  sel.innerHTML = CO_LIST.map(function(r,i){ return '<option value="'+i+'">'+labelModel(r.model)+' · '+r.set+' · '+r.cohort+' · '+r.subset+' · '+r.featureset+'/'+r.datatype+' · v'+r.varcut+' · AUC '+fmt(r.auc)+'</option>'; }).join('');
  if(cur && +cur < CO_LIST.length) sel.value=cur;
  renderCutoff();
}
function renderCutoff(){
  var sel=document.getElementById('coCurve'); if(!sel) return;
  var ro=document.getElementById('coReadout'), oo=document.getElementById('coOpt');
  if(!CO_LIST.length){ if(ro) ro.textContent='No curve matches the current ROC filters — widen the ROC filters above.'; if(oo) oo.textContent=''; if(plotReady('coChart')) Plotly.newPlot('coChart',[],{xaxis:{range:[0,1]},yaxis:{range:[0,1]}},{responsive:true,displayModeBar:false}); return; }
  var r=CO_LIST[Math.min(+sel.value||0, CO_LIST.length-1)];
  var fpr=parseCurve(r.fpr), tpr=parseCurve(r.tpr);
  var p=parseFloat(document.getElementById('coSlider').value);
  document.getElementById('coCut').textContent=p.toFixed(2);
  var k=cutoffIndexFromP(p), sens=tpr[k], spec=1-fpr[k];
  if(ro) ro.innerHTML='at cutoff <b>'+p.toFixed(2)+'</b>: sensitivity <b>'+sens.toFixed(3)+'</b>, specificity <b>'+spec.toFixed(3)+'</b>';
  var bestk=1,bestJ=-2; for(var i=1;i<=101;i++){ var J=tpr[i]-fpr[i]; if(J>bestJ){bestJ=J;bestk=i;} }
  var optp=(101-bestk)/100, osens=tpr[bestk], ospec=1-fpr[bestk];
  if(oo) oo.innerHTML='<b>Youden-optimal cutoff = '+optp.toFixed(2)+'</b> → sensitivity '+osens.toFixed(3)+', specificity '+ospec.toFixed(3)+' (Youden J = '+bestJ.toFixed(3)+'). AUC '+fmt(r.auc)+', n='+r.n_pred+(r.n_pos_pred!=null?(' ('+r.n_pos_pred+' recur / '+r.n_neg_pred+' non-recur)'):'')+'. Pooled CV predictions — illustrative operating point, not a validated threshold.';
  if(!plotReady('coChart')) return;
  var curve={x:fpr,y:tpr,mode:'lines',name:'ROC',line:{width:2,color:'#245f9d'},hovertemplate:'FPR %{x:.3f}<br>TPR %{y:.3f}<extra></extra>'};
  var chance={x:[0,1],y:[0,1],mode:'lines',name:'Chance',line:{color:'#999',dash:'dot',width:1},hoverinfo:'skip'};
  var cur={x:[fpr[k]],y:[tpr[k]],mode:'markers',name:'cutoff '+p.toFixed(2),marker:{size:12,color:'#e67e22'},hovertemplate:'cutoff '+p.toFixed(2)+'<br>sens %{y:.3f}<br>spec '+spec.toFixed(3)+'<extra></extra>'};
  var opt={x:[fpr[bestk]],y:[tpr[bestk]],mode:'markers',name:'Youden optimal',marker:{size:16,color:'#c0392b',symbol:'star'},hovertemplate:'optimal cutoff '+optp.toFixed(2)+'<br>sens '+osens.toFixed(3)+'<br>spec '+ospec.toFixed(3)+'<extra></extra>'};
  Plotly.newPlot('coChart',[curve,chance,cur,opt],{xaxis:{title:'False positive rate (1 − specificity)',range:[0,1]},yaxis:{title:'True positive rate (sensitivity)',range:[0,1]},margin:{t:20},legend:{orientation:'h'}},{responsive:true,displayModeBar:false});
}
