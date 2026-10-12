// ---- cutoff panel v2: toggle between Pooled and Per-fold vertical-averaged ----
function coMode(){ var e=document.getElementById('coMode'); return e? e.value : 'pooled'; }
function renderCutoff(){
  var sel=document.getElementById('coCurve'); if(!sel) return;
  var ro=document.getElementById('coReadout'), oo=document.getElementById('coOpt'), lab=document.getElementById('coSlideLab');
  var empty=function(msg){ if(ro) ro.textContent=msg; if(oo) oo.textContent=''; if(plotReady('coChart')) Plotly.newPlot('coChart',[],{xaxis:{range:[0,1]},yaxis:{range:[0,1]}},{responsive:true,displayModeBar:false}); };
  if(!CO_LIST.length){ empty('No curve matches the current ROC filters — widen the ROC filters above.'); return; }
  var r=CO_LIST[Math.min(+sel.value||0, CO_LIST.length-1)];
  var p=parseFloat(document.getElementById('coSlider').value);
  document.getElementById('coCut').textContent=p.toFixed(2);
  var chance={x:[0,1],y:[0,1],mode:'lines',name:'Chance',line:{color:'#999',dash:'dot',width:1},hoverinfo:'skip'};
  if(coMode()==='vavg'){
    if(lab) lab.textContent='Specificity';
    var key=geneKeyFor(r);
    var e=(typeof FOLDAVG!=='undefined')?FOLDAVG[key]:null;
    if(!e){ empty('Per-fold averaged curve not available for this config (gene-based feature sets only: fpkm / fpkm_clinics).'); return; }
    var auc=e[0], tpr=e[1];                 // fpr grid is fixed: i/100, i=0..100
    var spec=p, fprv=1-spec, idx=Math.max(0,Math.min(100,Math.round(fprv*100))), sens=tpr[idx];
    if(ro) ro.innerHTML='at specificity <b>'+spec.toFixed(2)+'</b>: sensitivity <b>'+sens.toFixed(3)+'</b> (FPR '+fprv.toFixed(2)+')';
    var bestk=0,bestJ=-2; for(var i=0;i<=100;i++){ var J=tpr[i]-(i/100); if(J>bestJ){bestJ=J;bestk=i;} }
    var ofpr=bestk/100, ospec=1-ofpr, osens=tpr[bestk];
    if(oo) oo.innerHTML='<b>Youden-optimal operating point: specificity '+ospec.toFixed(2)+', sensitivity '+osens.toFixed(3)+'</b> (J = '+bestJ.toFixed(3)+'). Vertical-averaged across the 10×100 CV folds; averaged AUC '+auc.toFixed(3)+'. No single probability cutoff — the threshold that gives this specificity differs per fold.';
    if(!plotReady('coChart')) return;
    var gx=[]; for(var i=0;i<=100;i++) gx.push(i/100);
    var curve={x:gx,y:tpr,mode:'lines',name:'Per-fold averaged ROC',line:{width:2,color:'#1e7a3f'},hovertemplate:'FPR %{x:.2f}<br>sens %{y:.3f}<extra></extra>'};
    var cur={x:[fprv],y:[sens],mode:'markers',name:'spec '+spec.toFixed(2),marker:{size:12,color:'#e67e22'}};
    var opt={x:[ofpr],y:[osens],mode:'markers',name:'Youden optimal',marker:{size:16,color:'#c0392b',symbol:'star'}};
    Plotly.newPlot('coChart',[curve,chance,cur,opt],{xaxis:{title:'False positive rate (1 − specificity)',range:[0,1]},yaxis:{title:'Sensitivity',range:[0,1]},margin:{t:20},legend:{orientation:'h'}},{responsive:true,displayModeBar:false});
  } else {
    if(lab) lab.textContent='Cutoff (probability)';
    var fpr=parseCurve(r.fpr), tp=parseCurve(r.tpr);
    var k=cutoffIndexFromP(p), sens=tp[k], spec=1-fpr[k];
    if(ro) ro.innerHTML='at cutoff <b>'+p.toFixed(2)+'</b>: sensitivity <b>'+sens.toFixed(3)+'</b>, specificity <b>'+spec.toFixed(3)+'</b>';
    var bestk=1,bestJ=-2; for(var i=1;i<=101;i++){ var J=tp[i]-fpr[i]; if(J>bestJ){bestJ=J;bestk=i;} }
    var optp=(101-bestk)/100, osens=tp[bestk], ospec=1-fpr[bestk];
    if(oo) oo.innerHTML='<b>Youden-optimal cutoff = '+optp.toFixed(2)+'</b> → sensitivity '+osens.toFixed(3)+', specificity '+ospec.toFixed(3)+' (J = '+bestJ.toFixed(3)+'). Pooled CV predictions on a single global probability threshold — can sit below the diagonal for small-event configs; use the per-fold averaged basis for those.';
    if(!plotReady('coChart')) return;
    var curve={x:fpr,y:tp,mode:'lines',name:'Pooled ROC',line:{width:2,color:'#245f9d'},hovertemplate:'FPR %{x:.3f}<br>TPR %{y:.3f}<extra></extra>'};
    var cur={x:[fpr[k]],y:[tp[k]],mode:'markers',name:'cutoff '+p.toFixed(2),marker:{size:12,color:'#e67e22'}};
    var opt={x:[fpr[bestk]],y:[tp[bestk]],mode:'markers',name:'Youden optimal',marker:{size:16,color:'#c0392b',symbol:'star'}};
    Plotly.newPlot('coChart',[curve,chance,cur,opt],{xaxis:{title:'False positive rate',range:[0,1]},yaxis:{title:'True positive rate',range:[0,1]},margin:{t:20},legend:{orientation:'h'}},{responsive:true,displayModeBar:false});
  }
}
