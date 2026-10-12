// ---- cutoff panel v4: four bases -- calibrated CV (deployable) / vertical-avg / pooled / in-sample ----
function renderCutoff(){
  var sel=document.getElementById('coCurve'); if(!sel) return;
  var ro=document.getElementById('coReadout'), oo=document.getElementById('coOpt'), lab=document.getElementById('coSlideLab');
  var empty=function(msg){ if(ro) ro.textContent=msg; if(oo) oo.textContent=''; if(plotReady('coChart')) Plotly.newPlot('coChart',[],{xaxis:{range:[0,1]},yaxis:{range:[0,1]}},{responsive:true,displayModeBar:false}); };
  if(!CO_LIST.length){ empty('No curve matches the current ROC filters — widen the ROC filters above.'); return; }
  var r=CO_LIST[Math.min(+sel.value||0, CO_LIST.length-1)];
  var p=parseFloat(document.getElementById('coSlider').value);
  document.getElementById('coCut').textContent=p.toFixed(2);
  var chance={x:[0,1],y:[0,1],mode:'lines',name:'Chance',line:{color:'#999',dash:'dot',width:1},hoverinfo:'skip'};
  var mode=coMode();
  function probBasis(data, title, color, note){
    if(lab) lab.textContent='Cutoff (probability)';
    var key=geneKeyFor(r); var e=(typeof data!=='undefined')?data[key]:null;
    if(!e){ empty('Curve not available for this config in this basis (gene-based feature sets only).'); return; }
    var auc=e[0], fpr=e[1], tp=e[2];
    var k=cutoffIndexFromP(p), sens=tp[k], spec=1-fpr[k];
    if(ro) ro.innerHTML='Sensitivity <b>'+sens.toFixed(3)+'</b> &nbsp;|&nbsp; Specificity <b>'+spec.toFixed(3)+'</b> &nbsp;&middot; at '+(mode==='calib'?'calibrated ':'')+'cutoff '+p.toFixed(2);
    var bestk=1,bestJ=-2; for(var i=1;i<=101;i++){ var J=tp[i]-fpr[i]; if(J>bestJ){bestJ=J;bestk=i;} }
    var optp=(101-bestk)/100, osens=tp[bestk], ospec=1-fpr[bestk];
    if(oo) oo.innerHTML='<b>Youden-optimal cutoff = '+optp.toFixed(2)+'</b> → sensitivity '+osens.toFixed(3)+', specificity '+ospec.toFixed(3)+' (J = '+bestJ.toFixed(3)+'; AUC '+auc.toFixed(3)+'). '+note;
    if(!plotReady('coChart')) return;
    var curve={x:fpr,y:tp,mode:'lines',name:title,line:{width:2,color:color},hovertemplate:'FPR %{x:.3f}<br>TPR %{y:.3f}<extra></extra>'};
    var cur={x:[fpr[k]],y:[tp[k]],mode:'markers',name:'cutoff '+p.toFixed(2),marker:{size:12,color:'#e67e22'},hovertemplate:'sensitivity '+sens.toFixed(3)+'<br>specificity '+spec.toFixed(3)+'<extra></extra>'};
    var opt={x:[fpr[bestk]],y:[tp[bestk]],mode:'markers',name:'Youden optimal',marker:{size:16,color:'#c0392b',symbol:'star'}};
    Plotly.newPlot('coChart',[curve,chance,cur,opt],{xaxis:{title:'False positive rate',range:[0,1]},yaxis:{title:'True positive rate',range:[0,1]},margin:{t:30,b:70,l:60,r:20},legend:{orientation:'h',y:1.15,yanchor:'bottom',x:0}},{responsive:true,displayModeBar:false});
  }
  if(mode==='vavg'){
    if(lab) lab.textContent='Specificity';
    var key=geneKeyFor(r); var e=(typeof FOLDAVG!=='undefined')?FOLDAVG[key]:null;
    if(!e){ empty('Per-fold averaged curve not available for this config (gene-based feature sets only).'); return; }
    var auc=e[0], tpr=e[1];
    var spec=p, fprv=1-spec, idx=Math.max(0,Math.min(100,Math.round(fprv*100))), sens=tpr[idx];
    if(ro) ro.innerHTML='Sensitivity <b>'+sens.toFixed(3)+'</b> &nbsp;|&nbsp; Specificity <b>'+spec.toFixed(2)+'</b> &nbsp;&middot; FPR '+fprv.toFixed(2);
    var bestk=0,bestJ=-2; for(var i=0;i<=100;i++){ var J=tpr[i]-(i/100); if(J>bestJ){bestJ=J;bestk=i;} }
    var ofpr=bestk/100, ospec=1-ofpr, osens=tpr[bestk];
    if(oo) oo.innerHTML='<b>Youden-optimal operating point: specificity '+ospec.toFixed(2)+', sensitivity '+osens.toFixed(3)+'</b> (J = '+bestJ.toFixed(3)+'; averaged AUC '+auc.toFixed(3)+'). Honest trade-off, but NO single probability cutoff (per-fold thresholds differ) — use the Calibrated CV basis to get a deployable cutoff.';
    if(!plotReady('coChart')) return;
    var gx=[]; for(var i=0;i<=100;i++) gx.push(i/100);
    var curve={x:gx,y:tpr,mode:'lines',name:'Per-fold averaged ROC',line:{width:2,color:'#1e7a3f'},hovertemplate:'FPR %{x:.2f}<br>sens %{y:.3f}<extra></extra>'};
    var cur={x:[fprv],y:[sens],mode:'markers',name:'spec '+spec.toFixed(2),marker:{size:12,color:'#e67e22'},hovertemplate:'sensitivity '+sens.toFixed(3)+'<br>specificity '+spec.toFixed(2)+'<extra></extra>'};
    var opt={x:[ofpr],y:[osens],mode:'markers',name:'Youden optimal',marker:{size:16,color:'#c0392b',symbol:'star'}};
    Plotly.newPlot('coChart',[curve,chance,cur,opt],{xaxis:{title:'False positive rate (1 − specificity)',range:[0,1]},yaxis:{title:'Sensitivity',range:[0,1]},margin:{t:30,b:70,l:60,r:20},legend:{orientation:'h',y:1.15,yanchor:'bottom',x:0}},{responsive:true,displayModeBar:false});
  } else if(mode==='insample'){
    probBasis(typeof INSAMPLE!=='undefined'?INSAMPLE:undefined,'Full-model in-sample ROC','#8e44ad','<b>IN-SAMPLE / OPTIMISTIC:</b> full-data model scored on its OWN training data. Shows the deployable probability scale, but real-world sens/spec are LOWER — use the Calibrated CV basis; validate externally.');
  } else if(mode==='calib'){
    probBasis(typeof CALIB!=='undefined'?CALIB:undefined,'Calibrated CV ROC','#16a085','<b>RECOMMENDED / DEPLOYABLE:</b> honest out-of-fold CV performance on a CALIBRATED probability scale. Apply this cutoff to the full-data model with the calibration baked into predict_new_cohort.R, then validate externally.');
  }
  if(mode==='pooled'){
    if(lab) lab.textContent='Cutoff (probability)';
    var fpr=parseCurve(r.fpr), tp=parseCurve(r.tpr);
    var k=cutoffIndexFromP(p), sens=tp[k], spec=1-fpr[k];
    if(ro) ro.innerHTML='Sensitivity <b>'+sens.toFixed(3)+'</b> &nbsp;|&nbsp; Specificity <b>'+spec.toFixed(3)+'</b> &nbsp;&middot; at probability cutoff '+p.toFixed(2);
    var bestk=1,bestJ=-2; for(var i=1;i<=101;i++){ var J=tp[i]-fpr[i]; if(J>bestJ){bestJ=J;bestk=i;} }
    var optp=(101-bestk)/100, osens=tp[bestk], ospec=1-fpr[bestk];
    if(oo) oo.innerHTML='<b>Youden-optimal cutoff = '+optp.toFixed(2)+'</b> → sensitivity '+osens.toFixed(3)+', specificity '+ospec.toFixed(3)+' (J = '+bestJ.toFixed(3)+'). Pooled CV on a single global RAW probability threshold — can sit below the diagonal for small-event configs; prefer the Calibrated CV basis.';
    if(!plotReady('coChart')) return;
    var curve={x:fpr,y:tp,mode:'lines',name:'Pooled ROC',line:{width:2,color:'#245f9d'},hovertemplate:'FPR %{x:.3f}<br>TPR %{y:.3f}<extra></extra>'};
    var cur={x:[fpr[k]],y:[tp[k]],mode:'markers',name:'cutoff '+p.toFixed(2),marker:{size:12,color:'#e67e22'},hovertemplate:'sensitivity '+sens.toFixed(3)+'<br>specificity '+spec.toFixed(3)+'<extra></extra>'};
    var opt={x:[fpr[bestk]],y:[tp[bestk]],mode:'markers',name:'Youden optimal',marker:{size:16,color:'#c0392b',symbol:'star'}};
    Plotly.newPlot('coChart',[curve,chance,cur,opt],{xaxis:{title:'False positive rate',range:[0,1]},yaxis:{title:'True positive rate',range:[0,1]},margin:{t:30,b:70,l:60,r:20},legend:{orientation:'h',y:1.15,yanchor:'bottom',x:0}},{responsive:true,displayModeBar:false});
  }
}
