/* Saved R600 joint-state transport results. No calculation is fitted in the page. */
(function () {
  'use strict';
  const data = JSON.parse(document.getElementById('joint-state-data').textContent);
  const red = data.reduction, ox = data.reoxidation;
  const oxProfiles = ox.filter(q => q.depth_nm);
  const color = ['#0072B2', '#D55E00'];
  const el = id => document.getElementById(id);
  const esc = s => String(s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
  function plot(id, letter, xlabel, ylabel, series, opt={}) {
    // Coordinates are points on a 5 × 5 inch (360 pt) publication panel.
    const L=77.54,T=37.40,W=262.80,H=262.80;
    const tx=x=>opt.depth ? Math.log10(1+x) : x;
    const ty=y=>opt.log ? Math.log10(Math.max(1e-10,y)) : y;
    const xx=series.flatMap(s=>s.x).filter(Number.isFinite);
    const yy=series.flatMap(s=>s.y).filter(Number.isFinite);
    const xmin=0,xmax=tx(Math.max(...xx)), ymin=opt.log ? Math.floor(Math.min(...yy.map(ty))) : 0;
    const ymax=opt.log ? Math.ceil(Math.max(...yy.map(ty))) : Math.ceil(Math.max(...yy)*1.05/25)*25;
    const X=x=>L+W*(tx(x)-xmin)/(xmax-xmin||1),Y=y=>T+H-H*(ty(y)-ymin)/(ymax-ymin||1);
    let svg='<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 360 360" width="5in" height="5in" role="img" aria-label="'+esc(ylabel)+'" style="background:white;font-family:Helvetica,Arial,sans-serif;max-width:100%;height:auto">';
    svg+='<defs><clipPath id="clip-'+id+'"><rect x="'+L+'" y="'+T+'" width="'+W+'" height="'+H+'"/></clipPath></defs><text x="12" y="52" font-size="18" font-weight="bold">'+letter+'</text><path d="M'+L+' '+T+'v'+H+'h'+W+'v-'+H+'Z" fill="none" stroke="black" stroke-width="2"/>';
    const xt=opt.depth ? [0,1,10,100,450] : [0,.25,.5,.75,1].map(q=>Math.round(q*Math.max(...xx)));
    for(const x of xt) svg+='<path d="M'+X(x)+' '+(T+H)+'v-4" stroke="black" stroke-width="2"/><text x="'+X(x)+'" y="'+(T+H+21)+'" text-anchor="middle" font-size="15">'+x+'</text>';
    const ystep=opt.log ? Math.max(1,Math.ceil((ymax-ymin)/6)) : Math.max(25,Math.ceil(ymax/100)*25);
    const yt=Array.from({length:Math.floor((ymax-ymin)/ystep)+1},(_,i)=>ymin+i*ystep);
    for(const y of yt) { const raw=opt.log?10**y:y;svg+='<path d="M'+L+' '+Y(raw)+'h4" stroke="black" stroke-width="2"/><text x="'+(L-9)+'" y="'+(Y(raw)+5)+'" text-anchor="end" font-size="15">'+(opt.log?'10<tspan baseline-shift="super" font-size="10">'+Math.round(y)+'</tspan>':Number(y.toPrecision(2)))+'</text>'; }
    series.forEach((s,k)=>{
      if(s.point) svg+='<circle data-measured="94" cx="'+X(s.x[0])+'" cy="'+Y(s.y[0])+'" r="6" fill="white" stroke="black" stroke-width="1.2"/>';
      else {
        // The wiki plots distinct atomic planes as points, not a continuous profile.
        const first=opt.depth?12:0;
        for(let i=0;i<first;i++) svg+='<circle clip-path="url(#clip-'+id+')" data-atomic-plane="'+i+'" cx="'+X(s.x[i])+'" cy="'+Y(s.y[i])+'" r="3" fill="'+(s.color||color[k])+'"/>';
        svg+='<path clip-path="url(#clip-'+id+')" data-series="'+esc(s.name)+'" d="'+s.x.slice(first).map((x,i)=>(i?'L':'M')+X(x).toFixed(2)+' '+Y(s.y[i+first]).toFixed(2)).join(' ')+'" fill="none" stroke="'+(s.color||color[k])+'" stroke-width="2"/>';
      }
      svg+='<text x="'+(L+(opt.depth?65:9))+'" y="'+(opt.depth?T+H-35+k*18:T+20+k*18)+'" font-size="14" fill="'+(s.point?'black':s.color||color[k])+'">'+esc(s.name)+'</text>';
    });
    svg+='<text x="'+(L+W/2)+'" y="348" text-anchor="middle" font-size="18">'+esc(xlabel)+'</text><text transform="translate(22 '+(T+H/2)+') rotate(-90)" text-anchor="middle" font-size="18">'+esc(ylabel)+'</text></svg>';
    el(id).innerHTML=svg;
  }
  function profiles() {
    const a=red[+el('jrReduction').value],b=oxProfiles[+el('jrReoxidation').value];
    for(const [id,letter,key,label] of [['jrB','b','vacancy_pct','Vacancy fraction (%)'],['jrC','c','Ti3_pct','Ti³⁺ fraction (%)']])
      plot(id,letter,'Depth (nm)',label,[{x:a.depth_nm,y:a[key],name:'Reduction, '+a.time_s+' s'},{x:b.depth_nm,y:b[key],name:'Reoxidation, '+b.time_s+' s'}],{depth:true,log:true});
  }
  for(const [id,rows] of [['jrReduction',red],['jrReoxidation',oxProfiles]]) {
    el(id).innerHTML=rows.map((q,i)=>'<option value="'+i+'">'+q.time_s+' s</option>').join('');
    el(id).value=rows.length-1;el(id).addEventListener('change',profiles);
  }
  plot('jrA','a','Reduction time (s)','Vacancies (µmol O g⁻¹)',[{x:[0,...red.map(q=>q.time_s)],y:[0,...red.map(q=>q.vacancy_total_umol_g)],name:'Calculated from imposed source'},{x:[1800],y:[94],point:true,name:'Measured R600 inventory'}]);
  plot('jrD','d','Reoxidation time (s)','Oxygen amount (µmol O g⁻¹)',[{x:ox.map(q=>q.time_s),y:ox.map(q=>q.vacancy_total_umol_g),name:'Remaining vacancies'},{x:ox.map(q=>q.time_s),y:ox.map(q=>ox[0].vacancy_total_umol_g-q.vacancy_total_umol_g),name:'Vacancies filled'}]);
  profiles();
  el('jrCases').innerHTML='<table><thead><tr><th>Case</th><th>Outer barrier (eV)</th><th>Interior barrier (eV)</th><th>Outer depth (nm)</th><th>Reoxidation end (s)</th><th>Status</th></tr></thead><tbody>'+data.cases.map(q=>'<tr>'+[q.tag,q.parameters.outer_barrier_eV,q.parameters.interior_barrier_eV,q.parameters.outer_depth_nm,Number(q.end_s.toFixed(2)),q.status].map(v=>'<td>'+esc(v)+'</td>').join('')+'</tr>').join('')+'</tbody></table>';
  const rows=[];
  for(const [history,checks] of Object.entries(data.checks)) for(const [key,value] of Object.entries(checks)) if(typeof value==='number' && !key.startsWith('relaxed_')) rows.push('<tr><td>'+esc(history)+'</td><td>'+esc(key.replace(/_/g,' '))+'</td><td>'+value.toExponential(3)+'</td></tr>');
  el('jrValidation').innerHTML='<div style="overflow:auto"><table><thead><tr><th>History</th><th>Check</th><th>Value</th></tr></thead><tbody>'+rows.join('')+'</tbody></table></div>';
  document.querySelectorAll('[data-jr-download]').forEach(b=>b.addEventListener('click',()=>{
    const blob=new Blob([el(b.dataset.jrDownload).innerHTML],{type:'image/svg+xml'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='R600-'+b.dataset.jrDownload+'.svg';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  }));
  document.querySelectorAll('[data-jr-download]').forEach(b=>{
    b.textContent='SVG';
    const p=document.createElement('button');p.type='button';p.textContent='PNG · 600 dpi';b.after(p);
    p.addEventListener('click',()=>{
      const url=URL.createObjectURL(new Blob([el(b.dataset.jrDownload).innerHTML],{type:'image/svg+xml'}));
      const img=new Image();img.onload=()=>{
        const canvas=document.createElement('canvas');canvas.width=canvas.height=3000;
        canvas.getContext('2d').drawImage(img,0,0,3000,3000);URL.revokeObjectURL(url);
        canvas.toBlob(blob=>{const u=URL.createObjectURL(blob),a=document.createElement('a');a.href=u;a.download='R600-'+b.dataset.jrDownload+'-600dpi.png';a.click();setTimeout(()=>URL.revokeObjectURL(u),1000);});
      };img.src=url;
    });
  });
  window.TimeDependentDistribution={data,profiles};
})();
