"""Build paired peak-search probes using the local deterministic audit setup.

Requires matrix-audit and convergence-probe/signed from the earlier audit
builders. Keeps convergence, covariance and complex matrix subtraction at the
selected baseline. Use the saved peak-only source after later numerical changes.
"""
import argparse,hashlib,json,shlex,subprocess
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__);p.add_argument('build',type=Path)
p.add_argument('--base-ref',default='485d761',help='Commit used to build the saved baseline audit')
p.add_argument('--source',type=Path,help='Saved peak-only production source, before later covariance changes')
a=p.parse_args()
build=a.build.resolve();repo=Path(__file__).resolve().parents[1]
baseline_commit=subprocess.check_output(['git','rev-parse',a.base_ref],text=True,cwd=repo).strip()
old=subprocess.check_output(['git','show',baseline_commit+':src/cal_ccf.cpp'],text=True,cwd=repo)
production_source=a.source.resolve() if a.source else repo/'src/cal_ccf.cpp'
current=production_source.read_text()
def search(source):return source.rindex('int search_max(')
assert old[:search(old)]==current[:search(current)],'Only peak search should change; use --source with the saved peak-only production snapshot'
base=(build/'convergence-probe/signed/cal_ccf.cpp').read_text()
assert base[search(base):]==old[search(old):],'Rebuild baseline probe from the current base commit'
marker='  ptime t_e = t0 + milliseconds(pos_end);'
assert base.count(marker)==1
base=base.replace(marker,marker+'''
  auto log_seed=[&](int c,double x,double y) {
    std::ostringstream s;s<<std::setprecision(17);
    s<<"#PeakSeed {\\"id\\":"<<convergence_audit_id+1<<",\\"component\\":"<<c
      <<",\\"start\\":\\""<<to_iso_string(t_s)<<"\\",\\"end\\":\\""<<to_iso_string(t_e)
      <<"\\",\\"px\\":"<<x<<",\\"py\\":"<<y<<"}";
    std::cerr<<s.str()<<'\\n';
  };
''')
for needle,addition in [('    prm.p = sqrt(pxV[i] * pxV[i] + pyV[i] * pyV[i]);','    log_seed(2,pxV[i],pyV[i]);\n'),('      prm.p = sqrt(px[i] * px[i] + py[i] * py[i]);','      log_seed(icmp,px[i],py[i]);\n')]:
 assert base.count(needle)==1; base=base.replace(needle,addition+needle)
entries=json.loads((build/'compile_commands.json').read_text());entry=next(e for e in entries if e['file'].endswith('/cal_ccf.cpp'))
cmd0=shlex.split(entry['command']);obj0=Path(cmd0[cmd0.index('-o')+1]);link0=shlex.split((Path(entry['directory'])/obj0.parent/'link.txt').read_text())
for variant in ['baseline','fixed']:
 out=build/'peak-search-probe'/variant;out.mkdir(parents=True,exist_ok=True)
 text=base if variant=='baseline' else base[:search(base)]+current[search(current):]
 source=out/'cal_ccf.cpp';source.write_text(text);obj=out/'cal_ccf.cpp.o';exe=out/'cal_ccf_peak_search'
 cmd=cmd0.copy();cmd[cmd.index('-o')+1]=str(obj);cmd[cmd.index('-c')+1]=str(source);cmd+=['-I'+str(repo/'src')]
 subprocess.run(cmd,cwd=entry['directory'],check=True)
 objects={'cal_ccf.cpp.o':str(obj),'station_info.cpp.o':str(build/'io-fixed-current/station_info.cpp.o')}
 link=[objects.get(Path(s).name,s) for s in link0];link[link.index('-o')+1]=str(exe);subprocess.run(link,cwd=entry['directory'],check=True)
 (out/'build.json').write_text(json.dumps({'base_commit':baseline_commit,'variant':variant,'compile':cmd,'link':link,'sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [production_source,source,exe,Path(__file__),build/'io-fixed-current/station_info.cpp.o']}},indent=2)+'\n')
 print(exe,flush=True)
