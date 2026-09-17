"""Build signed/absolute convergence probes from the existing matrix audit.

Only the second convergence predicate changes; line-search behavior and all
other algorithms remain fixed. Production files are not edited.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('build',type=Path)
a=p.parse_args();build=a.build.resolve();repo=Path(__file__).resolve().parents[1]
source=(build/'matrix-audit/cal_ccf.cpp').read_text()
assert '1837u + 104729u * i' in source
start=source.rindex('static int est_dist_grad(');body=source.index('{',start)
depth=1;end=body+1
while depth:
 depth+=(source[end]=='{')-(source[end]=='}');end+=1
function=source[start:end]
function=function.replace('  cpu_profile::Timer cpu_timer(cpu_profile::fitting);','''  const int audit_id=++convergence_audit_id;
  auto audit_finish=[&](int result,const char *reason) {
    std::ostringstream s; s<<std::setprecision(17);
    s<<"#ConvergenceResult {\\"id\\":"<<audit_id<<",\\"result\\":"<<result
      <<",\\"reason\\":\\""<<reason<<"\\",\\"p\\":"<<prm.p<<",\\"theta\\":"<<prm.θ
      <<",\\"delta\\":"<<prm.Δ<<",\\"dp_delta\\":"<<prm.dp_Δ<<"}";
    std::cerr<<s.str()<<'\\n';return result;
  };
  { std::ostringstream s;s<<std::setprecision(17);
    s<<"#ConvergenceStart {\\"id\\":"<<audit_id<<",\\"p\\":"<<prm.p
      <<",\\"theta\\":"<<prm.θ<<",\\"delta\\":"<<prm.Δ<<",\\"dp_delta\\":"<<prm.dp_Δ<<"}";
    std::cerr<<s.str()<<'\\n'; }
  cpu_profile::Timer cpu_timer(cpu_profile::fitting);''')
function=function.replace('      bool flag_loop = 0;', '      bool flag_loop = 0;\n      double accepted_r=0.;')
function=function.replace('          flag_loop = 1;', '          flag_loop = 1;\n          accepted_r=r;')
old='abs(ε) < 1E-9 && dprm(2) * 180. / M_PI < 0.1'
assert function.count(old)==1
function=function.replace('      if ('+old+') {','''      { std::ostringstream s;s<<std::setprecision(17);
        s<<"#ConvergenceStep {\\"id\\":"<<audit_id<<",\\"iteration\\":"<<i+1
          <<",\\"epsilon\\":"<<ε<<",\\"raw_delta_degrees\\":"<<dprm(2)*180./M_PI
          <<",\\"accepted_delta_degrees\\":"<<-accepted_r*dprm(2)*180./M_PI
          <<",\\"r\\":"<<accepted_r<<",\\"old_stop\\":"<<(abs(ε)<1E-9 && dprm(2)*180./M_PI<.1)
          <<",\\"absolute_stop\\":"<<(abs(ε)<1E-9 && std::abs(dprm(2))*180./M_PI<.1)
          <<",\\"negative_eigenvalues\\":"<<num_eig<<"}";
        std::cerr<<s.str()<<'\\n'; }
      if ('''+old+''') {''')
function=function.replace('return (-1);','return audit_finish(-1,"line_search_failed");',1)
function=function.replace('return (-1);','return audit_finish(-1,num_loop==0 ? "iteration_limit" : "hessian_not_negative");',1)
function=function.replace('return (num_loop);','return audit_finish(num_loop,"converged");')
function=function.replace('return (-1);','return audit_finish(-1,"invalid_initial_distance");')
assert 'return (-1);' not in function
source=source[:start]+function+source[end:]
source=source.replace('static int est_dist_grad(', 'static int convergence_audit_id=0;\nstatic int est_dist_grad(',1)
# Record which optimizer call produced an emitted row (original matrix recursion is excluded).
marker='  if (!matrix_audit_writing) {'
assert source.count(marker)==1
source=source.replace(marker,marker+'''\n    std::cerr<<"#ConvergenceEvent {\\"id\\":"<<convergence_audit_id
      <<",\\"component\\":"<<icmp<<",\\"start\\":\\""<<to_iso_string(t_s)
      <<"\\",\\"end\\":\\""<<to_iso_string(t_e)<<"\\"}\\n";''')
entries=json.loads((build/'compile_commands.json').read_text())
entry=next(e for e in entries if e['file'].endswith('/cal_ccf.cpp'))
basecmd=shlex.split(entry['command']);original=Path(basecmd[basecmd.index('-o')+1])
linkbase=shlex.split((Path(entry['directory'])/original.parent/'link.txt').read_text())
root=build/'convergence-probe';root.mkdir(exist_ok=True)
for variant in ['signed','absolute']:
 out=root/variant;out.mkdir(exist_ok=True)
 text=source if variant=='signed' else source.replace(old,'abs(ε) < 1E-9 && std::abs(dprm(2)) * 180. / M_PI < 0.1')
 generated=out/'cal_ccf.cpp';generated.write_text(text)
 obj=out/'cal_ccf.cpp.o';cmd=basecmd.copy();cmd[cmd.index('-o')+1]=str(obj);cmd[cmd.index('-c')+1]=str(generated);cmd+=['-I'+str(repo/'src')]
 subprocess.run(cmd,cwd=entry['directory'],check=True)
 objects={'cal_ccf.cpp.o':str(obj),'station_info.cpp.o':str(build/'io-fixed-current/station_info.cpp.o')}
 link=[objects.get(Path(s).name,s) for s in linkbase];exe=out/'cal_ccf_convergence';link[link.index('-o')+1]=str(exe)
 subprocess.run(link,cwd=entry['directory'],check=True)
 (out/'build.json').write_text(json.dumps({'variant':variant,'compile':cmd,'link':link,'sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [generated,exe,build/'matrix-audit/cal_ccf.cpp',build/'io-fixed-current/station_info.cpp.o',Path(__file__)]}},indent=2)+'\n')
 print(exe,flush=True)
