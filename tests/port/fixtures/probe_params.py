import sys, os, json, hashlib
sys.path.insert(0, r'C:\work\satpipe_bashport')
from satpipe.config import load_config
from satpipe import run as R
def flat(node, path, out):
    if isinstance(node, dict):
        if not node: out.append((path, 'emptymap', '{}'))
        for k, v in node.items(): flat(v, f'{path}/{k}' if path else k, out)
    elif isinstance(node, list):
        if not node: out.append((path, 'emptylist', '[]'))
        for i, v in enumerate(node): flat(v, f'{path}/{i}', out)
    else:
        out.append((path, type(node).__name__, str(node)))
for cfgp in sys.argv[1:]:
    cfg = load_config(cfgp)
    print('=====', cfgp)
    rows = []; flat(cfg.data, '', rows)
    for p, t, v in rows: print('FLAT', p, t, repr(v), sep='\t')
    ctx = R.Ctx(cfg=cfg, run=cfg.data['out_dir'])
    for fn in R.STAGES:
        try:
            ins, params = fn.plan(ctx) if hasattr(fn, 'plan') else ([], {})
        except Exception as e:
            print('PLAN', fn.sid, 'ERROR', repr(e), sep='\t'); continue
        js = json.dumps(params, sort_keys=True)
        print('PLAN', fn.sid, hashlib.sha256(js.encode()).hexdigest()[:16], js, sep='\t')
        print('INPUTS', fn.sid, json.dumps(ins), sep='\t')
