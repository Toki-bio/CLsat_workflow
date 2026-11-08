"""WP9 test oracle: run Python's stage 21_anchor (satpipe/run.py st_anchor and its plan) on a config given as JSON
(the validated config data) with a stub Ctx, for tests/sh/test_anchor.sh.

  python tests/port/anchor_stage_oracle.py CONFIG.json RUN_DIR
stdout: 'PLAN_HASH\t<hash>', 'INPUT\t<path>'..., then 'SAY\t<text>', 'MANIFEST\t<json sort_keys>', 'OUTPUT\t<path>'...
or 'ERROR\t<class>\t<message with \\n escaped>'. Files land in RUN_DIR/anchored as the stage writes them.
"""
import hashlib
import json
import os
import sys
import types

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from satpipe.run import st_anchor  # noqa: E402


def main():
    with open(sys.argv[1], encoding='utf-8') as fh:
        data = json.load(fh)
    run = sys.argv[2]
    out = []
    ctx = types.SimpleNamespace(cfg=types.SimpleNamespace(data=data, families=data['families']), manifest={}, run=run)

    def d(*parts):
        p = os.path.join(run, *parts)
        os.makedirs(p, exist_ok=True)
        return p
    ctx.d = d
    ctx.say = lambda m: out.append('SAY\t' + m)
    inputs, params = st_anchor.plan(ctx)
    print('PLAN_HASH\t' + hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()[:16])
    for p in inputs:
        print('INPUT\t' + p)
    try:
        outputs = st_anchor(ctx)
    except Exception as e:  # noqa: BLE001
        msg = str(e).replace(os.sep, '/')        # Windows: os.path.join writes '\'; the test compares '/' paths
        print('ERROR\t%s\t%s' % (type(e).__name__, msg.replace('\\', '\\\\').replace('\n', '\\n')))
        return
    for line in out:
        print(line)
    print('MANIFEST\t' + json.dumps(ctx.manifest['anchor'], sort_keys=True))
    for p in outputs:
        print('OUTPUT\t' + p.replace(os.sep, '/'))


if __name__ == '__main__':
    sys.stdout.reconfigure(newline='\n')
    main()
