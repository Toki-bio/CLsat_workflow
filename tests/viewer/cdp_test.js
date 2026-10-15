// Functional test of the config-driven viewer page in headless Chrome via CDP.
// Usage: node cdp_test.js <dir> <http port>
const { exec, execFile, spawn } = require('child_process');
const os = require('os'); const path = require('path'); const fs = require('fs');
const DIR = process.argv[2]; const PORT = +process.argv[3] || 8765;
const CHROME = 'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe';
const chromePort = 9224 + Math.floor(Math.random() * 1000);
const profileDir = path.join(os.tmpdir(), 'cdp-test-profile-' + Date.now());
fs.mkdirSync(profileDir, { recursive: true });
const sleep = ms => new Promise(r => setTimeout(r, ms));

function killThisInstanceOnly() {
  return new Promise((resolve) => {
    execFile('wmic', ['process', 'where', `(name='chrome.exe') and (CommandLine like '%${profileDir.replace(/\\/g, '\\\\')}%')`, 'get', 'ProcessId'], (err, stdout) => {
      if (err || !stdout) return resolve();
      const pids = stdout.split('\n').map(s => s.trim()).filter(s => /^\d+$/.test(s));
      if (!pids.length) return resolve();
      let remaining = pids.length;
      pids.forEach(pid => execFile('taskkill', ['/PID', pid, '/T', '/F'], () => { if (--remaining === 0) resolve(); }));
    });
  });
}

(async () => {
  // DIR may be a directory (served locally) or an http(s) URL of a deployed page
  const remote = /^https?:\/\//.test(DIR);
  const server = remote ? { kill() {} } : spawn('python', ['-m', 'http.server', String(PORT), '--bind', '127.0.0.1'], { cwd: DIR, stdio: 'ignore' });
  await sleep(1500);
  const url = remote ? DIR : `http://127.0.0.1:${PORT}/index.html`;
  exec(`"${CHROME}" --headless=new --disable-gpu --no-sandbox --user-data-dir="${profileDir}" --remote-debugging-port=${chromePort} --window-size=1500,1000 "${url}"`, { timeout: 120000 });
  await sleep(3000);
  let targets = [];
  for (let i = 0; i < 10 && !targets.length; i++) {
    try { targets = (await (await fetch(`http://127.0.0.1:${chromePort}/json`)).json()).filter(t => t.type === 'page'); } catch (e) { await sleep(500); }
  }
  if (!targets.length) { console.log('no chrome target'); server.kill(); await killThisInstanceOnly(); process.exit(2); }
  const ws = new WebSocket(targets[0].webSocketDebuggerUrl);
  await new Promise(r => ws.onopen = r);
  let id = 0; const pending = {}; const console_msgs = []; const errors = [];
  ws.onmessage = ev => {
    const m = JSON.parse(ev.data);
    if (m.id && pending[m.id]) { pending[m.id](m); delete pending[m.id]; }
    if (m.method === 'Runtime.consoleAPICalled') console_msgs.push(m.params.type + ': ' + m.params.args.map(a => a.value ?? a.description ?? '').join(' '));
    if (m.method === 'Runtime.exceptionThrown') errors.push(m.params.exceptionDetails.text + ' ' + (m.params.exceptionDetails.exception?.description || ''));
    if (m.method === 'Log.entryAdded' && m.params.entry.level === 'error') errors.push('LOG ' + m.params.entry.text + ' ' + (m.params.entry.url || ''));
  };
  const send = (method, params = {}) => new Promise(res => { const i = ++id; pending[i] = res; ws.send(JSON.stringify({ id: i, method, params })); });
  const evalJs = async (expr) => { const r = await send('Runtime.evaluate', { expression: expr, returnByValue: true, awaitPromise: true }); return r.result?.result?.value ?? r.result?.exceptionDetails?.text; };
  await send('Runtime.enable'); await send('Page.enable'); await send('Log.enable');
  await send('Page.navigate', { url });
  // wait for the data load
  let status = '';
  for (let i = 0; i < 60; i++) {
    await sleep(1000);
    status = await evalJs('document.getElementById("status")?.textContent || document.querySelector("#statusEl,#status,pre")?.textContent || ""');
    if (/Select a scaffold|Per-window data|missing|invalid|Could not/.test(status || '')) break;
  }
  const info = await evalJs(`JSON.stringify({
    title: document.title, h1: document.getElementById('pageTitle')?.textContent,
    families: typeof FAMILIES !== 'undefined' ? FAMILIES : null,
    steps: typeof CLSAT_STEP !== 'undefined' ? CLSAT_STEP : null,
    colours: typeof CLSAT_COLORS !== 'undefined' ? CLSAT_COLORS : null,
    asm: typeof ASM_ORDER !== 'undefined' ? ASM_ORDER : null,
    lists: typeof ORTHO_LOCI !== 'undefined' ? [ORTHO_LOCI.length, FLAGGED_LOCI.length, FULL_TABLE.length] : null,
    loci: typeof loci !== 'undefined' ? loci.length : null,
    assignments: typeof assignments !== 'undefined' ? Object.keys(assignments).length : null,
    bitsSample: typeof bits !== 'undefined' ? bits['L1'] : null,
    legend: document.getElementById('familyLegend')?.children.length,
    filterOptions: [...document.querySelectorAll('#subtypeFilter option')].map(o => o.value),
    subfam: typeof CLSAT_TO_CLUSTERS !== 'undefined' ? CLSAT_TO_CLUSTERS : null,
    speciesButtons: [...document.querySelectorAll('.spBtn')].map(b => b.textContent),
    lociTableRows: document.querySelectorAll('#locusTable tbody tr, table tbody tr').length
  })`);
  console.log('STATUS:', (status || '').split('\n').slice(0, 4).join(' | '));
  console.log('INFO:', info);
  // exercise: select species, click first locus row, open the three lists
  const steps = [
    `(() => { const b=[...document.querySelectorAll('.spBtn')].find(x=>x.textContent==='arm'); b && b.click(); return !!b; })()`,
    `(() => { const r=document.querySelector('tr.lrow, tbody tr'); r && r.click(); return r ? r.textContent.slice(0,40) : null; })()`,
    `(() => { document.getElementById('orthoListBtn').click(); const n=document.querySelectorAll('#listBody tbody tr').length; const t=document.getElementById('listTitle').textContent; const p=[...document.querySelectorAll('#listPresets button')].map(b=>b.textContent); document.getElementById('listClose').click(); return {rows:n, title:t, presets:p}; })()`,
    `(() => { document.getElementById('flagListBtn').click(); const n=document.querySelectorAll('#listBody tbody tr').length; const first=document.querySelector('#listBody tbody tr')?.textContent.slice(0,60); document.getElementById('listClose').click(); return {rows:n, first}; })()`,
    `(() => { document.getElementById('fullTableBtn').click(); const n=document.querySelectorAll('#listBody tbody tr').length; const note=document.getElementById('listNote').textContent.slice(0,120); document.getElementById('listClose').click(); return {rows:n, note}; })()`,
    `(() => { const s=document.getElementById('subtypeFilter'); [...s.options].forEach(o=>o.selected=(o.value==='CLsat3.1')); s.dispatchEvent(new Event('change')); return document.querySelectorAll('tr.lrow, #locusTable tbody tr').length; })()`,
    `document.getElementById('status')?.textContent?.slice(0,200) || (document.querySelector('pre')||{}).textContent?.slice(0,200)`
  ];
  for (const st of steps) { await sleep(600); console.log('STEP:', JSON.stringify(await evalJs(`JSON.stringify(${st})`))); }
  await sleep(1000);
  const shot = await send('Page.captureScreenshot', { format: 'png' });
  fs.writeFileSync(path.join(remote ? os.tmpdir() : DIR, 'screenshot.png'), Buffer.from(shot.result.data, 'base64'));
  console.log('CONSOLE (' + console_msgs.length + '):'); console_msgs.slice(0, 25).forEach(m => console.log('  ' + m.slice(0, 200)));
  console.log('ERRORS (' + errors.length + '):'); errors.forEach(e => console.log('  ' + e.slice(0, 300)));
  ws.close(); server.kill(); await killThisInstanceOnly();
  process.exit(errors.length ? 1 : 0);
})().catch(async e => { console.error('TEST FAILURE', e); await killThisInstanceOnly(); process.exit(2); });
