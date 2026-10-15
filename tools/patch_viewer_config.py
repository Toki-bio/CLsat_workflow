#!/usr/bin/env python3
"""Turn the CLsat viewer page (v2.25) into a config-driven page: every Darevskia-specific
constant (family names, monomer lengths, colours, consensus sequences, sub-subfamily
consensuses, species order and labels, the three inspection lists, titles and notes) is read
from viewer_config.json at load time instead of being written into index.html.

Usage: patch_viewer_config.py <index.html in> <index.html out>
Every replacement must match exactly once; the script fails otherwise, so a changed upstream
page is noticed instead of silently half-patched. Used once to produce
satpipe/viewer_html/index.html from the viewer repository's audit-fixes page (109b38c).
"""
import re
import sys

src_path, out_path = sys.argv[1], sys.argv[2]
s = open(src_path, encoding='utf-8').read()
n_edits = 0


def rep(old, new, count=1):
    global s, n_edits
    k = s.count(old)
    if k != count:
        sys.exit('pattern found %d times (expected %d):\n%s' % (k, count, old[:300]))
    s = s.replace(old, new)
    n_edits += 1


def rep_re(pattern, new, flags=re.S):
    global s, n_edits
    m = list(re.finditer(pattern, s, flags))
    if len(m) != 1:
        sys.exit('regex found %d times (expected 1): %s' % (len(m), pattern[:200]))
    s = s[:m[0].start()] + new + s[m[0].end():]
    n_edits += 1


# ---- head: title ------------------------------------------------------------------------
rep('<title>CLsat Locus Detail Viewer – Per-window Bitscores</title>',
    '<title>Satellite locus viewer</title>')
rep('<h1>CLsat Locus Detail Viewer <span id="versionTag"',
    '<h1><span id="pageTitle">Satellite locus viewer</span> <span id="versionTag"')

# ---- filters: family / sub-subfamily options generated from the config ----------------------
rep_re(r'<select id="subtypeFilter" multiple size="8" style="height:auto; min-width:120px;">\s*<option value="ALL" selected>ALL</option>.*?</select>',
       '<select id="subtypeFilter" multiple size="8" style="height:auto; min-width:120px;">\n        <option value="ALL" selected>ALL</option>\n      </select>')
rep('<label title="CLsat major family (CLsat1-4) or sub-subfamily (CLsat1.1-1.5, CLsat3.1-3.2) assigned to the locus. Selecting a major family includes all of its sub-subfamilies.">Subtype</label>',
    '<label title="Family or sub-subfamily assigned to the locus. Selecting a family includes all of its sub-subfamilies.">Subtype</label>')

# ---- legend generated from the config --------------------------------------------------------
rep_re(r'<h2>CLsat colors</h2>\s*<div class="legend">.*?</div>\s*</div>',
       '<h2>Family colors</h2>\n    <div class="legend" id="familyLegend"></div>')

# ---- manual intro -------------------------------------------------------------------------
rep('<p>A locus-resolved browser for the CLsat satellite arrays catalogued across six <i>Darevskia</i>\n  genome assemblies (<span class="k">mix, arm, val, unp, unm, nai</span>). Every array in the catalogue is a',
    '<p id="manualIntro">A locus-resolved browser for the satellite arrays catalogued across the assemblies of this run. Every array in the catalogue is a')

# ---- file names and constants ---------------------------------------------------------------
rep_re(r'// Per-window bitscore files \(one per CLsat\)\nconst FILE_BITS_CLSAT = \{.*?\n\};',
       '// Per-window bitscore files (one per family) -- filled from viewer_config.json\nlet FILE_BITS_CLSAT = {};')
rep('const SAT_BITS    = 50;', 'let SAT_BITS    = 50;      // from viewer_config.json (expansion min_bits)')
rep('const EDGE_TOL    = 50;', 'let EDGE_TOL    = 50;      // from viewer_config.json')
rep('const FLANK_VIEW  = 5000;', 'let FLANK_VIEW  = 5000;    // from viewer_config.json')

rep_re(r'/\* =+\n   CONSENSUS SEQUENCES \(embedded so no server is needed\)\n=+ \*/\nconst CONSENSUS_SEQS = \{.*?\n\};\n',
       '''/* ============================================================================
   RUN CONFIGURATION (viewer_config.json, written by `satpipe viewer`)
   Families, monomer lengths, colours, consensus and sub-subfamily sequences, species order and
   labels, inspection lists and notes all come from that file; nothing run-specific is embedded.
============================================================================ */
let VIEWER_CONFIG = null;
let FAMILIES = [];            // family names, in config order (= column order of MASTER.locus_bits.tsv)
let CONSENSUS_SEQS = {};      // family -> monomer consensus
''')
rep_re(r'// Cluster consensus sequences \(cl1–cl9\) and their CLsat subfamily mapping\nconst CLUSTER_SEQS = \{.*?\n\};\nconst CLUSTER_TO_CLSAT = \{.*?\n\};\n',
       '// Sub-subfamily consensus sequences (optional, config key "subfamilies") and their family\nlet CLUSTER_SEQS = {};\nlet CLUSTER_TO_CLSAT = {};\n')
rep_re(r'// Sub-subfamily consensus sequences in the pipeline monomer phase.*?const SUBSUBFAM_SEQS = \{.*?\n\};\n',
       'let SUBSUBFAM_SEQS = CLUSTER_SEQS;\n')
rep('''// Reverse: CLsat type → list of cluster names
const CLSAT_TO_CLUSTERS = {};
for (const [cl, ct] of Object.entries(CLUSTER_TO_CLSAT)) {
  if (!CLSAT_TO_CLUSTERS[ct]) CLSAT_TO_CLUSTERS[ct] = [];
  CLSAT_TO_CLUSTERS[ct].push(cl);
}
''', '''// Reverse: family -> list of sub-subfamily names
let CLSAT_TO_CLUSTERS = {};
let CLSAT_COLORS = {};
let CLSAT_STEP = {};
let WINDOW_STEP = 146;        // display grid for grouping bars of different families (max monomer length)
let ORTHO_LOCI = [], FLAGGED_LOCI = [], FULL_TABLE = [], ASM_ORDER = [], SPECIES_LABEL = {};
let LIST_NOTES = {ortho: "", flag: "", full: ""};
let LIST_COUNTS = {};

function esc_(t){ return String(t).replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c])); }

/** Apply viewer_config.json: fill every run-specific table and the config-driven DOM parts. */
function applyViewerConfig(cfg){
  VIEWER_CONFIG = cfg;
  FAMILIES = cfg.families.map(f => f.name);
  CONSENSUS_SEQS = {}; CLSAT_COLORS = {}; CLSAT_STEP = {}; FILE_BITS_CLSAT = {};
  for (const f of cfg.families){
    CONSENSUS_SEQS[f.name] = (f.consensus || "").toUpperCase();
    CLSAT_COLORS[f.name] = f.colour || "#9ca3af";
    CLSAT_STEP[f.name] = f.mono;
    FILE_BITS_CLSAT[f.name] = f.bits_file || ("MASTER_" + f.name + ".bits.tsv");
  }
  WINDOW_STEP = cfg.window_step || Math.max(...FAMILIES.map(f => CLSAT_STEP[f]), 1);
  if (cfg.sat_bits != null) SAT_BITS = +cfg.sat_bits;
  if (cfg.edge_tol != null) EDGE_TOL = +cfg.edge_tol;
  if (cfg.flank_view != null) FLANK_VIEW = +cfg.flank_view;
  CLUSTER_SEQS = {}; CLUSTER_TO_CLSAT = {}; CLSAT_TO_CLUSTERS = {};
  for (const [name, rec] of Object.entries(cfg.subfamilies || {})){
    CLUSTER_SEQS[name] = (rec.seq || "").toUpperCase();
    CLUSTER_TO_CLSAT[name] = rec.family;
    (CLSAT_TO_CLUSTERS[rec.family] = CLSAT_TO_CLUSTERS[rec.family] || []).push(name);
  }
  SUBSUBFAM_SEQS = CLUSTER_SEQS;
  ASM_ORDER = (cfg.species || []).map(x => x.code);
  SPECIES_LABEL = {}; for (const x of (cfg.species || [])) SPECIES_LABEL[x.code] = x.label || x.code;
  const L = cfg.lists || {};
  ORTHO_LOCI = L.ortho || []; FLAGGED_LOCI = L.flagged || []; FULL_TABLE = L.full || [];
  FLAGGED_LOCI.sort((a,b)=>(a.cat==="UNIQUE_SEED_CONTRADICTION")===(b.cat==="UNIQUE_SEED_CONTRADICTION")?0:(a.cat==="UNIQUE_SEED_CONTRADICTION"?-1:1));
  LIST_NOTES = Object.assign({ortho:"", flag:"", full:""}, cfg.notes || {});
  LIST_COUNTS = cfg.counts || {};
  // DOM: title, legend, family / sub-subfamily filter
  const title = cfg.title || "Satellite locus viewer";
  document.title = title;
  const pt = document.getElementById("pageTitle"); if (pt) pt.textContent = title;
  const leg = document.getElementById("familyLegend");
  if (leg) leg.innerHTML = FAMILIES.map(f => '<div class="legendItem"><span class="swatch" style="background:' + CLSAT_COLORS[f] + '"></span> ' + esc_(f) + ' <span style="color:#64748b">' + CLSAT_STEP[f] + ' bp</span></div>').join("");
  const sel = document.getElementById("subtypeFilter");
  if (sel){
    sel.innerHTML = '<option value="ALL" selected>ALL</option>' + FAMILIES.map(f =>
      '<option value="' + esc_(f) + '">' + esc_(f) + '</option>' +
      (CLSAT_TO_CLUSTERS[f] || []).filter(sname => sname !== f).map(sname => '<option value="' + esc_(sname) + '">' + esc_(sname) + '</option>').join("")).join("");
  }
  const mi = document.getElementById("manualIntro");
  if (mi && ASM_ORDER.length) mi.innerHTML = 'A locus-resolved browser for the satellite arrays catalogued across ' + ASM_ORDER.length +
    ' assemblies (<span class="k">' + ASM_ORDER.map(esc_).join(", ") + '</span>) in run <span class="k">' + esc_(cfg.run_name || "") + '</span>. Every array in the catalogue is a ' +
    '<b>segmented MASTER locus</b> &mdash; a contiguous satellite-similar domain, split where the dominant family changes along the array.';
  console.log("[viewer] config applied: " + FAMILIES.length + " families, " + ASM_ORDER.length + " assemblies, lists " +
              ORTHO_LOCI.length + "/" + FLAGGED_LOCI.length + "/" + FULL_TABLE.length);
}
''')

# ---- colours / step / grid blocks removed (now set by applyViewerConfig) --------------------
rep_re(r'const CLSAT_COLORS = \{\n  CLsat1: "#60a5fa",.*?\n\};\n', '// CLSAT_COLORS: see applyViewerConfig()\n')
rep_re(r'// Per-CLsat monomer lengths \(= scan window step for each consensus\)\nconst CLSAT_STEP = \{.*?\n\};\n',
       '// CLSAT_STEP (monomer length per family): see applyViewerConfig()\n')
rep('const WINDOW_STEP = 146;  // uniform display grid for grouping bars from different CLsats',
    '// WINDOW_STEP: see applyViewerConfig()')

# ---- family loops -----------------------------------------------------------------------------
rep('      byBin.set(binPos, {pos: binPos, CLsat1: 0, CLsat2: 0, CLsat3: 0, CLsat4: 0});',
    '      const b0 = {pos: binPos}; for (const f of FAMILIES) b0[f] = 0;\n      byBin.set(binPos, b0);')
rep('    w.best  = Math.max(w.CLsat1, w.CLsat2, w.CLsat3, w.CLsat4);',
    '    w.best  = Math.max(0, ...FAMILIES.map(f => w[f] || 0));')
rep('      byBin.set(binPos, {pos: binPos, CLsat1: 0, CLsat2: 0, CLsat3: 0, CLsat4: 0,\n                         hasMain: false, hasFlank: false});',
    '      const b0 = {pos: binPos, hasMain: false, hasFlank: false}; for (const f of FAMILIES) b0[f] = 0;\n      byBin.set(binPos, b0);')
rep('    w.best      = Math.max(w.CLsat1, w.CLsat2, w.CLsat3, w.CLsat4);',
    '    w.best      = Math.max(0, ...FAMILIES.map(f => w[f] || 0));')
rep('  const CLSATS = ["CLsat1", "CLsat2", "CLsat3", "CLsat4"];', '  const CLSATS = FAMILIES;', count=2)
rep('  s += `Max bits: CLsat1=${fmt(b.b1||0)} CLsat2=${fmt(b.b2||0)} CLsat3=${fmt(b.b3||0)} CLsat4=${fmt(b.b4||0)}`;',
    '  s += "Max bits: " + FAMILIES.map(f => f + "=" + fmt(b[f] || 0)).join(" ");')
rep('    const consensusName = assignment?.seed || "CLsat1";', '    const consensusName = assignment?.seed || FAMILIES[0];')
rep("    const mono = (typeof CLSAT_STEP !== 'undefined' && CLSAT_STEP[a.seed]) ? CLSAT_STEP[a.seed] : 146;",
    "    const mono = CLSAT_STEP[a.seed] || WINDOW_STEP;")
rep("    consensus_length: CONSENSUS_SEQS[usedClsats[0]]?.length || 146,",
    "    consensus_length: CONSENSUS_SEQS[usedClsats[0]]?.length || WINDOW_STEP,")
# bits matrix keyed by family name, one column per configured family
rep('''    obj[p[0]] = { b1: +p[1], b2: +p[2], b3: +p[3], b4: +p[4] };''',
    '''    const rec = {};
    FAMILIES.forEach((f, i) => { rec[f] = +p[i + 1]; });
    obj[p[0]] = rec;''')
rep('    if (p.length < 5) continue;\n    if (p[0].toLowerCase() === "locus_id" || p[0].toLowerCase() === "locus") continue;',
    '    if (p.length < 1 + FAMILIES.length) continue;\n    if (p[0].toLowerCase() === "locus_id" || p[0].toLowerCase() === "locus") continue;')

# ---- version tag: show the generator too -----------------------------------------------------
rep('const VIEWER_VERSION = "v2.25-audit-20261015";', 'const VIEWER_VERSION = "v3.0-config-20261015";')

# ---- loading: fetch the config first ---------------------------------------------------------
rep('''async function loadViaFetch(){
  setStatus("Loading files via fetch()…");
''', '''async function loadViaFetch(){
  setStatus("Loading viewer_config.json…");
  try {
    const cfgText = await fetchText("viewer_config.json");
    applyViewerConfig(JSON.parse(cfgText));
  } catch (e) {
    setStatus("viewer_config.json is missing or invalid: " + e.message +
              "\\nThis page is run-specific only through that file; write it with `satpipe viewer <run_dir>`.");
    fileOk.textContent = "\\u2717 no viewer_config.json";
    fileOk.className = "pill bad";
    throw e;
  }
  setStatus("Loading files via fetch()…");
''')

# ---- inspection lists: constants removed, presets/notes from the config ------------------------
rep_re(r'const ORTHO_LOCI = \[.*?\];\nconst FLAGGED_LOCI = \[.*?\];\nFLAGGED_LOCI\.sort\(.*?\);\nconst FULL_TABLE = \[.*?\];\nconst ASM_ORDER = \[.*?\];\n',
       '// ORTHO_LOCI, FLAGGED_LOCI, FULL_TABLE, ASM_ORDER: filled from viewer_config.json by applyViewerConfig()\n')
rep('''      ["all","All flagged (28 + 7 candidates)", d=>true],''',
    '''      ["all","All flagged", d=>true],''')
rep('''      ["cand","Uncatalogued candidates (C1–C7)", d=>d.cat==="CANDIDATE_UNCATALOGUED"],''',
    '''      ["cand","Uncatalogued candidates", d=>d.cat==="CANDIDATE_UNCATALOGUED"],''')
rep('''      ["all","All loci (627)", d=>!String(d.id).startsWith("C")],''',
    '''      ["all","All loci", d=>!String(d.id).startsWith("C")],''')
rep_re(r'  const NOTES=\{\n    ortho:".*?",\n    flag:".*?",\n    full:".*?"\n  \};\n',
       '  // list notes come from viewer_config.json (LIST_NOTES), generated from the run\'s own numbers\n')
rep('    noteEl.textContent=NOTES[m];', '    noteEl.textContent=LIST_NOTES[m]||"";')
rep('''    title.textContent = m==="ortho" ? "Orthologous loci — positional conservation spectrum"
                       : m==="flag" ? "Flagged loci — arbitration margins"
                                    : "Flank orthology table — every locus, all six assemblies";''',
    '''    title.textContent = m==="ortho" ? "Orthologous loci — positional conservation spectrum"
                       : m==="flag" ? "Flagged loci — arbitration margins"
                                    : "Flank orthology table — every locus, all assemblies";''')
rep('''      const b=document.createElement("button");
      b.className="preset"+(k===preset?" active":""); b.textContent=lab;''',
    '''      const b=document.createElement("button");
      const nAll = (m==="ortho" ? ORTHO_LOCI : m==="flag" ? FLAGGED_LOCI : FULL_TABLE.filter(d=>!String(d.id).startsWith("C"))).length;
      b.className="preset"+(k===preset?" active":""); b.textContent = (k==="all") ? lab+" ("+nAll+")" : lab;''')

open(out_path, 'w', encoding='utf-8', newline='\n').write(s)
print('patched: %d edits -> %s' % (n_edits, out_path))
