#!/usr/bin/env python3
"""Fuzz documents for the yaml.awk differential (PORT_SPEC.md section 6, WP2).

    python3 tests/port/gen_yaml_fuzz.py --seed 1 --n 200 OUTDIR     # N documents inside the accepted subset
    python3 tests/port/gen_yaml_fuzz.py --refused OUTDIR              # one file per refused feature

Subset documents use: nested block mappings with a random but consistent indent step per mapping, one-line flow
mappings/sequences nested to depth 3 (random spacing, trailing commas, empties), plain scalars of every kept type and
string forms that look numeric to a human but are strings for PyYAML (1e3, 09, -.5), bool/null words, double-quoted
scalars with the five allowed escapes and literal tabs, single-quoted with '', typed and quoted keys (incl. '/', ': ',
'#', non-ASCII), duplicate keys, full-line and trailing comments, blank lines, trailing spaces, CRLF, BOM, '---'.
Refused files carry their expected reason in the name; each must make yaml.awk exit 1 with a 'YAML subset:' line."""
import argparse
import os
import random
import re

BRK = "\x00"          # possible line break inside a flow collection, resolved by breaks()
INDICATORS ="-?:,[]{}#&*!|>'\"%@`"
REFUSED_RE = [r"[-+]?0b[01_]+", r"[-+]?0x[0-9a-fA-F_]+", r"[-+]?0[0-7_]+", r"[-+]?(0|[1-9][0-9_]*)",
              r"[-+]?[1-9][0-9_]*(:[0-5]?[0-9])+", r"[-+]?[0-9][0-9_]*(:[0-5]?[0-9])+\.[0-9_]*",
              r"[-+]?[0-9][0-9_]*\.[0-9_]*([eE][-+][0-9]+)?", r"\.[0-9][0-9_]*([eE][-+][0-9]+)?",
              r"[-+]?\.(inf|Inf|INF)", r"\.(nan|NaN|NAN)", r"[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}.*", r"=", r"<<",
              r"~|null|Null|NULL|true|True|TRUE|false|False|FALSE|yes|Yes|YES|no|No|NO|on|On|ON|off|Off|OFF"]
WORDS = ["CLsat1", "spA", "SatX.1", "x-y", "node_7", "Sat X", "ref/SatX.cons.fa", "../ref/x.fa", "/staging/a b/c.fna",
         "x,y", "f(x)", "a#b", "a:b", "http://x.org/a?b=1", "é", "中文", "über alles", "plus", "both", "drop",
         "2 5 7 80 10 50 2000", "#1f77b4", "P. formosa (GCA_013036135.2)", "a;b", "50%", "~user/x", "C:/x/y",
         "-a", "-1x", "?x", ":x", "a b c", "it's", 'say "hi"', "a\\b", "1e3", "1.0e3", "09", "-.5", "+.5", "1.e5",
         "0o10", "1,000", "12abc", "NaN", "Infinity", "y", "n", "tRue", "0.1.2", "1-2", "+-1", "--x", "a.-"]
INTS = ["0", "-0", "+0", "7", "-12", "+5", "146", "100000000", "123456789012345678901234567890", "-99999999999999999999"]


def floats(rng):
    d = lambda k: "".join(rng.choice("0123456789") for _ in range(k))
    mant = rng.choice([f"{rng.randint(0, 999)}.{d(rng.randint(0, 6))}", f"{rng.randint(0, 99)}.", f".{d(rng.randint(1, 5))}",
                       "00.5", f"{d(rng.randint(1, 9))}.{d(rng.randint(1, 9))}", "0.1", "0.90", "80.0", "1.5"])
    if rng.random() < 0.3:
        mant += rng.choice("eE") + rng.choice("+-") + str(rng.randint(0, 320))
    if rng.random() < 0.3 and not mant.startswith("."):
        mant = rng.choice("+-") + mant
    if float(mant) in (float("inf"), float("-inf")):        # overflow is refused, tested apart
        mant = mant.split("e")[0].split("E")[0]
    return mant


def plain_ok(s, flow):
    if not s or s != s.strip() or "\t" in s or ": " in s or s.endswith(":") or " #" in s:
        return False
    c, nx = s[0], s[1:2]
    if c in INDICATORS and not (nx not in ("", " ") and (c == "-" or (not flow and c in "?:"))):
        return False
    if s.startswith("---") or s.startswith("..."):
        return False
    if flow and (any(ch in s for ch in ",?[]{}") or re.search(r":[,\[\]{}]", s)):
        return False
    return not any(re.fullmatch(r, s) for r in REFUSED_RE)


def rand_plain(rng, flow):
    for _ in range(100):
        if rng.random() < 0.6:
            s = rng.choice(WORDS)
        else:
            s = "".join(rng.choice("abcXYZ019._-/() ,#:é'\"") for _ in range(rng.randint(1, 14)))
        if plain_ok(s, flow):
            return s
    return "fallback"


def dq(rng):
    parts = []
    for _ in range(rng.randint(0, 6)):
        parts.append(rng.choice(["a", "B c", "\\\"", "\\\\", "\\n", "\\t", "\\u00e9", "\\u4e2d", "\\u0041", "\t", "#",
                                 ": ", " # x", "'", ",[]{}", "é", "/", "  "]))
    return '"' + "".join(parts) + '"'


def sq(rng):
    parts = [rng.choice(["a", "''", "\\n", "#", ": ", '"', "x y", "é", "\t", "[1]"]) for _ in range(rng.randint(0, 5))]
    return "'" + "".join(parts) + "'"


def scalar(rng, flow):
    r = rng.random()
    if r < 0.15:
        return rng.choice(INTS)
    if r < 0.32:
        return floats(rng)
    if r < 0.40:
        return rng.choice(["true", "True", "TRUE", "false", "yes", "No", "on", "OFF", "~", "null", "Null", "NULL"])
    if r < 0.55:
        return dq(rng)
    if r < 0.65:
        return sq(rng)
    return rand_plain(rng, flow)


def flow(rng, depth, force_map=False):
    sp = lambda: rng.choice(["", " ", "  "])
    if not force_map and rng.random() < 0.5:
        items = [flow(rng, depth + 1) if depth < 3 and rng.random() < 0.25 else scalar(rng, True)
                 for _ in range(rng.randint(0, 4))]
        tail = "," if items and rng.random() < 0.15 else ""
        return "[" + BRK + sp() + ("," + BRK + sp()).join(items) + tail + sp() + BRK + "]"
    ents = []
    for _ in range(rng.randint(0, 4)):
        k = key(rng, True)
        r = rng.random()
        if r < 0.1:
            ents.append(k)                                    # key without value -> null
        elif r < 0.15:
            ents.append(k + ":")                              # 'k:' then ',' or '}' -> null
        else:
            v = flow(rng, depth + 1) if depth < 3 and rng.random() < 0.25 else scalar(rng, True)
            ents.append(k + ":" + rng.choice([" ", " ", BRK + " "]) + v)
    tail = "," if ents and rng.random() < 0.15 else ""
    return "{" + BRK + sp() + ("," + BRK + sp()).join(ents) + tail + sp() + BRK + "}"


def key(rng, in_flow):
    r = rng.random()
    if r < 0.1:
        return rng.choice(["1", "0", "-3", "1.5", "2.0", "yes", "off", "null", "~", "True"])
    if r < 0.25:
        return rng.choice(['"a/b"', '"x: y"', '"#k"', "'it''s'", '"é/中"', '""', '"tab\\there"', '"1"', "'-'"])
    k = rand_plain(rng, in_flow)
    while k.startswith("-") or ":" in k:                      # keep plain keys simple to read
        k = rand_plain(rng, in_flow)
    return k


def comment(rng):
    return rng.choice(["# c", "#x: y", "# a # b", "#", "#\t tab in comment", "# é"])


def block(rng, indent, depth, lines):
    step = rng.randint(1, 4)
    keys = []
    for _ in range(rng.randint(1, 6)):
        if rng.random() < 0.15:
            lines.append(" " * rng.randint(0, 8) + comment(rng))
        if rng.random() < 0.1:
            lines.append(" " * rng.randint(0, 3))
        k = rng.choice(keys) if keys and rng.random() < 0.08 else key(rng, False)
        keys.append(k)
        pre = " " * indent + k + rng.choice([":", ":", " :"])
        r = rng.random()
        if depth < 4 and r < 0.25:
            lines.append(pre + rng.choice(["", "", "  ", " " + comment(rng)]))
            block(rng, indent + step, depth + 1, lines)
        elif r < 0.32:
            lines.append(pre + rng.choice(["", " ", " " + comment(rng)]))   # null
        else:
            v = flow(rng, 1) if rng.random() < 0.3 else scalar(rng, False)
            tail = rng.choice(["", "", "", "  ", " " + comment(rng), "   " + comment(rng)])
            if v[0] in "\"'[{" and rng.random() < 0.1:
                tail = comment(rng)                           # comment glued to a closing quote/bracket
            lines.append(pre + " " + v + tail)


def breaks(rng, line):
    """Resolve BRK markers: in 30 % of lines some markers become a line break (optionally after a comment, optionally
    followed by blank or comment lines) with a random continuation indent, as configs/darevskia_clsat.yaml does on the
    flank-single-copy branch; all other markers vanish."""
    multi = rng.random() < 0.3
    parts = line.split(BRK)
    out = [parts[0]]
    for part in parts[1:]:
        if multi and rng.random() < 0.3:
            out.append(rng.choice(["", "", " # c", "  #x"]) + "\n" + rng.choice(["", "\n", "   # mid\n"])
                       + " " * rng.randint(0, 10))
        out.append(part)
    return "".join(out)


def document(rng):
    lines = []
    if rng.random() < 0.1:
        lines.append(rng.choice(["---", "--- # doc"]))
    if rng.random() < 0.2:
        lines.append(comment(rng))
    if rng.random() < 0.05:
        lines.append(flow(rng, 1, force_map=True) + rng.choice(["", " " + comment(rng)]))     # top-level flow mapping
    else:
        block(rng, rng.choice([0, 0, 0, 2]), 1, lines)
    nl = "\r\n" if rng.random() < 0.2 else "\n"
    lines = [breaks(rng, ln) for ln in lines]
    text = nl.join(ln.replace("\n", nl) for ln in lines) +(nl if rng.random() < 0.9 else "")
    if rng.random() < 0.05:
        text = "\ufeff" + text
    return text


REFUSED = {
    "octal": "a: 010\n", "underscore_int": "a: 1_000\n", "underscore_float": "a: 1_0.5\n", "sexagesimal": "a: 1:30\n",
    "sexagesimal_float": "a: 1:30.5\n", "hex": "a: 0x1F\n", "binary": "a: 0b101\n", "inf": "a: .inf\n",
    "neg_inf": "a: -.inf\n", "nan": "a: .nan\n", "timestamp": "a: 2026-11-08\n", "timestamp_time": "a: 2026-11-08 10:00:00\n",
    "timestamp_key": "2026-11-08: x\n", "float_overflow": "a: 1.0e+400\n", "octal_in_flow": "a: [1, 010]\n",
    "underscore_key": "1_0: x\n",
    "anchor": "a: &x 1\nb: 2\n", "alias": "a: 1\nb: *x\n", "tag": "a: !!str 1\n", "merge_key": "b: 1\n<<: {a: 1}\n",
    "value_tag": "a: =\n", "block_seq": "a:\n  - x\n", "block_seq_same_indent": "a:\n- x\n", "literal": "a: |\n  x\n",
    "folded": "a: >\n  x\n", "multiline_plain": "a: foo\n  bar\n", "multiline_plain_nullkey": "a:\n  foo\n",
    "multiline_dquote": 'a: "foo\n  bar"\n', "multiline_squote": "a: 'foo\n  bar'\n", "dquote_line_continuation": 'a: "foo\\\n  bar"\n',
    "flow_plain_folded": "a: [x\n  y]\n", "flow_key_colon_next_line": "a: {b\n  : 1}\n", "flow_unclosed": "a: [1,\n",
    "flow_doc_marker": "a: [x,\n---\n]\n", "flow_tab_continuation": "a: [x,\n\t y]\n", "tab_indent": "a:\n\tb: 1\n",
    "tab_separator": "a:\t1\n", "tab_in_plain": "a: x\ty\n", "several_docs": "a: 1\n---\nb: 2\n", "doc_end": "a: 1\n...\n",
    "directive": "%YAML 1.1\n---\na: 1\n", "complex_key": "? a\n: b\n", "escape_x": 'a: "\\x41"\n', "escape_r": 'a: "\\r"\n',
    "escape_slash": 'a: "\\/"\n', "escape_surrogate": 'a: "\\ud800"\n', "flow_pair_in_seq": "a: [b: c]\n",
    "empty_key": ": x\n", "nested_map_one_line": "a: b: c\n", "bad_dedent": "a:\n    b: 1\n  c: 2\n",
    "unexpected_indent": "a: 1\n  b: 2\n", "flow_key": "a: {[x]: 1}\n", "empty_flow_entry": "a: [1, , 2]\n",
    "text_after_quote": "a: 'x' y\n", "flow_value_own_line": "a:\n  [1, 2]\n", "content_after_flow_root": "{a: 1}\nb: 2\n",
    "reserved_at": "a: @x\n",
}
REFUSED_BYTES = {
    "lone_cr": b"a: 1\rb: 2\n", "control_char": b"a: x\x01y\n", "invalid_utf8": b"a: \xff\n", "nel": b"a: x\xc2\x85y\n",
    "line_separator": b"a: x\xe2\x80\xa8y\n",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--refused", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.refused:
        for name, text in REFUSED.items():
            with open(os.path.join(a.out, f"{name}.yaml"), "wb") as fh:
                fh.write(text.encode("utf-8"))
        for name, data in REFUSED_BYTES.items():
            with open(os.path.join(a.out, f"{name}.yaml"), "wb") as fh:
                fh.write(data)
        return
    rng = random.Random(a.seed)
    for i in range(a.n):
        with open(os.path.join(a.out, f"fuzz_{a.seed}_{i:04d}.yaml"), "wb") as fh:
            fh.write(document(rng).encode("utf-8"))


if __name__ == "__main__":
    main()
