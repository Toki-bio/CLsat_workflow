# schema_doc.awk -- `satpipe config --schema-doc`: render lib/satpipe/schema.tsv as docs/config.md.
# Port of schema.py to_markdown/_type_str (schema.py:346-383); output is byte-identical to the Python one.
# Usage: LC_ALL=C gawk -f lib/satpipe/schema_doc.awk lib/satpipe/schema.tsv
# Rows come in schema definition order (depth first). A row is rendered when every ancestor is a block or the
# <key> value block of a map of blocks; blocks themselves print nothing; a map of blocks prints one
# `map of blocks` row and then its value block's fields; lists and other maps print one row and hide their
# <item>/<key> descendants.
BEGIN {
    FS = "\t"; OFS = ""; skip = ""; n = 0
}
/^#/ {                                   # header: column name -> index
    sub(/^#/, "")
    for (i = 1; i <= NF; i++) col[$i] = i
    next
}
{
    path = $col["path"]
    if (path == ".") next                  # the root block
    if (skip != "" && substr(path, 1, length(skip)) == skip) next
    skip = ""
    t = $col["type"]
    if (t == "block") next
    doc_path = path; gsub(/\//, ".", doc_path)
    req = ($col["required"] == "True") ? "yes" : ""
    if (t == "map" && $col["elem_type"] == "block") {
        rows[++n] = "| `" doc_path ".<key>` | map of blocks | " req " | | " $col["doc"] " |"
        next                               # the <key> block row is skipped (block); its fields follow
    }
    if (t == "map" || t == "list") skip = path "/"
    dflt = ($col["default_type"] == "null") ? "" : "`" $col["default_str"] "`"
    rows[++n] = "| `" doc_path "` | " type_str() " | " req " | " dflt " | " $col["doc"] " |"
}
function type_str(   s, a, k, i) {
    if (t == "enum") {                     # choices cell is a Python list repr of plain strings: ['plus', 'both']
        s = $col["choices"]; gsub(/^\[|\]$/, "", s)
        k = split(s, a, /, /)
        s = "one of "
        for (i = 1; i <= k; i++) { gsub(/^'|'$/, "", a[i]); s = s (i > 1 ? ", " : "") "`" a[i] "`" }
        return s
    }
    if (t == "list") return "list of " $col["elem_type"] ($col["length"] != "" ? " x" $col["length"] : "")
    if (t == "map") return "map -> " $col["elem_type"]
    return t
}
END {
    printf "%s", "# satpipe configuration reference\n\n" \
        "Generated from `satpipe/schema.py` by `satpipe config --schema-doc`; do not edit by hand.\n\n" \
        "Relative paths are resolved against the directory that holds the config file. Unknown keys " \
        "are errors. Defaults are the values of the published Darevskia run.\n\n" \
        "| key | type | required | default | meaning |\n|---|---|---|---|---|\n"
    for (i = 1; i <= n; i++) printf "%s%s", rows[i], (i < n ? "\n" : "")
    printf "\n"
}
