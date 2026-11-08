# normpath.awk - Python posixpath path functions in gawk (PORT_SPEC section 3, "Paths"; schema.py:240-245).
# Library: load with `gawk -f normpath.awk ...`. Pure text, no symlink resolution (realpath is never used).
#
#   py_normpath(p)            = posixpath.normpath(p): '' -> '.'; keeps exactly two leading slashes, 3+ -> one;
#                               drops '' and '.' components; '..' pops textually, a leading '..' of a relative path
#                               is kept, '..' directly under the root is dropped
#   py_expanduser(p)          = posixpath.expanduser(p): '~' -> $HOME (getent passwd of the uid when HOME is unset),
#                               '~user' -> pw_dir from `getent passwd user`; unknown user -> p unchanged
#   py_isabs(p)               = p starts with '/'
#   py_join(a, b)             = posixpath.join(a, b) for two components
#   py_dirname(p)             = posixpath.dirname(p)
#   py_abspath(p, cwd)        = normpath(join(cwd, p)) for a relative p (cwd = os.getcwd(), i.e. `pwd -P`)
#   satpipe_norm_path(v, base)= schema.py _norm_path: expanduser, join to base when relative, normpath, '\' -> '/'

function py_isabs(p) { return substr(p, 1, 1) == "/" }

function py_join(a, b) {
    if (py_isabs(b) || a == "") return b
    return (substr(a, length(a)) == "/") ? a b : a "/" b
}

function py_normpath(p,   init, comps, n, i, c, st, top, out) {
    if (p == "") return "."
    init = 0
    if (substr(p, 1, 1) == "/") init = (substr(p, 1, 2) == "//" && substr(p, 1, 3) != "///") ? 2 : 1
    n = split(p, comps, "/")
    top = 0
    for (i = 1; i <= n; i++) {
        c = comps[i]
        if (c == "" || c == ".") continue
        if (c != ".." || (!init && top == 0) || (top > 0 && st[top] == "..")) st[++top] = c
        else if (top > 0) top--
    }
    out = ""
    for (i = 1; i <= top; i++) out = out (i > 1 ? "/" : "") st[i]
    if (init == 1) out = "/" out
    else if (init == 2) out = "//" out
    return (out == "") ? "." : out
}

function py_dirname(p,   i, head) {
    i = match(p, /\/[^\/]*$/)
    if (!i) return ""
    head = substr(p, 1, i)
    if (head !~ /^\/+$/) sub(/\/+$/, "", head)
    return head
}

function py_abspath(p, cwd) {
    if (!py_isabs(p)) p = py_join(cwd, p)
    return py_normpath(p)
}

function _np_shq(s) { gsub(/'/, "'\\''", s); return "'" s "'" }

# home directory of a user (name "" = the current uid); sets _NP_OK = 0 when the lookup fails
function _np_pwdir(name,   cmd, line, f, n) {
    _NP_OK = 0
    if (name == "") cmd = "getent passwd \"$(id -u)\" 2>/dev/null"
    else cmd = "getent passwd -- " _np_shq(name) " 2>/dev/null"
    line = ""
    if ((cmd | getline line) > 0) {
        n = split(line, f, ":")
        if (n >= 7) _NP_OK = 1
    }
    close(cmd)
    return _NP_OK ? f[6] : ""
}

function py_expanduser(p,   i, home, name) {
    if (substr(p, 1, 1) != "~") return p
    i = index(substr(p, 2), "/"); i = (i == 0) ? length(p) + 1 : i + 1      # 1-based index of the first '/' after '~'
    if (i == 2) {
        if ("HOME" in ENVIRON) home = ENVIRON["HOME"]
        else { home = _np_pwdir(""); if (!_NP_OK) return p }
    } else {
        name = substr(p, 2, i - 2)
        if (index(name, "\000")) return p
        home = _np_pwdir(name); if (!_NP_OK) return p
    }
    sub(/\/+$/, "", home)                                                    # userhome.rstrip('/')
    home = home substr(p, i)
    return (home == "") ? "/" : home
}

function satpipe_norm_path(v, base,   x) {
    x = py_expanduser(v)
    if (base != "" && !py_isabs(x)) x = py_join(base, x)
    x = py_normpath(x)
    gsub(/\\/, "/", x)
    return x
}
