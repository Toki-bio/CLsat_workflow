# st_02_trf.sh - stage 02_trf (run.py:344-410): TRF per scaffold in parallel, scaffolds above
# finalisation.large_scaffold_bp left to 02b; given TRF outputs (inputs.trf_dat) linked in.
. "$SATPIPE_LIB/st_common_wp6.sh"

plan_02_trf() {                                                                # run.py:405-410
    local -a k=()
    local c
    if [ "$(cfg_type inputs/trf_dat 2>/dev/null)" = map ]; then
        mapfile -t k < <(cfg_keys inputs/trf_dat)
        for c in "${k[@]}"; do emit_inputs "$(cfg "inputs/trf_dat/$c")"; done
    fi
    mapfile -t k < <(cfg_keys genomes)
    for c in "${k[@]}"; do emit_inputs "$(cfg "genomes/$c").fai"; done
    genome_files
    cfg_subtree discovery/trf trf
    cfg_subtree finalisation/large_scaffold_bp large                           # int(...) of an int field
    cfg_subtree genomes genomes
}

st_02_trf() {
    local trf_dir trfp thr jobs code genome dat split f b n rc cmd big msg
    local -a outs=() dats=() failed=()
    local -A given=() isg=()
    trf_dir=$(ctx_d trf)
    trfp=$(cfg discovery/trf)
    thr=$(cfg finalisation/large_scaffold_bp)
    jobs=$(cfg expansion/jobs)
    wp6_codes
    if [ "$(cfg_type inputs/trf_dat 2>/dev/null)" = map ]; then
        local -a gk=()
        mapfile -t gk < <(cfg_keys inputs/trf_dat)
        for code in "${gk[@]}"; do given[$code]=$(cfg "inputs/trf_dat/$code"); done
    fi
    for code in "${WP6_CODES[@]}"; do isg[$code]=1; done
    # stale .dat of a genome no longer in the config (run.py:355-358), in os.listdir (= readdir, find) order like Python
    local old
    mapfile -d '' -t dats < <(find "$trf_dir" -mindepth 1 -maxdepth 1 -name '*.dat' -print0)
    for f in "${dats[@]}"; do
        b=${f##*/}; b=${b%.dat}
        if [ -z "$b" ] || [ -z "${isg[$b]:-}" ]; then
            rm -f -- "$f" || stage_die "cannot remove $f"
            say "  removed ${f##*/}: not a genome of this config"
        fi
    done
    big=$RUN_TMP/wp6_big_scaffolds
    : > "$big"
    for code in "${WP6_CODES[@]}"; do
        genome=$(cfg "genomes/$code")
        dat=$trf_dir/$code.dat
        # scaffolds longer than the threshold, from the .fai (run.py:349-352, 361-364)
        C=$code T=$thr "$GAWK" -F '\t' '{ if ($2 + 0 > ENVIRON["T"] + 0) printf "%s\t%s\t%d\n", ENVIRON["C"], $1, $2 }' "$genome.fai" >> "$big" \
            || stage_die "[Errno 2] No such file or directory: '$genome.fai'"
        if [ -n "${given[$code]+x}" ]; then                                    # run.py:365-371
            rm -f -- "$dat" || stage_die "cannot remove $dat"
            ln -s -- "$(F=${given[$code]} CWD=$(pwd -P) "$GAWK" -f "$SATPIPE_LIB/normpath.awk" -e 'BEGIN { printf "%s", py_abspath(ENVIRON["F"], ENVIRON["CWD"]) }')" "$dat" \
                || stage_die "cannot link $dat"
            say "  $code: using given TRF output ${given[$code]}"
            outs+=("$dat")
            continue
        fi
        split=$(ctx_d trf "$code.split")
        # split the genome per scaffold, large ones skipped; Python's text mode reads CRLF as LF (run.py:374-392)
        n=$(C=$code BIG=$big D=$split "$GAWK" '
            BEGIN { RS = "\n"; while ((getline l < ENVIRON["BIG"]) > 0) { split(l, p, "\t"); if (p[1] == ENVIRON["C"]) L[p[2]] = 1 }
                    close(ENVIRON["BIG"]); out = ""; nn = 0 }
            /^>/ {
                h = substr($0, 2); sub(/^[ \t\n\r\f\v\034-\037]+/, "", h)
                if (h == "") { print "IndexError: list index out of range" > "/dev/stderr"; exit 3 }
                match(h, /[ \t\n\r\f\v\034-\037]/); name = RSTART ? substr(h, 1, RSTART - 1) : h
                skip = (name in L)
                if (out != "") close(out)
                out = ""
                if (!skip) { N[++nn] = name; out = ENVIRON["D"] "/" name ".fa"; printf ">%s\n", name > out }
                next
            }
            out != "" { if (sub(/\r$/, "")) printf "%s\n", $0 > out; else printf "%s%s", $0, RT > out }
            END { if (out != "") close(out); j = ENVIRON["D"] "/jobs.txt"; printf "" > j
                  for (i = 1; i <= nn; i++) print N[i] ".fa" > j
                  close(j); print nn }' "$genome") || stage_die "02_trf $code: cannot split $genome"
        # run.py:383-386: the identical parallel command string
        cmd=$(printf 'cat jobs.txt | %s --jobs %d "%s {} %s -h > /dev/null 2> {}.err; rc=\\$?; if [ \\$rc -ge 125 ]; then echo \\$rc > {}.failed; fi; exit 0"' \
              "$(tool parallel)" "$jobs" "$(tool trf)" "$trfp")
        sh_run 02_trf "$split" log_name="02_trf_$code" -- bash -c "$cmd"
        # crashed TRF runs (run.py:387-390), names sorted as Python strings
        old=$(shopt -p nullglob dotglob || true)
        shopt -s nullglob dotglob
        failed=("$split"/*.fa.failed)
        eval "$old"
        if [ ${#failed[@]} -gt 0 ]; then
            msg=$(for f in "${failed[@]}"; do b=${f##*/}; printf '%s\0' "${b%.fa.failed}"; done | "$GAWK" -f "$PYJSON" -e '
                BEGIN { RS = "\0" } { S[$0] = 1; n++ }
                END { PROCINFO["sorted_in"] = "@ind_str_asc"; for (k in S) if (++i <= 5) o = o (i > 1 ? ", " : "") py_repr_str(k)
                      printf "%d scaffolds (exit >= 125), e.g. [%s]", n, o }')
            stage_die "02_trf $code: TRF crashed on $msg; see $split/<scaffold>.fa.err"
        fi
        # concatenate the parts in jobs order; missing parts fail the stage after the others are written (run.py:391-401)
        rc=0
        msg=$(D=$split SUF=".$(printf '%s' "$trfp" | tr ' ' '.').dat" OUT=$dat "$GAWK" -f "$PYJSON" -e '
            BEGIN { RS = "\n"; out = ENVIRON["OUT"]; printf "" > out; j = ENVIRON["D"] "/jobs.txt"; nm = 0
                    while ((getline l < j) > 0) {
                        part = ENVIRON["D"] "/" l ENVIRON["SUF"]
                        r = (getline x < part)
                        if (r < 0) { if (++nm <= 5) ms = ms (nm > 1 ? ", " : "") py_repr_str(substr(l, 1, length(l) - 3)); continue }
                        while (r > 0) { if (sub(/\r$/, "", x)) printf "%s\n", x > out; else printf "%s%s", x, RT > out; r = (getline x < part) }
                        close(part)
                    }
                    close(out)
                    if (nm) { printf "%d scaffolds, e.g. [%s]", nm, ms; exit 3 } }') || rc=$?
        if [ $rc -ne 0 ]; then stage_die "02_trf $code: TRF produced no output for $msg"; fi
        rm -rf -- "$split" || true
        say "  $code: TRF on $n scaffolds -> $dat"
        outs+=("$dat")
    done
    {   printf '# scaffolds longer than finalisation.large_scaffold_bp=%d; scanned in chunks by stage 02b\n' "$thr"
        cat -- "$big"
    } > "$trf_dir/large_scaffolds.tsv"
    emit_outputs "${outs[@]}" "$trf_dir/large_scaffolds.tsv"
}
