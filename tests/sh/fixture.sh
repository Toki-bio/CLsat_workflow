# fixture.sh - the `project` fixture of tests/conftest.py:24-48 as a bash function (sourced by tests/sh/test_check.sh and
# tests/sh/test_schema.sh). make_project DIR writes the same files byte for byte (textwrap.dedent of the conftest text);
# the config is DIR/configs/run.yaml. conftest.py:12-17 lifts check.MAX_BASE_PATH for every test: SATPIPE_MAX_BASE_PATH=10000.

export SATPIPE_MAX_BASE_PATH=10000

_fx_rep() { local s=$1 n=$2 o=""; while [ "$n" -gt 0 ]; do o+=$s; n=$((n - 1)); done; printf '%s' "$o"; }

make_project() {
    local d=$1 code
    mkdir -p "$d/genomes" "$d/trf" "$d/ref" "$d/configs"
    for code in spa spb; do
        printf '>%s_scf1\n%s\n' "$code" "$(_fx_rep ACGT 50)" > "$d/genomes/$code.fna"
        printf 'Sequence: scf1\n\n1 500 146 3.4 146 95 2 300 ACGT ACGTACGT\n' > "$d/trf/$code.dat"
    done
    printf '>SatX\n%sACGT\n' "$(_fx_rep ACGTTGCA 18)" > "$d/ref/SatX.cons.fa"          # 148 bp
    printf '# cluster\tfamily\n0\tSatX\n5\tSatX\tnote\n' > "$d/ref/cluster_family.tsv"
    printf 'spa_scf1\t0\t1000\tspa|spa_scf1:1-1000|MANUAL\t0\t.\n' > "$d/ref/manual.bed"
    cat > "$d/configs/run.yaml" <<'EOF'

run_name: toy
out_dir: ../runs/toy
genomes:
  spa: ../genomes/spa.fna
  spb: ../genomes/spb.fna
inputs:
  trf_dat:
    spa: ../trf/spa.dat
    spb: ../trf/spb.dat
families:
  SatX: {consensus: ../ref/SatX.cons.fa, mono: 148}
family_assignment: ../ref/cluster_family.tsv
manual_loci: ../ref/manual.bed
EOF
}

# edit FILE OLD NEW      cfg.write_text(cfg.read_text().replace(OLD, NEW)) (every occurrence, like str.replace)
edit() { local t; t=$(cat "$1"; printf x); t=${t%x}; printf '%s' "${t//"$2"/"$3"}" > "$1"; }
# append FILE TEXT       cfg.write_text(cfg.read_text() + TEXT)
append() { printf '%s' "$2" >> "$1"; }
