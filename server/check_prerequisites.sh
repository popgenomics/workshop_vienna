#!/usr/bin/env bash
# Non-destructive prerequisite check for the Vienna workshop server.
# Does not write into the repository. Does not run IQ-TREE, ASTRAL, or a DILS analysis.

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
TMPDIR_CHECK=""
blocking=0
warnings=0

note() { printf '%s\n' "$*"; }
fail() { blocking=$((blocking + 1)); printf 'ERROR: %s\n' "$*"; }
warn() { warnings=$((warnings + 1)); printf 'WARNING: %s\n' "$*"; }

cleanup() {
  if [[ -n "${TMPDIR_CHECK}" && -d "${TMPDIR_CHECK}" ]]; then
    rm -rf "${TMPDIR_CHECK}"
  fi
}
trap cleanup EXIT

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

note "Repository root: ${REPO_ROOT}"
note "Host: $(uname -s) $(uname -m)"
note "Kernel: $(uname -r)"

if have_cmd bash; then
  note "bash: $(bash --version | head -n 1)"
else
  fail "bash is not available"
fi

if have_cmd git; then
  note "git: $(git --version)"
else
  fail "git is not available"
fi

if have_cmd R; then
  note "R: $(R --vanilla --version | head -n 1)"
else
  fail "R is not available"
fi

if have_cmd Rscript; then
  note "Rscript: $(command -v Rscript)"
else
  fail "Rscript is not available"
fi

if have_cmd gcc; then
  note "gcc: $(gcc --version | head -n 1)"
else
  fail "gcc is not available"
fi

if have_cmd python3; then
  note "python3: $(python3 --version 2>&1)"
else
  fail "python3 is not available"
fi

if have_cmd python3; then
  py_st=0
  py_out="$(python3 - <<'PY' 2>&1
import sys
try:
    import pandas as pd
except Exception as exc:
    sys.stderr.write("ERROR: cannot import pandas: %s\n" % exc)
    sys.exit(2)
print("pandas: %s" % pd.__version__)
parts = []
for token in pd.__version__.split("."):
    digits = "".join(ch for ch in token if ch.isdigit())
    parts.append(int(digits) if digits else 0)
while len(parts) < 3:
    parts.append(0)
if tuple(parts[:3]) > (2, 1, 0):
    sys.stderr.write("WARNING: pandas %s is greater than 2.1.0\n" % pd.__version__)
    sys.exit(1)
sys.exit(0)
PY
)" || py_st=$?
  printf '%s\n' "${py_out}"
  if [[ "${py_st}" -eq 0 ]]; then
    note "pandas import and version check: ok"
  elif [[ "${py_st}" -eq 2 ]]; then
    fail "pandas is not importable"
  else
    warn "pandas is installed but its version is greater than 2.1.0"
  fi
fi

arch="$(uname -m)"
if [[ "${arch}" == "x86_64" || "${arch}" == "amd64" ]]; then
  note "Architecture: ${arch} (Linux x86-64 expected)"
else
  warn "Architecture is ${arch}; Aphid Day 2 expects Linux x86-64."
fi

if have_cmd getconf; then
  note "glibc (getconf GNU_LIBC_VERSION): $(getconf GNU_LIBC_VERSION 2>/dev/null || echo unknown)"
elif have_cmd ldd; then
  note "glibc (ldd --version): $(ldd --version 2>&1 | head -n 1)"
else
  warn "Could not determine glibc version."
fi

if have_cmd iqtree3; then
  iqtree3_bin="$(command -v iqtree3)"
  note "iqtree3: ${iqtree3_bin}"
  note "IQ-TREE was not executed (presence check only)."
  iqtree3_dir="$(dirname "${iqtree3_bin}")"
  if [[ -f "${iqtree3_bin}" ]] && grep -q 'iqtree3_intel\|iqtree3_arm' "${iqtree3_bin}" 2>/dev/null; then
    if [[ ! -e "${iqtree3_dir}/iqtree3_intel" && ! -e "${iqtree3_dir}/iqtree3_arm" ]]; then
      warn "iqtree3 looks like a wrapper, but iqtree3_intel / iqtree3_arm was not found next to it."
    fi
  fi
else
  warn "iqtree3 is not in PATH. Required for Day 2."
fi

if have_cmd astral; then
  astral_bin="$(command -v astral)"
  note "astral: ${astral_bin}"
  note "Exact ASTRAL implementation and version: pending confirmation from Arthur Boddaert."
  if [[ -f "${astral_bin}" ]] && grep -q -E -- '-i|-o' "${astral_bin}" 2>/dev/null; then
    note "The astral command file mentions -i/-o flags (not executed)."
  else
    warn "Could not confirm -i/-o support without running ASTRAL. Expected syntax: astral -i <input> -o <output>."
  fi
  warn "Presence of 'astral' is not a scientific validation of the Day 2 species-tree step."
else
  warn "astral is not in PATH. Exact ASTRAL implementation and version: pending confirmation from Arthur Boddaert."
fi

TMPDIR_CHECK="$(mktemp -d "${TMPDIR:-/tmp}/workshop_vienna_check.XXXXXX")"
note "Temporary directory: ${TMPDIR_CHECK}"
if [[ ! -w "${TMPDIR_CHECK}" ]]; then
  fail "Temporary directory is not writable: ${TMPDIR_CHECK}"
else
  probe="${TMPDIR_CHECK}/write_test"
  if printf 'ok\n' > "${probe}" && [[ -f "${probe}" ]]; then
    note "Write test in temporary directory: ok"
  else
    fail "Cannot write to temporary directory ${TMPDIR_CHECK}"
  fi
fi

R_CHECK="${SCRIPT_DIR}/check_R_packages.R"
if [[ ! -f "${R_CHECK}" ]]; then
  fail "Missing ${R_CHECK}"
elif have_cmd Rscript; then
  note "Running ${R_CHECK}"
  if Rscript --vanilla "${R_CHECK}"; then
    note "R package check: ok"
  else
    fail "R package check failed"
  fi
fi

SRC="${REPO_ROOT}/Aphid/software/aphid.0.11.c"
if [[ ! -f "${SRC}" ]]; then
  fail "Missing Aphid source: ${SRC}"
elif have_cmd gcc; then
  note "Compiling Aphid in the temporary directory (not in the repository)"
  if gcc -O2 -std=c11 -Wall -Wextra "${SRC}" -lm -o "${TMPDIR_CHECK}/aphid" 2>"${TMPDIR_CHECK}/gcc.err"; then
    warn_n="$(grep -c 'warning:' "${TMPDIR_CHECK}/gcc.err" || true)"
    note "Aphid test compile: ok (${warn_n} compiler warning lines; expected, source not modified)"
    if [[ -x "${TMPDIR_CHECK}/aphid" ]]; then
      usage="$("${TMPDIR_CHECK}/aphid" 2>&1 || true)"
      note "Aphid usage: ${usage}"
      if printf '%s\n' "${usage}" | grep -q 'usage: aphid'; then
        note "Aphid prints a usage message without arguments."
      else
        warn "Aphid did not print the expected usage line."
      fi
      arg_n="$(printf '%s\n' "${usage}" | grep -o 'usage: aphid .*' | awk '{print NF-2}')"
      if [[ "${arg_n}" == "4" ]]; then
        note "Aphid usage lists 4 arguments."
      else
        warn "Could not confirm four-argument usage from: ${usage}"
      fi
    else
      fail "Aphid test compile did not produce an executable"
    fi
  else
    fail "Aphid test compile failed"
  fi
fi

note ""
note "Blocking errors: ${blocking}; warnings: ${warnings}"
note "This script does not run IQ-TREE, ASTRAL, or a DILS random forest."
if [[ "${blocking}" -gt 0 ]]; then
  exit 1
fi
exit 0
