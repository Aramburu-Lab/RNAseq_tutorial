#!/usr/bin/env bash
# /// codebase-meta
# exempt: >-
#   Sourced helper library, not a pipeline step. It defines functions and writes
#   nothing of its own, so it has no outputs to describe and no place in the DAG.
# ///
# Shared bash helpers. Source it, never execute it:
#
#     source "$(dirname "$0")/../../bin/lib_slurm.sh"
#
# This file is copied into the project rather than imported, because bash has no
# import mechanism and cannot reach a pixi-installed package. ADR-4's "logic is
# never copied" holds for every language that can import; for this one it cannot,
# so the copy is deliberate and arrives through Copier like any other scaffold
# file. Do not fork it inside a script — fix it here and `copier update`.
#
# Everything below is POSIX-portable on purpose: awk only, no `du -sb`, no
# `stat -c`, no `numfmt`. Those are GNU-only and fail on macOS, which is where
# half of this code gets debugged.
#
# Conventions it encodes, lifted from the four SHARED_*.sh scripts that already
# work on Kebnekaise:
#   - long --snake_case flags only, no short forms
#   - exit 2 = usage error, exit 1 = runtime failure
#   - one flag declaration drives usage text, parsing, validation and re-submission

# One record per option: --flag|VAR_NAME|default|req|opt|Description
# Declare it in the calling script before calling anything here. Three separate
# copies of the flag list — header comment, usage heredoc, sbatch re-invocation —
# is how they drift, and the lab already has that bug in production.
OPT_SPEC=()

# Print each spec row belonging to one class ("req" or "opt"), with the current
# value of the variable interpolated so the help text cannot go stale.
_spec_rows() {
    local want="$1" rec flag var def req desc
    for rec in ${OPT_SPEC[@]+"${OPT_SPEC[@]}"}; do
        IFS='|' read -r flag var def req desc <<<"$rec"
        [ "$req" = "$want" ] || continue
        if [ "$want" = "req" ]; then
            printf '  %-22s %s\n' "$flag" "$desc"
        else
            printf '  %-22s %s (default: %s)\n' "$flag" "$desc" "${!var:-$def}"
        fi
    done
}

# Render the whole usage message from OPT_SPEC.
usage_from_spec() {
    local me
    me="$(basename "$0")"
    printf 'Usage: %s [OPTIONS]\n' "$me"
    printf '       sbatch %s [OPTIONS]\n\n' "$me"
    printf 'Required:\n'
    _spec_rows req
    printf '\nOptional:\n'
    _spec_rows opt
    printf '  %-22s %s\n' "-h, --help" "Show this help and exit"
}

# Look up the variable a flag writes to. Prints it; returns 1 for an unknown flag.
_spec_var() {
    local rec flag var def req desc
    for rec in ${OPT_SPEC[@]+"${OPT_SPEC[@]}"}; do
        IFS='|' read -r flag var def req desc <<<"$rec"
        if [ "$flag" = "$1" ]; then
            printf '%s' "$var"
            return 0
        fi
    done
    return 1
}

# Seed every spec variable with its declared default, unless already set.
_apply_defaults() {
    local rec flag var def req desc
    for rec in ${OPT_SPEC[@]+"${OPT_SPEC[@]}"}; do
        IFS='|' read -r flag var def req desc <<<"$rec"
        [ -n "${!var:-}" ] || printf -v "$var" '%s' "$def"
    done
}

# Parse "$@" against OPT_SPEC. Unknown flag or missing value exits 2.
parse_args() {
    _apply_defaults
    local var
    while [ "$#" -gt 0 ]; do
        case "$1" in
            -h|--help) usage_from_spec; exit 0 ;;
            --) shift; break ;;
            --*)
                var="$(_spec_var "$1")" || {
                    echo "ERROR: unknown argument: $1" >&2
                    usage_from_spec >&2
                    exit 2
                }
                [ "$#" -ge 2 ] || { echo "ERROR: $1 needs a value" >&2; exit 2; }
                printf -v "$var" '%s' "$2"
                shift 2
                ;;
            *) echo "ERROR: unexpected argument: $1" >&2; usage_from_spec >&2; exit 2 ;;
        esac
    done
}

# Report every missing required flag at once, not one per run.
require_args() {
    local rec flag var def req desc
    local missing=()
    for rec in ${OPT_SPEC[@]+"${OPT_SPEC[@]}"}; do
        IFS='|' read -r flag var def req desc <<<"$rec"
        [ "$req" = "req" ] || continue
        [ -n "${!var:-}" ] || missing+=("$flag")
    done
    if [ ${#missing[@]} -gt 0 ]; then
        echo "ERROR: missing required argument(s): ${missing[*]}" >&2
        usage_from_spec >&2
        exit 2
    fi
}

# Re-submit this script to Slurm with the flags it was given, unless it is
# already running under Slurm or there is no scheduler here.
#
# The flag list is rebuilt from OPT_SPEC, so adding an option never leaves a
# stale re-invocation behind. `mkdir -p logs/slurm` runs BEFORE the exec on
# purpose: Slurm opens the --output file before the first line of the script
# runs, so a committed logs/.gitkeep only helps inside the generated tree — a
# script invoked from an arbitrary data directory would fail before bash starts.
slurm_selfsubmit() {
    local account="${1:-}" rec flag var def req desc value
    local args=()
    [ -n "${SLURM_JOB_ID:-}" ] && return 0
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "No sbatch here — running locally." >&2
        return 0
    fi
    [ -n "$account" ] || {
        echo "ERROR: no Slurm account. Pass --use_this_account, or set one in the script." >&2
        exit 2
    }
    for rec in ${OPT_SPEC[@]+"${OPT_SPEC[@]}"}; do
        IFS='|' read -r flag var def req desc <<<"$rec"
        value="${!var:-}"
        [ -n "$value" ] && args+=("$flag" "$value")
    done
    mkdir -p logs/slurm
    echo "Not under Slurm — self-submitting (account: $account)" >&2
    exec sbatch -A "$account" "$0" ${args[@]+"${args[@]}"}
}

# Make a tool available, or fail. The module system is OPTIONAL — a laptop has
# none — but the tool is NOT. Never `module load` unguarded: it is the single
# most common reason a lab script cannot be run locally.
#
#   need_tool zstd GCCcore/13.2.0 zstd
need_tool() {
    local tool="$1"
    shift
    if [ "$#" -gt 0 ] && command -v module >/dev/null 2>&1; then
        module load "$@" 2>/dev/null || true
    fi
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: $tool not found in PATH" >&2
        exit 1
    }
}

# POSIX du in kilobytes; silenced stderr for unreadable subtrees.
du_kb() { du -sk "$1" 2>/dev/null | awk '{print $1}'; }

# Human-readable from KB: K / M / G / T
human() {
    awk -v k="$1" 'BEGIN{
        if (k+0 == 0)        { print "0"; exit }
        if (k>=1073741824)     printf "%.1fT\n", k/1073741824;
        else if (k>=1048576)   printf "%.1fG\n", k/1048576;
        else if (k>=1024)      printf "%.1fM\n", k/1024;
        else                   printf "%dK\n",    k;
    }'
}
