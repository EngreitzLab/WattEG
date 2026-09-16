#!/usr/bin/env bash
#
# Activate .githooks/pre-commit, by whichever mechanism this git supports.
#
# WHY THIS IS NOT ONE LINE. `git config core.hooksPath .githooks` is the modern answer and it is
# what every README suggests -- but core.hooksPath arrived in **git 2.9**, and Sherlock's default
# git is **1.8.3.1**, which parses the setting, stores it, reports it back from `git config --get`,
# and then ignores it completely. So the one-liner leaves you with a repo that looks hooked and
# silently is not: `git config --get core.hooksPath` answers `.githooks`, and nothing runs.
#
# That is exactly how this repo's hook went unused from the day it was written: `git config --get
# core.hooksPath` answered `.githooks` while no hook had ever run. Found 2026-09-16.
#
# So: set core.hooksPath for the gits that honour it, AND symlink into .git/hooks, which every
# version has honoured since git existed. The symlink is what actually works here.
#
# Idempotent. Safe to re-run, and safe to run before or after `ml system git/2.45.1`.

set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

git_version="$(git --version | awk '{print $3}')"
git_major="${git_version%%.*}"
rest="${git_version#*.}"
git_minor="${rest%%.*}"

chmod +x .githooks/pre-commit

# 1. core.hooksPath, for git >= 2.9.
git config core.hooksPath .githooks

# 2. The symlink, which works everywhere. Relative, so the repo stays movable.
hooks_dir="$(git rev-parse --git-dir)/hooks"
mkdir -p "$hooks_dir"
target="../../.githooks/pre-commit"
link="${hooks_dir}/pre-commit"
if [ -e "$link" ] && [ ! -L "$link" ]; then
    echo "WARNING: ${link} exists and is not a symlink. Leaving it alone." >&2
    echo "         Move it aside and re-run if you want the shared hook instead." >&2
else
    ln -sfn "$target" "$link"
fi

echo "pre-commit hook installed:"
echo "  git             ${git_version}"
if [ "$git_major" -gt 2 ] || { [ "$git_major" -eq 2 ] && [ "$git_minor" -ge 9 ]; }; then
    echo "  core.hooksPath  .githooks   (honoured by this git)"
else
    echo "  core.hooksPath  .githooks   (SET BUT IGNORED -- needs git >= 2.9, you have ${git_version})"
fi
echo "  symlink         ${link} -> ${target}   (what actually runs here)"
echo
echo "Verify with:  git commit --allow-empty -m 'hook check'   (you should see the checks run)"
