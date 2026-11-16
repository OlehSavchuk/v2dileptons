#!/usr/bin/env bash
# Creates the private GitHub repo "v2dileptons" from this folder and pushes it.
# Requires the GitHub CLI, already logged in:  gh auth login
set -euo pipefail
cd "$(dirname "$0")"
git rev-parse --git-dir >/dev/null            # must be run inside the repository
gh repo create v2dileptons --private --source=. --remote=origin --push
