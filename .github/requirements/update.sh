#!/usr/bin/env bash
# Writes the hash-locked requirements of the CI and CD workflows, and the checksum of the
# actionlint download. Needs uv (https://docs.astral.sh/uv/) and curl.
#
#   bash .github/requirements/update.sh             lock what the .in files ask for
#   bash .github/requirements/update.sh --upgrade   and move everything to the newest versions
#
# requirements-X.in lists what a job installs. requirements-X.txt is its lock: every
# package, transitive ones included, pinned to a version and to the hashes of its files,
# so that `pip install --require-hashes -r requirements-X.txt` installs exactly that, or
# nothing. The locks are "universal": the same file for every platform, and for every
# Python from the one given below up (where the right version differs, a marker says so).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

upgrade=()
if [ "${1:-}" = "--upgrade" ]; then upgrade=(--upgrade); fi

lock() {  # lock NAME LOWEST_PYTHON [FOLDER]
  (cd "${3:-.}" && uv pip compile "requirements-$1.in" --output-file "requirements-$1.txt" \
    --universal --generate-hashes --python-version "$2" ${upgrade[@]+"${upgrade[@]}"})
}
lock lint 3.12
lock test 3.11   # the Pythons still supported: one version per package, which Dependabot needs
lock test 3.9 ../requirements-frozen   # the end-of-life Pythons, where Dependabot does not look
lock cd-build 3.12

# actionlint is not a Python package. Its version is ACTIONLINT_VERSION in ci.yml; the sha256
# of its release archive is kept here, and ci.yml refuses a download that does not match it.
sha256() { sha256sum "$@" 2>/dev/null || shasum -a 256 "$@"; }
version=$(sed -n 's/^ *ACTIONLINT_VERSION: "\(.*\)"$/\1/p' ../workflows/ci.yml)
[ -n "$version" ] || { echo "no ACTIONLINT_VERSION found in ci.yml" >&2; exit 1; }
base="https://github.com/rhysd/actionlint/releases/download/v${version}"
archive="actionlint_${version}_linux_amd64.tar.gz"
sum=$(curl -sSfL "${base}/${archive}" | sha256 | cut -d' ' -f1)
# the release publishes its own checksums too: the two have to agree
curl -sSfL "${base}/actionlint_${version}_checksums.txt" | grep -q "^${sum}  ${archive}$" \
  || { echo "the checksum of ${archive} is not the one the release publishes" >&2; exit 1; }
echo "$sum" > actionlint.sha256
echo "actionlint ${version}: ${sum}"
