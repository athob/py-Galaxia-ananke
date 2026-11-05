# Security policy

## Supported versions

Security fixes go into the latest release only.

## Reporting a vulnerability

Please **do not open a public issue** for a security problem. Use GitHub's
private reporting instead: on the repository page, open the **Security** tab and
choose **Report a vulnerability**. Say what you found, how to reproduce it, and
which version you were using. If you cannot use GitHub for this, email the
maintainer at <adrien.thob@gmail.com>.

You can expect an acknowledgement within a few days. A fix is released as a new
version, and the report is credited in the release notes unless you prefer
otherwise.

## What counts

Installing py-Galaxia-ananke fetches the Galaxia-ananke C++ source (a git
submodule of this repository, pinned to one commit) and builds it, and using the
package runs the resulting `galaxia` executable on your data, so these are the
areas where a bug is a security matter:

- the installation fetching or building something other than the intended
  Galaxia-ananke source, which is the commit this repository pins;
- input that reaches a shell, or the executable's command line, without being
  escaped;
- files being read or written outside of the directories the wrapper was pointed
  at, and its own cache directory;
- a published release (the sdist on PyPI) that does not match the tagged source.

A crash, a wrong result or a slow run is an ordinary bug: open an issue.
