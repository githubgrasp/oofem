#!/usr/bin/env bash
#
# Build the HTML documentation for Peter Grassl's OOFEM fork and stage it as a
# single tree ready to publish under petergrassl.com/oofem/.
#
# Usage:
#   ./build-fork-docs.sh [STAGE_DIR]
#
#   STAGE_DIR  where to write the site (default: doc/_fork-site, a local preview).
#              Point it straight at your Hugo static path to stage for the site:
#                  ./build-fork-docs.sh ~/petergrassl.com/static/oofem
#
# Each Sphinx manual is built with `sphinx-build -b html` into STAGE_DIR/<name>,
# and a landing index.html (clearly marked as the fork) links them together.
# Afterwards deploy with your usual rclone ... --checksum (many small files).

set -euo pipefail

DOC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="${1:-$DOC_DIR/_fork-site}"

command -v sphinx-build >/dev/null 2>&1 || {
    echo "ERROR: sphinx-build not found. Install with: pip install sphinx" >&2
    exit 1
}

# <source-path (relative to doc/)>|<output-name>|<title on the landing page>
# (order = order on the page). Fork tool manuals live with their code under src/.
MANUALS=(
    "usermanual|usermanual|User Manual"
    "oofemInput|oofemInput|Input Manual"
    "matlibmanual|matlibmanual|Material Library Manual"
    "elementlibmanual|elementlibmanual|Element Library Manual"
    "theorymanual2|theorymanual2|Theory Manual"
    "developerGuide|developerGuide|Developer Guide"
    "python|python|Python Bindings"
    "../src/aggregate/doc|aggregate|Aggregate Tool"
    "../src/converter/doc|converter|Converter Tool"
    "../src/generator/doc|generator|Generator Tool"
)

echo "Staging OOFEM fork docs into: $STAGE"
mkdir -p "$STAGE"

links=""
for entry in "${MANUALS[@]}"; do
    srcpath="${entry%%|*}"
    rest="${entry#*|}"
    outname="${rest%%|*}"
    title="${rest#*|}"
    src="$DOC_DIR/$srcpath"
    if [ ! -f "$src/conf.py" ]; then
        echo "  skip  $outname (no conf.py at $src)"
        continue
    fi
    echo "  build $outname -> $STAGE/$outname"
    rm -rf "${STAGE:?}/$outname"
    sphinx-build -q -b html "$src" "$STAGE/$outname"
    links+="        <li><a href=\"$outname/index.html\">$title</a></li>"$'\n'
done

built="$(date +%Y-%m-%d)"

# Landing page. Unquoted heredoc so $links/$built expand; the template has no
# other $ or backticks. Keep the fork statement prominent and link upstream.
cat > "$STAGE/index.html" <<HTML
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>OOFEM Documentation &mdash; Peter Grassl's fork</title>
<style>
  body { font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
         max-width: 52rem; margin: 3rem auto; padding: 0 1.2rem; line-height: 1.55; color: #222; }
  h1 { margin-bottom: .2rem; }
  .sub { color: #666; margin-top: 0; }
  .fork { background: #fff8e1; border: 1px solid #f0d890; border-radius: 6px;
          padding: .8rem 1rem; margin: 1.4rem 0; }
  ul { list-style: none; padding: 0; }
  li { margin: .55rem 0; font-size: 1.08rem; }
  a { color: #1565c0; text-decoration: none; }
  a:hover { text-decoration: underline; }
  footer { margin-top: 2.5rem; color: #888; font-size: .9rem; }
</style>
</head>
<body>
  <h1>OOFEM Documentation</h1>
  <p class="sub">Peter Grassl's fork</p>

  <div class="fork">
    <strong>This is the documentation for Peter Grassl's fork of OOFEM.</strong>
    It includes additional lattice and discrete material/element models (and related
    tooling) that are not part of upstream OOFEM. For the upstream project see
    <a href="https://www.oofem.org/">oofem.org</a>
    (<a href="https://github.com/oofem/oofem">github.com/oofem/oofem</a>).
    Fork source: <a href="https://github.com/githubgrasp/oofem">github.com/githubgrasp/oofem</a>.
  </div>

  <h2>Manuals</h2>
  <ul>
$links  </ul>

  <footer>Built $built with Sphinx. Maintained by Peter Grassl &middot;
    <a href="https://www.petergrassl.com/">petergrassl.com</a>.</footer>
</body>
</html>
HTML

echo
echo "Done. Landing page: $STAGE/index.html"
echo "Preview locally:    open \"$STAGE/index.html\""
echo
echo "To publish: point STAGE_DIR at your Hugo static/oofem path (or copy $STAGE there),"
echo "then deploy with rclone ... --checksum."
