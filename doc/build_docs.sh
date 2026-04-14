#!/bin/bash
# Convert rst documentation to html and markdown.
#
# Usage:
#   bash build_docs.sh
#
# Input:  *.rst files in the same directory as this script
# Output: html/ and md/ subdirectories

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

mkdir -p html md

for rst in *.rst; do
    [ -f "$rst" ] || continue
    base="${rst%.rst}"

    # Preprocess: convert :doc:`text <target>` to local links
    # For html: link to target.html
    # For md: link to target.md
    sed_html='s/:doc:`\([^<]*\) <\([^>]*\)>`/`\1 <\2.html>`__/g; s/:doc:`\([^`]*\)`/`\1 <\1.html>`__/g'
    sed_md='s/:doc:`\([^<]*\) <\([^>]*\)>`/`\1 <\2.md>`__/g; s/:doc:`\([^`]*\)`/`\1 <\1.md>`__/g'

    # Generate preprocessed tmp files
    sed "$sed_html" "$rst" > "/tmp/${base}_html.rst"
    sed "$sed_md" "$rst" > "/tmp/${base}_md.rst"

    # Convert to HTML (use MathJax CDN for proper math rendering)
    pandoc --from rst --to html5 --standalone \
           --mathjax="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js" \
           --metadata title="$base" \
           "/tmp/${base}_html.rst" -o "html/${base}.html" 2>/dev/null
    echo "  html/${base}.html"

    # Convert to Markdown (GitHub Flavored)
    pandoc --from rst --to gfm \
           "/tmp/${base}_md.rst" -o "md/${base}.md" 2>/dev/null
    echo "  md/${base}.md"

    rm -f "/tmp/${base}_html.rst" "/tmp/${base}_md.rst"
done

echo "Done. Output in html/ and md/"
