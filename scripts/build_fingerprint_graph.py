"""Generate (or verify with --check) fccgroup/patterns/_dependency_graph.py."""

import argparse
import sys
from pathlib import Path

from fccgroup.patterns.dependencies import (
    extract_dependencies,
    find_cycle,
    misordered,
    render_graph_module,
)
from fccgroup.patterns.library import fingerprints

TARGET = Path(__file__).parent.parent / "fccgroup" / "patterns" / "_dependency_graph.py"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail if the file is out of date")
    args = parser.parse_args()

    graph = extract_dependencies(fingerprints)
    cycle = find_cycle(graph)
    if cycle:
        print("Dependency cycle: " + " -> ".join(cycle), file=sys.stderr)
        return 1
    out_of_order = misordered(fingerprints, graph)
    if out_of_order and not args.check:
        print("Note: defined before their dependencies (evaluated in dependency order): "
              + ", ".join(out_of_order))

    rendered = render_graph_module(graph)
    current = TARGET.read_text(encoding="utf-8") if TARGET.exists() else None
    if args.check:
        if current != rendered:
            print(f"{TARGET.name} is out of date. Run scripts/build_fingerprint_graph.py.",
                  file=sys.stderr)
            return 1
        return 0

    if current != rendered:
        TARGET.write_text(rendered, encoding="utf-8", newline="\n")
        print(f"Wrote {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
