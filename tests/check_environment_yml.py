#!/usr/bin/env python3
"""Assert environment.yml lists exactly the packages the modules ask for.

environment.yml duplicates information that lives in the per-process `conda` directives, so it
will drift the moment someone adds a tool to a module and forgets the aggregate file. Compare the
two package sets and fail on any difference.

Versions are deliberately *not* compared: three packages are requested at different versions by
different modules, and environment.yml collapses those to the newest, which is documented in its
header. Only the set of packages has to agree.

Deliberately parses environment.yml by hand rather than importing yaml, so the check has no
dependency beyond the standard library.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def package_name(spec):
    """'bioconda::samtools=1.24' -> 'samtools'; also strips <, >, <=, >= bounds."""
    spec = spec.split("::", 1)[1] if "::" in spec else spec
    return re.split(r"[=<>!]", spec, maxsplit=1)[0].strip()


def packages_from_modules():
    found = {}
    for nf in sorted((ROOT / "modules" / "local").glob("*.nf")):
        for directive in re.findall(r'conda\s+"([^"]+)"', nf.read_text()):
            for token in directive.split():
                if "::" in token:
                    found.setdefault(package_name(token), set()).add(nf.name)
    return found


def packages_from_environment():
    found = set()
    in_deps = False
    for line in (ROOT / "environment.yml").read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if re.match(r"^dependencies:\s*$", stripped):
            in_deps = True
            continue
        # any other top-level key ends the dependencies block
        if in_deps and not line.startswith((" ", "\t", "-")):
            in_deps = False
        if in_deps and stripped.startswith("- "):
            found.add(package_name(stripped[2:]))
    return found


def main():
    modules = packages_from_modules()
    env = packages_from_environment()

    missing = sorted(set(modules) - env)
    extra = sorted(env - set(modules))

    for name in missing:
        print(
            f"  MISSING from environment.yml: {name} "
            f"(required by {', '.join(sorted(modules[name]))})",
            file=sys.stderr,
        )
    for name in extra:
        print(f"  in environment.yml but no module asks for it: {name}", file=sys.stderr)

    if missing or extra:
        print(
            "\nenvironment.yml is out of step with the module conda directives.",
            file=sys.stderr,
        )
        return 1

    print(f"OK: environment.yml lists the same {len(env)} packages as the module specs")
    return 0


if __name__ == "__main__":
    sys.exit(main())
