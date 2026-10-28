"""The package stays free of web-server code, and its version numbers agree."""
import ast
import re
import subprocess
import sys
from pathlib import Path

import support  # noqa: F401

root = Path(__file__).resolve().parent.parent
FORBIDDEN = {'nicegui', 'starlette', 'fastapi', 'uvicorn', 'httpx', 'authlib', 's3fs'}

# 1. No module imports a web framework (or anything else the server alone needs).
for path in sorted((root / 'timetag_explorer').glob('*.py')):
    for node in ast.walk(ast.parse(path.read_text())):
        names = [a.name for a in node.names] if isinstance(node, ast.Import) else (
            [node.module or ''] if isinstance(node, ast.ImportFrom) and node.level == 0 else [])
        for name in names:
            assert name.split('.')[0] not in FORBIDDEN, f'{path.name} imports {name}'

# 2. Importing the package does not pull in a web framework or pyplot, and nothing runs at import time.
out = subprocess.run([sys.executable, '-c',
                      'import sys, timetag_explorer\n'
                      'bad = [m for m in sys.modules if m.split(".")[0] in %r or m == "matplotlib.pyplot"]\n'
                      'print(",".join(bad))' % sorted(FORBIDDEN)],
                     capture_output=True, text=True, cwd=root)
assert out.returncode == 0, out.stderr
assert out.stdout.strip() == '', out.stdout

# 3. One version number, everywhere.
import timetag_explorer  # noqa: E402

version = timetag_explorer.__version__
assert re.fullmatch(r'\d+\.\d+\.\d+', version), version
cff = (root / 'CITATION.cff').read_text()
assert re.search(rf'^version: "?{re.escape(version)}"?\s*$', cff, re.M), 'CITATION.cff version differs from __version__'
assert f'## [{version}]' in (root / 'CHANGELOG.md').read_text(), 'CHANGELOG.md has no entry for this version'

# 4. No private data in the package: no personal paths or addresses.
for path in sorted((root / 'timetag_explorer').glob('*.py')):
    text = path.read_text()
    assert '/home/' not in text and '@gmail.com' not in text, path.name
print('package ok')
