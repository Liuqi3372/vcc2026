import ast
import hashlib
import io
import json
from pathlib import Path
import tokenize

root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / 'results/delivery_manifest.json').read_text(encoding='utf-8'))
for name, expected in manifest['files'].items():
    path = root / name
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError('Checksum mismatch: ' + name)
count = 0
for path in root.rglob('*.py'):
    text = path.read_text(encoding='utf-8')
    compile(text, str(path), 'exec')
    tree = ast.parse(text)
    if any(isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant) and isinstance(node.value.value, str) for node in ast.walk(tree)):
        raise ValueError('Unexpected standalone string: ' + str(path))
    if any(token.type == tokenize.COMMENT for token in tokenize.generate_tokens(io.StringIO(text).readline)):
        raise ValueError('Unexpected comment: ' + str(path))
    count += 1
print(json.dumps({'checksummed_files': len(manifest['files']), 'python_files': count, 'syntax_valid': True, 'comments_removed': True, 'full_inference_rerun': False}))
