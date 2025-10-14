
# Contributing

- Discuss major changes in an issue before opening a PR.
- Add tests for new equations and keep units explicit in docstrings.
- Update docs in `documentation/` when public APIs change.

## Dev quickstart
```bash
python -m venv .venv && source .venv/bin/activate
pip install -e .[yaml]
pip install -U pytest black ruff pre-commit
pre-commit install
pytest -q
```
