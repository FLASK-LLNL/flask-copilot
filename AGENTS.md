# Repository Guidelines

## Project Structure & Module Organization

This repository combines a Python backend with a React frontend. Backend source lives in `charge_backend/`, with FastAPI/WebSocket entry points in `charge_backend/charge_server.py`, shared helpers such as `backend_manager.py`, and focused packages like `pdf/`, `retrosynthesis/`, `lmo/`, and `moleculedb/`. Tests are in `charge_backend/tests/`.

Frontend code lives in `flask-app/src/`, with reusable UI in `components/`, hooks in `hooks/`, and RDKit code in `rdkit/`. Build output goes to `flask-app/dist/`.

## External Submodules

`externals/ChARGe` is the `charge` Python agent framework used by the backend for tasks, clients, and experiments. `externals/lc_conductor` provides LC-Conductor Python utilities plus the `lcc_ui_components` React library consumed by `flask-app/package.json`. Treat both as separate repositories: make intentional commits inside the submodule, then update the parent repo's submodule pointer.

## Reuse-First Development

Before adding code, inspect `ARCHITECTURE.md` and search for existing helpers, models, hooks, components, and submodule APIs. Reuse `charge_backend` domain helpers, LC-Conductor orchestration/tooling, ChARGe task/backend abstractions, and existing React components before creating new primitives. Add abstractions only when existing ones cannot represent the behavior cleanly.

## Local Skills

Shared Codex skills live in `.codex/skills/`. To enable them locally, symlink or copy each skill directory into `$CODEX_HOME/skills` or `~/.codex/skills`; agents should still follow the reuse-first guidance here even when the skill is not installed.

## Build, Test, and Development Commands

- `pip install -e .[all]`: install the backend package and optional dependencies.
- `flask-copilot-install --extras all`: run package setup after installation.
- `python mock_server.py`: run a mock backend for frontend development.
- `cd flask-app && npm install`: install frontend dependencies.
- `cd flask-app && npm run dev`: start Vite; this first builds `externals/lc_conductor/lcc_ui_components`.
- `cd flask-app && npm run build`: create a production frontend build.
- `pytest charge_backend/tests`: run backend tests.
- `pre-commit run --all-files`: run formatting and lint hooks before submitting.

## Coding Style & Naming Conventions

Python uses Black via pre-commit. Keep functions and modules in `snake_case`, classes in `PascalCase`, and tests named `test_*.py`. Prefer small, typed helpers.

TypeScript and CSS use Prettier: semicolons, single quotes, 2-space indentation, trailing commas where valid, and 100-character print width. React components use `PascalCase`; hooks use `useCamelCase`.

## Testing Guidelines

Use `pytest` for backend coverage. Add tests near changed behavior in `charge_backend/tests/`, following names like `test_pdf_registry.py` and `test_websocket_callbacks.py`. For frontend changes, at minimum run `npm run build`; add React Testing Library tests for UI state or data transformations.

## Commit & Pull Request Guidelines

Recent commits use concise, imperative summaries, often with PR numbers, such as `Improve npm build (#199)` or `Add image uploads, fix imports (#204)`. Keep commits focused.

Pull requests should include a summary, test results, linked issues when available, and screenshots or recordings for UI changes. Call out configuration, model asset, MCP tool, or dependency changes.

## Security & Configuration Tips

Do not commit secrets, API keys, local server caches, virtualenvs, build artifacts, or generated logs. Be cautious with files such as `flask_copilot_active_tool_servers.json`, `lightning_logs/`, and local `.env`-style configuration.
