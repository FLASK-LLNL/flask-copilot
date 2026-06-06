# Architecture Specification

## Purpose

FLASK Copilot is a molecular discovery agentic framework realized as a web application. It combines a React/Vite frontend with a Python FastAPI backend that runs LLM-backed chemistry workflows through the ChARGe agent framework and LC-Conductor orchestration utilities.

This document is for agents developing or debugging code. Use it to identify the owning module before changing behavior.

## Repository Boundaries

- `charge_backend/`: FLASK Copilot-specific backend behavior, domain workflows, WebSocket actions, PDF support, and chemistry tools.
- `flask-app/`: React frontend for graph-based molecule workflows, project state, settings, attachments, and result visualization.
- `externals/charge`: submodule providing the generic `ChARGe` task, client, backend, and experiment framework.
- `externals/lc_conductor`: submodule providing shared backend orchestration utilities and reusable React components consumed by the frontend.

Changes to generic agent execution usually belong in `externals/ChARGe`. Changes to tool registration, selected MCP runtime, local MCP proxying, or shared UI components usually belong in `externals/lc_conductor`. Changes to FLASK-specific chemistry workflows belong in `charge_backend/` or `flask-app/src/`.

## Runtime Flow

The backend entry point is `charge_backend/charge_server.py`. It creates a FastAPI app, registers REST endpoints for MCP server management, and exposes the main `/ws` WebSocket endpoint. Each WebSocket connection creates connection-scoped state: a `TaskManager`, an `Experiment`, a `FlaskActionManager`, and optional registries such as `PdfDocumentRegistry`.

Frontend actions are sent over the WebSocket as typed JSON messages. `FlaskActionManager` in `charge_backend/backend_manager.py` routes compute actions by `problemType`:

- `optimization`: runs lead molecule optimization through `charge_backend/lmo/`.
- `retrosynthesis`: runs template, database, and AI retrosynthesis paths through `charge_backend/retrosynthesis/`.
- `custom`: runs custom user prompts through `charge_backend/charge_backend_custom.py`.

The action manager delegates generic task execution, tool runtime selection, save/load state, and callback messaging to LC-Conductor's base `ActionManager`.

## Backend Domain Modules

`charge_backend/backend_helper_funcs.py` defines shared data shapes used in WebSocket payloads, including molecule graph nodes, reactions, alternatives, callbacks, and run settings.

`charge_backend/retrosynthesis/` owns reaction planning. `template.py` wraps AiZynthFinder/template planning, `ai.py` builds ChARGe tasks for AI retrosynthesis, `database.py` queries known reaction data, `alternatives.py` mutates selected reaction alternatives, and `context.py` tracks retrosynthesis state across a run.

`charge_backend/lmo/` owns lead molecule optimization tasks. `lmo_task.py` defines structured output schemas and task prompts; `lmo_charge_backend_funcs.py` turns agent output into graph nodes.

`charge_backend/pdf/` owns uploaded PDF reference handling. The registry stores per-user active documents; `pdf_document.py` implements page reading, search, and formatting for the `consult_with_document` built-in tool.

## Frontend Architecture

`flask-app/src/App.tsx` is the main coordinator. It owns WebSocket lifecycle, selected problem type, run settings, tool selection, graph state, PDF/image attachments, metrics, reasoning sidebar state, and modal state.

Important supporting modules:

- `components/graph.tsx`: molecule graph rendering and graph state helpers.
- `tree_utils.ts`: tree layout, descendant traversal, and graph mutation helpers.
- `hooks/useProjectData.ts`: project persistence abstraction backed by local storage.
- `components/project_sidebar.tsx`: project list and saved context workflows.
- `components/reaction_alternatives.tsx`: retrosynthesis alternative selection UI.
- `components/metrics.tsx`: run and token metrics display.
- `rdkit/reactionPayload.ts`: RDKit-related reaction payload construction.

Shared UI from `lc-conductor` is imported for settings, attachment upload, markdown rendering, reasoning sidebar, and local MCP proxy handling.

## Message and State Model

The frontend sends action messages such as compute, cancel, save/load context, MCP server updates, and PDF reference configuration. The backend sends response messages for reasoning text, molecule graph updates, metrics, save/load responses, PDF reference status, and completion.

Connection-local backend state should stay inside `TaskManager`, `FlaskActionManager`, `Experiment`, or registries. Frontend session/project state should go through `useProjectData` or React state in `App.tsx`.

## Debugging Guide

- WebSocket routing or missing completion: inspect `charge_server.py`, `FlaskActionManager.handle_compute`, and LC-Conductor `TaskManager._handle_task_done`.
- Wrong tool availability: inspect `selected_tool_runtime`, `builtin_tools.py`, and LC-Conductor `tool_registration.py` / `local_mcp_proxy.py`.
- Bad LLM task behavior: inspect FLASK task construction first, then ChARGe backend implementation.
- Graph mismatch: inspect backend node payload construction and frontend `graph.tsx` / `tree_utils.ts`.
- PDF reference issues: inspect `PdfDocumentRegistry`, `SubagentDocumentPdf`, and frontend attachment state.
