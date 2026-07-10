import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from charge_backend.charge_server import (
    DEFAULT_DATA_CLASSIFICATION,
    resolve_data_classification,
)


def test_unset_env_returns_default(monkeypatch):
    monkeypatch.delenv("FLASK_DATA_CLASSIFICATION", raising=False)
    assert resolve_data_classification() == DEFAULT_DATA_CLASSIFICATION


def test_valid_json_object_is_parsed(monkeypatch):
    config = {
        "fallbackLevel": "UNKNOWN",
        "rules": [{"backend": "livai", "level": "UNCLASSIFIED data"}],
    }
    monkeypatch.setenv("FLASK_DATA_CLASSIFICATION", json.dumps(config))
    assert resolve_data_classification() == config


def test_invalid_json_falls_back_to_default(monkeypatch):
    monkeypatch.setenv("FLASK_DATA_CLASSIFICATION", "{not valid json")
    assert resolve_data_classification() == DEFAULT_DATA_CLASSIFICATION


def test_non_object_json_falls_back_to_default(monkeypatch):
    # Valid JSON, but a list rather than the expected object.
    monkeypatch.setenv("FLASK_DATA_CLASSIFICATION", "[1, 2, 3]")
    assert resolve_data_classification() == DEFAULT_DATA_CLASSIFICATION


def test_resolved_config_embeds_cleanly_in_app_config(monkeypatch):
    # The result is injected into the HTML template via json.dumps; confirm it
    # serializes to the expected JS object literal shape.
    config = {"fallbackLevel": "UNKNOWN", "rules": []}
    monkeypatch.setenv("FLASK_DATA_CLASSIFICATION", json.dumps(config))
    snippet = f"DATA_CLASSIFICATION: {json.dumps(resolve_data_classification())}"
    assert snippet.startswith("DATA_CLASSIFICATION: {")
    assert json.loads(snippet[len("DATA_CLASSIFICATION: ") :]) == config
