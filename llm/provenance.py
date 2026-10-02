import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

_MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_LOG_PATH = _MODULE_DIR / "logs" / "provenance.jsonl"


class ProvenanceLogger:
    """Appends one JSON line per run_task() call to a configurable log file, for curation
    provenance and future evaluation. Logs both successful and failed runs."""

    def __init__(self, log_path=None):
        if log_path is None:
            log_path = os.environ.get("LLM_LOG_PATH", DEFAULT_LOG_PATH)
        self.log_path = Path(log_path)
        self.prompts_dir = self.log_path.parent / "prompts"

    def log_run(self, task, profile: str, client, result, prompt_text: str, error: Optional[Exception] = None) -> None:
        """
        task: registry.TaskDefinition. profile: the profile name used. client: the LLMClient
        instance used (for provider_name/model). result: the LLMResult on success, or None if
        `error` is set. prompt_text: the full assembled prompt (system_text + user_text), stored
        once as a sha256-named sidecar file - only the hash + relative path are written to the
        JSONL line, keeping it small and letting identical-input runs be diffed by hash.
        """
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        self.prompts_dir.mkdir(parents=True, exist_ok=True)

        prompt_hash = hashlib.sha256(prompt_text.encode("utf-8")).hexdigest()
        prompt_path = self.prompts_dir / f"{prompt_hash}.txt"
        if not prompt_path.exists():
            prompt_path.write_text(prompt_text)

        record = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "task_name": task.name,
            "task_version": task.version,
            "profile": profile,
            "provider": getattr(client, "provider_name", None),
            "model": getattr(client, "model", None),
            "prompt_sha256": prompt_hash,
            "prompt_path": str(prompt_path.relative_to(self.log_path.parent)),
            "status": "error" if error else "ok",
            "raw_output": self._raw_output(result, error),
            "parsed_output": result.parsed.model_dump(mode="json") if (result and result.parsed) else None,
            "usage": {
                "input_tokens": result.input_tokens if result else None,
                "output_tokens": result.output_tokens if result else None,
                "cache_read_input_tokens": result.cache_read_input_tokens if result else None,
                "cache_creation_input_tokens": result.cache_creation_input_tokens if result else None,
            },
            "latency_ms": result.latency_ms if result else None,
            "stop_reason": result.stop_reason if result else None,
            "retries": result.retries if result else None,
            "error": str(error) if error else None,
        }

        with self.log_path.open("a") as f:
            f.write(json.dumps(record) + "\n")

    @staticmethod
    def _raw_output(result, error) -> Optional[str]:
        if result is not None:
            return result.raw_text
        return getattr(error, "raw_output", None) or getattr(error, "raw_text", None)
