from pathlib import Path
from typing import Optional

from . import prompt as prompt_module
from . import provenance as provenance_module
from . import registry
from .providers import profiles as profiles_module
from .providers.base import LLMClient, LLMOutputValidationError, LLMRefusalError

YELLOW = "\033[33m"
RESET = "\033[0m"


def run_task(
    task_name: str,
    inputs: dict,
    profile: str,
    *,
    content_root=None,
    profiles_path=None,
    log_path=None,
    client: Optional[LLMClient] = None,
    **complete_opts,
):
    """
    The public entry point: from llm import run_task; result = run_task("review_evidence_item",
    inputs={...}, profile="claude").

    Loads the task, warns (doesn't error) if the profile doesn't advertise support for what the
    task requires, assembles the prompt, calls the client, logs provenance for both success and
    failure, and returns the LLMResult - or lets LLMOutputValidationError/LLMRefusalError
    propagate after logging the failed attempt.

    `client` is a test seam: pass an already-constructed LLMClient (e.g. testing.FakeClient) to
    bypass get_client()/profile-based construction entirely. `profile` is still used to label the
    provenance record and, where the named profile happens to exist, to run the
    requires-vs-supports warning check.
    """
    content_root = Path(content_root) if content_root else registry.DEFAULT_CONTENT_ROOT
    task = registry.load_task(task_name, content_root)

    try:
        profile_obj = profiles_module.get_profile(profile, profiles_path)
        _warn_on_unsupported_requirements(task, profile_obj)
    except (FileNotFoundError, KeyError):
        pass  # unknown/test-only profile name - nothing to warn about

    if client is None:
        client = profiles_module.get_client(profile, profiles_path)

    system_text, user_text = prompt_module.assemble_prompt(task, inputs)
    prompt_text = f"{system_text}\n\n{user_text}"
    messages = [{"role": "user", "content": user_text}]

    logger = provenance_module.ProvenanceLogger(log_path=log_path)

    try:
        result = client.complete(system=system_text, messages=messages, output_schema=task.schema, **complete_opts)
    except (LLMOutputValidationError, LLMRefusalError) as error:
        logger.log_run(task, profile, client, None, prompt_text, error=error)
        raise

    logger.log_run(task, profile, client, result, prompt_text)

    return result


def _warn_on_unsupported_requirements(task: registry.TaskDefinition, profile_obj) -> None:
    missing = [r for r in task.requires if r not in profile_obj.supports]
    if missing:
        print(
            f"{YELLOW}WARNING: profile '{profile_obj.name}' does not advertise support for "
            f"{missing}, required by task '{task.name}'. Proceeding anyway - results may be "
            f"degraded.{RESET}"
        )
