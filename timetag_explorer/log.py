"""One logger for the package and a helper that writes key=value event lines to it.

The logger is called 'tag_explorer', the same one the web app configures, so the app's log file also receives
what the science code reports. Without any configuration Python prints warnings to stderr and nothing else."""
from __future__ import annotations

import json
import logging

logger = logging.getLogger('tag_explorer')


def _value(v) -> str:
    text = str(v)
    return json.dumps(text) if (not text or any(c.isspace() or c in '"=' for c in text)) else text


def audit(event: str, level: int = logging.INFO, **fields) -> None:
    """One log line: the event name, then key=value pairs (values with spaces are quoted)."""
    logger.log(level, '%s %s', event, ' '.join(f'{k}={_value(v)}' for k, v in fields.items() if v is not None))
