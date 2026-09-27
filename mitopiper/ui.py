"""Small helpers for the interactive (menu driven) interface."""

from __future__ import annotations

from typing import Callable, Optional, Sequence

from tabulate import tabulate

from .metadata import parse_selection

EXIT_WORDS = ("exit", "q", "quit")


def banner(text: str) -> None:
    width = max(len(line) for line in text.split("\n"))
    print(f"{'#' * width}\n{text}\n{'#' * width}\n")


def table(rows, headers="keys") -> str:
    return tabulate(rows, headers=headers, tablefmt="psql")


def ask(prompt: str, default: Optional[str] = None, input_fn: Callable = input) -> str:
    suffix = f" [{default}]" if default not in (None, "") else ""
    answer = input_fn(f"{prompt}{suffix}\n> ").strip()
    return answer if answer else (default or "")


def yes_no(prompt: str, default: bool = False, input_fn: Callable = input) -> bool:
    hint = "Y/n" if default else "y/N"
    answer = input_fn(f"{prompt} ({hint})\n> ").strip().lower()
    if not answer:
        return default
    return answer.startswith("y")


def choose(
    prompt: str,
    options: Sequence[str],
    allow_exit: bool = True,
    input_fn: Callable = input,
    descriptions: Optional[Sequence[str]] = None,
) -> Optional[int]:
    """Show a numbered list and return the chosen index (``None`` on exit)."""
    if not options:
        print("Nothing to choose from.\n")
        return None
    while True:
        print(prompt)
        for i, option in enumerate(options):
            desc = f"\t{descriptions[i]}" if descriptions else ""
            print(f"  {i}\t{option}{desc}")
        if allow_exit:
            print("  (type 'exit' to go back)")
        answer = input_fn("> ").strip()
        if allow_exit and answer.lower() in EXIT_WORDS:
            return None
        if answer.isdigit() and int(answer) < len(options):
            return int(answer)
        print(f"Error: '{answer}' is not one of the listed options.\n")


def choose_many(prompt: str, size: int, input_fn: Callable = input, allow_exit: bool = True) -> Optional[list]:
    """Ask for a selection like ``all``, ``3``, ``1-4`` (inclusive) or ``0 2 5``."""
    while True:
        answer = input_fn(f"{prompt}\n> ").strip()
        if allow_exit and answer.lower() in EXIT_WORDS:
            return None
        try:
            return parse_selection(answer, size)
        except ValueError as exc:
            print(f"Error: {exc}\n")
