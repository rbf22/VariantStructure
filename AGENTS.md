This project uses a suite of tools to ensure code quality. Before submitting any changes, please run the following checks and address any issues they report.

## Installation

To run the checks, you first need to install the development dependencies:

```bash
pip install -e .[dev]
```

## Ruff

Ruff is a fast Python linter, written in Rust. It can be used to check for a variety of issues and also to format the code.

To run ruff and automatically fix any issues it finds, run:

```bash
ruff check --fix .
```

To format the code with ruff, run:

```bash
ruff format .
```

## Pylint

Pylint is a static code analysis tool which looks for programming errors, helps enforcing a coding standard and sniffs for some code smells.

To run pylint, use the following command:

```bash
pylint protein_rebuilder
```

## MyPy

MyPy is a static type checker for Python. It helps to catch type errors before they cause issues at runtime.

To run mypy, use the following command:

```bash
mypy protein_rebuilder
```

## Deptry

Deptry is a command line tool to check for obsolete, missing, transitive and unused dependencies in a Python project.

To run deptry, use the following command:

```bash
deptry .
```

## Running all checks

You can run all the checks in one go with the following command:

```bash
ruff check . && ruff format . && pylint protein_rebuilder && mypy protein_rebuilder && deptry .
```
