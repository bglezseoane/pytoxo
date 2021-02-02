PyToxo
======

*A Python library for calculating penetrance tables of any bivariate epistasis model.*

PyToxo is a library for calculating penetrance tables of any bivariate epistasis model, writing in Python. It is a renewed version of the original [Toxo](https://github.com/UDC-GAC/toxo) library.

This work is part of the final project of Borja González Seoane's Bachelor Computer Science Engineering studies.


## Running the tests

PyToxo project uses the Python `unittest` — Unit testing framework to all its associated tests routines. PyToxo has unit and integration level tests. All of them are stored into the `test` directory and respect the format `test_*.py`.

To run the tests, the easiest way is using the helper scripts `test/run_all_tests.py`, `test/run_all_tests_unit.py` and `test/run_all_tests_integration.py` to, respectively, run all, the unit ones or the integration ones tests. This script solve by their own the execution working directory, which must be the project home one.

It is also possible to run the tests since the command line. E.g.:

```sh
python -m unittest discover  # Runs all
python -m unittest discover test/unit  # Runs unit tests
python -m unittest discover test/integration  # Runs integration tests
```

The tests use material generated with Toxo to compare the two programs behaviour.
