# anti-malarial-MCMC-bayesian-algorithm

Georgia Tech Junior Design 2020 project; built in collaboration with Dr. Mateusz Plucinski and the Georgia Tech 0100 team.

## To Install Locally

1.  Clone the git repo to your computer
2.  In addition to git, ensure you have the following programs installed:

-   [Python 3.8+](https://www.python.org/downloads/)
-   [pip](https://pip.pypa.io/en/stable/installing/)

3.  Install the following Python packages via pip (if you aren't sure what this does, see this [pip installer guide](https://packaging.python.org/tutorials/installing-packages/))

```bash
pip install numpy
pip install pandas
pip install scipy
```

4.  (Optional) If you want to develop for this project, make sure to also install [pytest](https://docs.pytest.org/en/stable/getting-started.html) to run our unit tests and verify the code is working.

```bash
pip install pytest
```

## To Run

1.  Open `main.py` in your text editor of choice (e.g. Notepad or TextEdit); if needed, change the line `inputfile = "Angola2017_example.xlsx"` to the path to your `.xlsx` input file
    -   You can adjust other settings in this file, such as the number of iterations
2.  In your command line, run the main program via `python main.py` in the project's root folder. You should see the program begin running and printing to your console, and generating `.csv` files upon completion.
    -   If the program stops running with an error, double-check that the path to your input file is correct

## To Run Unit Tests

1.  In your command line, run `pytest tests` in the project's root folder to execute all currently existing unit tests. Upon completion, pytest should print if any tests failed.
