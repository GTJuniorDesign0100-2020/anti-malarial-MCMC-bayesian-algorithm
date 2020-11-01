This project was bootstrapped with [Create React App](https://github.com/facebook/create-react-app).

## Available Scripts

In the project directory, you can run:

### `npm start`

Runs the app in the development mode.<br />
Open [http://localhost:3000](http://localhost:3000) to view it in the browser.

The page will reload if you make edits.<br />
You will also see any lint errors in the console.

### `npm test`

Launches the test runner in the interactive watch mode.<br />
See the section about [running tests](https://facebook.github.io/create-react-app/docs/running-tests) for more information.

### `npm run build`

Builds the app for production to the `build` folder.<br />
It correctly bundles React in production mode and optimizes the build for the best performance.

The build is minified and the filenames include the hashes.<br />
Your app is ready to be deployed!

See the section about [deployment](https://facebook.github.io/create-react-app/docs/deployment) for more information.

### `npm run eject`

**Note: this is a one-way operation. Once you `eject`, you can’t go back!**

If you aren’t satisfied with the build tool and configuration choices, you can `eject` at any time. This command will remove the single build dependency from your project.

Instead, it will copy all the configuration files and the transitive dependencies (webpack, Babel, ESLint, etc) right into your project so you have full control over them. All of the commands except `eject` will still work, but they will point to the copied scripts so you can tweak them. At this point you’re on your own.

You don’t have to ever use `eject`. The curated feature set is suitable for small and middle deployments, and you shouldn’t feel obligated to use this feature. However we understand that this tool wouldn’t be useful if you couldn’t customize it when you are ready for it.

## Learn More

You can learn more in the [Create React App documentation](https://facebook.github.io/create-react-app/docs/getting-started).

To learn React, check out the [React documentation](https://reactjs.org/).

### Code Splitting

This section has moved here: https://facebook.github.io/create-react-app/docs/code-splitting

### Analyzing the Bundle Size

This section has moved here: https://facebook.github.io/create-react-app/docs/analyzing-the-bundle-size

### Making a Progressive Web App

This section has moved here: https://facebook.github.io/create-react-app/docs/making-a-progressive-web-app

### Advanced Configuration

This section has moved here: https://facebook.github.io/create-react-app/docs/advanced-configuration

### Deployment

This section has moved here: https://facebook.github.io/create-react-app/docs/deployment

### `npm run build` fails to minify

This section has moved here: https://facebook.github.io/create-react-app/docs/troubleshooting#npm-run-build-fails-to-minify
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
