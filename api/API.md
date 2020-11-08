# Anti-Malarial Web API

This is the code for the application's backend API server, which provides access to the rest of the backend's functionality.

## Getting Started (Local Development)

To run the code locally, install the following `pip` packages:

```bash
pip install flask
pip install flask-restful
```

Then, you can start the server just by running `api.py` from the root folder (*not* the API folder):

```bash
python -m api.api
```

...which should generate output like this:

```
 * Serving Flask app "api" (lazy loading)
 * Environment: production
   WARNING: This is a development server. Do not use it in a production deployment.
   Use a production WSGI server instead.
 * Debug mode: on
 * Restarting with stat
 * Debugger is active!
 * Debugger PIN: 320-502-271
 * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
```

The last line specifies the base URL that the API is running on.

## User API Guide

The root endpoint for all API requests for version 1.0 of the API is `/api/v1`

### RecrudescenceTest

Given an input `.xlsx` file with [drug testing data in the proper format (TODO)](), return estimates of which patients specified in the file are reinfections or recrudescences.

*Method(s):* `POST`

*URL Endpoint:* `/api/v1/recrudescences`

#### POST

**Input Parameters (URL):**

-   (optional) `iterations` (integer): The number of iterations to run the algorithm for on the data; higher values are more accurate but take longer to complete. Defaults to 10000.

**Input Parameters (Data):**

-   `file` (.xlsx file, FormData): The file containing the testing data to be processed.

**Input Parameters (URL):**

-   (optional) `advanced` (boolean): Boolean indicating whether or not to return advanced stats. Defaults to False.

**Input Parameters (URL):**

-   (optional) `repeats` (list of integers): List of integers representing loci repeats.

**Example Usage:**

In cURL:

```bash
curl -F "file=@example.xlsx" http://localhost:5000/api/v1/recrudescences?iterations=10
```

In Javascript (via `fetch`):

```javascript
const NUM_ITERATIONS = 100
let formData = new FormData();
formData.append('file', fileInputElement.files[0])
fetch(`/api/v1/recrudescences?iterations=${NUM_ITERATIONS}`, {
    method: 'POST',
    body: formData
})
.then(response => {
    // (...)
});
```

**Example Response:**

```json
{
    "runDate": "2020-09-25T18:14:44+00:00",     # ISO formatted date when request was completed
    "totalRunTime": 1806,                       # Computation time in seconds
    "samples": {
        "BQ17-269": {
            "isRecrudescence": false,
            "probability": 0.842,

        }
    }
}
```

**Potential Errors:**
-   `400` - No input file was provided, or an invalid number of iterations was given
-   `413` - The provided file was too large to be processed
-   `415` - The provided file wasn't a `.xlsx` file
