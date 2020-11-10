/**
 * Make a request to the recrudescence API, and return the result if successful
 * as a promise
 * @param inputfile The file data object to pass to the API for analysis
 * @param locirepeats The integer list of locirepeats to pass to the API
 * @param numIterations The number of iterations to run the algorithm for
 */
export function recrudescenceAPIRequest(inputFile, locirepeats, numIterations) {
  let formData = new FormData();
  formData.append('file', inputFile);

  return new Promise((resolve, reject) => {
    fetch(`/api/v1/recrudescences?iterations=${numIterations}`, {
      method: 'POST',
      body: formData
    })
    .then(response => {
      response.json().then(jsonData => {
        // Pass JSON API response to a callback function
        if (response.status < 400) {
          resolve(jsonData);
        } else {
          reject(jsonData);
        }
      });
    });
  });
}

/**
 * Returns the estimated time the algorithm will take to complete for the file
 * in seconds
 * @param inputFile The file to be run by the algorithm
 * @param numIterations The number of iterations the algorithm is set to run for
 */
export function estimateRunTime(inputFile, numIterations) {
  const inputFileSize = inputFile.size;
  return 2.0 + 70.0 * (numIterations / 1000.0) * (inputFileSize / 28500.0)
}
