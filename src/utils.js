/**
 * Make a request to the recrudescence API, and return the result if successful
 * as a promise
 * @inputfile The file data object to pass to the API for analysis
 * @locirepeats The integer list of locirepeats to pass to the API
 * @numIterations The number of iterations to run the algorithm for
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
        resolve(jsonData);
      });
    }, error => {
      console.error(error);
    });
  });
}
