import datetime
import math
import os
import time

from flask import Flask, request
from flask_restful import Resource, Api
import numpy as np
import pandas as pd
import werkzeug

from api.algorithm_instance import AlgorithmInstance


# =============================================================================
# Initialize Flask
# =============================================================================

app = Flask(__name__)

MAX_FILE_SIZE_MB = 50
app.config['MAX_CONTENT_LENGTH'] = MAX_FILE_SIZE_MB*1024**2
app.config['UPLOAD_EXTENSIONS'] = ['.xlsx']      # Valid filetypes

api = Api(app)

# =============================================================================
# Define actual functionality
# =============================================================================


def error_response(message, error_code=400):
    '''
    Generates an appropriate JSON error message and HTML error code that can be
    returned from a Flask API handler
    '''
    return {'message': message}, error_code


@app.errorhandler(413)
def too_large(e):
    return error_response(f'Input file is too large (over {MAX_FILE_SIZE_MB} MB)', 413)


class RecrudescenceTest(Resource):
    def __init__(self):
        # TODO: Remove!
        self.PLACEHOLDER_RESPONSE = {
            'runDate': '2020-09-25T18:14:44+00:00',
            'totalRunTime': 1806,
            'samples': {
                'BQ17-269': {
                    'isRecrudescence': False,
                    'probability': 0.842,

                }
            }
        }

    def post(self):
        '''
        Takes in a posted CSV data file, analyzes it to determine if each
        case is a recrudescence/reinfection, and returns the calculated results
        '''
        uploaded_file = request.files.get('file')
        iterations = request.args.get('iterations', default=10000, type=int)
        # TODO: Find cleaner way of doing validation?
        if not (uploaded_file and uploaded_file.filename):
            return error_response('No input file provided')

        file_extension = os.path.splitext(uploaded_file.filename)[-1].lower()
        if file_extension not in app.config['UPLOAD_EXTENSIONS']:
            return error_response('Provided file is not a .xlsx', 415)

        if iterations < 0:
            return error_response(f'{iterations} is an invalid iteration count')

        json_results = {}
        try:
            json_results = self._get_test_results_json(uploaded_file, iterations)
        except Exception as e:
            # TODO: Return more specific error message?
            return error_response('A problem occurred while the server was processing this data', 500)

        return json_results, 200

    def _get_test_results_json(
        self,
        uploaded_file: werkzeug.datastructures.FileStorage,
        iterations: int,
        locirepeats=[2,2,3,3,3,3,3]):
        '''
        Runs a recrudescence test on the given data and returns the test results
        via JSON (plus the status code of the test). If an error occurs,
        returns an error message and a 500 status code.
        '''
        run_start_datetime = datetime.datetime.utcnow()
        start_time = time.time()

        test_run = AlgorithmInstance(uploaded_file.stream, locirepeats)

        # calculate burnin (number of runs to discard) and record interval (which n_th iterations should be recorded)
        record_interval = math.ceil(iterations / 1000)
        burnin = math.ceil(iterations * 0.25)

        posterior_recrudescence_distribution_df, probability_of_recrudescence_df, run_posterior_dfs, run_summary_stat_dfs, sample_ids = test_run.run_algorithm(iterations, burnin, record_interval)

        end_time = time.time()
        json_response = {
            'runDate': run_start_datetime.isoformat(),
            'totalRunTime': int(end_time - start_time),
            'samples': self._get_sample_information(sample_ids, probability_of_recrudescence_df)
        }
        return json_response

    def _get_sample_information(self, sample_ids, probability_of_recrudescence_df: pd.DataFrame):
        '''
        Return a dictionary of the probability of each sample
        TODO: Modify format to better reflect sample information?
        '''
        samples = {}
        for i, sample_id in enumerate(sample_ids):
            recrud_prob = probability_of_recrudescence_df.iloc[i]

            sample_info = {}
            sample_info['isRecrudescence'] = False if recrud_prob <= 0.5 else True
            sample_info['probability'] = recrud_prob

            samples[sample_id] = sample_info
        return samples



# =============================================================================
# Assign URL routes/endpoints
# =============================================================================

api.add_resource(RecrudescenceTest, '/api/v1/recrudescences')

if __name__ == '__main__':
    app.run(port=5000, debug=True)
