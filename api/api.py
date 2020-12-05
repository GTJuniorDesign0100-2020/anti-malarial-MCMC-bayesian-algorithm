from collections import OrderedDict
from typing import List
import datetime
import math
import os
import time
import traceback

from flask import Flask, request
from flask_restful import Resource, Api
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
import numpy as np
import pandas as pd
import werkzeug

from api.algorithm_instance import AlgorithmInstance, AlgorithmResults
from api.algorithm_site_instance import LociRepeatError


# =============================================================================
# Initialize Flask
# =============================================================================

#application variable name is necessary for recognition by aws
application = app = Flask(__name__, static_folder='../build', static_url_path='/')

MAX_FILE_SIZE_MB = 50
app.config['MAX_CONTENT_LENGTH'] = MAX_FILE_SIZE_MB*1024**2
app.config['UPLOAD_EXTENSIONS'] = ['.xlsx']      # Valid filetypes

api = Api(app)

limiter = Limiter(app, key_func=get_remote_address)
limiter.init_app(app)

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

    decorators = [
     limiter.limit("250/hour", error_message="Requests per hour limit exceeded"),
     limiter.limit("5/minute", error_message="Requests per minute limit exceeded"),
     limiter.limit("1/second", error_message="Requests per second limit exceeded")
     ]
    def post(self):
        '''
        Takes in a posted xlsx data file, analyzes it to determine if each
        case is a recrudescence/reinfection, and returns the calculated results
        '''
        uploaded_file = request.files.get('file')
        iterations = request.args.get('iterations', default=10000, type=int)
        loci_repeats = request.args.getlist('locirepeat', type=int)
        if not loci_repeats:
            return error_response('Mandatory paramater locirepeats not included')
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
            json_results = self._get_test_results_json(uploaded_file, iterations, loci_repeats)
        except LociRepeatError as e:
            print(e)
            return error_response(str(e), 400)
        except Exception as e:
            # TODO: Return more specific error message?
            print(traceback.format_exc())
            return error_response('A problem occurred while the server was processing this data', 500)

        return json_results, 200

    def _get_test_results_json(
        self,
        uploaded_file: werkzeug.datastructures.FileStorage,
        iterations: int,
        locirepeats: List[int],
        ):
        '''
        Runs a recrudescence test on the given data and returns the test results
        via JSON (plus the status code of the test). If an error occurs, raises
        an exception.
        '''
        run_start_datetime = datetime.datetime.utcnow()
        start_time = time.time()

        test_run = AlgorithmInstance(uploaded_file.stream, locirepeats)

        # calculate burnin (number of runs to discard) and record interval (which n_th iterations should be recorded)
        record_interval = math.ceil(iterations / 1000)
        burnin = math.ceil(iterations * 0.25)

        try:
            results = test_run.run_algorithm(iterations, burnin, record_interval)
        except LociRepeatError as e:
            raise e

        posterior_recrudescence_distribution_df, probability_of_recrudescence_df = results.get_summary_stats()

        end_time = time.time()
        json_response = {
            'runDate': run_start_datetime.isoformat(),
            'totalRunTime': int(end_time - start_time),
            'site_samples': self._get_sample_information(results.site_sample_ids, probability_of_recrudescence_df),
            'output_file_text': results.get_output_file_text()
        }
        return json_response

    def _get_sample_information(self, site_sample_ids: OrderedDict, probability_of_recrudescence_df: pd.DataFrame):
        '''
        Return a dictionary of the probability of each sample, split by site
        :param sample_ids: A dictionary of site names with a list of sample IDs
        for each site, in order of their appearance in the dataset
        :param probability_of_recrudescence_df:
        :return: A dictionary of site names, each of which has a sub-dictionary
        containing the sites samples and
        '''
        samples = {}
        sample_index = 0
        for site_name, ids in site_sample_ids.items():
            site_samples = {}
            for sample_id in ids:
                sample_info = {}
                recrud_prob = probability_of_recrudescence_df.iloc[sample_index]
                sample_info['recrud_probability'] = recrud_prob.values[0]
                site_samples[sample_id] = sample_info
                sample_index += 1

            samples[site_name] = site_samples
        return samples


# =============================================================================
# Assign URL routes/endpoints
# =============================================================================

api.add_resource(RecrudescenceTest, '/api/v1/recrudescences')

@app.route('/')
@app.errorhandler(404)
def index():
    '''
    By default, serve the frontend webpage itself
    '''
    return app.send_static_file('index.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=False, port=os.environ.get('PORT', 5000), threaded=True)
