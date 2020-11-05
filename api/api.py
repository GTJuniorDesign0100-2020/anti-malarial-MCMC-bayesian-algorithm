from collections import OrderedDict
import datetime
import math
import os
import time

from flask import Flask, request
from flask_restful import Resource, Api
import numpy as np
import pandas as pd
import werkzeug

from api.algorithm_instance import AlgorithmInstance, AlgorithmResults


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
            return error_response(f'A problem occurred while the server was processing this data: {e}', 500)

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

        results = test_run.run_algorithm(iterations, burnin, record_interval)

        posterior_recrudescence_distribution_df, probability_of_recrudescence_df = results.get_summary_stats()

        end_time = time.time()
        json_response = {
            'runDate': run_start_datetime.isoformat(),
            'totalRunTime': int(end_time - start_time),
            'site_samples': self._get_sample_information(results.site_sample_ids, probability_of_recrudescence_df),
            'output_file_text': self._get_output_file_text(results)
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
                sample_info['recrud_probability'] = recrud_prob
                site_samples[sample_id] = sample_info
                sample_index += 1

            samples[site_name] = site_samples
        return samples

    def _get_output_file_text(self, results: AlgorithmResults):
        '''
        Gets the text of each output .csv file and stores it in a dictionary,
        with the filename as the key
        TODO: Near-duplicated with main.py?
        '''
        output_files = {}

        posterior_recrudescence_distribution_df, probability_of_recrudescence_df = results.get_summary_stats()
        for site_name, posterior_df in results.run_posterior_dfs.items():
            filename = f'{site_name}_posterior.csv'
            output_files[filename] = posterior_df.to_csv()
        for site_name, summary_df in results.run_summary_stat_dfs.items():
            filename = f'{site_name}_summary_statistics.csv'
            output_files[filename] = summary_df.to_csv()

        output_files['microsatellite_correction.csv'] = posterior_recrudescence_distribution_df.to_csv(index=False)
        output_files['probability_of_recrudescence.csv'] = probability_of_recrudescence_df.to_csv(index=False)

        return output_files


# =============================================================================
# Assign URL routes/endpoints
# =============================================================================

api.add_resource(RecrudescenceTest, '/api/v1/recrudescences')

if __name__ == '__main__':
    app.run(port=5000, debug=True)
