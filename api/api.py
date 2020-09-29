import os

from flask import Flask, request
from flask_restful import Resource, Api


# =============================================================================
# Initialize Flask
# =============================================================================

app = Flask(__name__)

MAX_FILE_SIZE_MB = 50
app.config['MAX_CONTENT_LENGTH'] = MAX_FILE_SIZE_MB*1024**2
app.config['UPLOAD_EXTENSIONS'] = ['.csv']      # Valid filetypes

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
            return error_response('Provided file is not a .csv', 415)

        if iterations < 0:
            return error_response(f'{iterations} is an invalid iteration count')

        # TODO: Actually process the file
        return self.PLACEHOLDER_RESPONSE, 200


# =============================================================================
# Assign URL routes/endpoints
# =============================================================================

api.add_resource(RecrudescenceTest, '/api/v1/recrudescences')

if __name__ == '__main__':
    app.run(port=5000, debug=True)
