import React from 'react';

import Help from './Help';
import Logout from './Logout';
import RunButton from './RunButton';
import Settings from './Settings';

export default class Board extends React.Component {
    constructor(props) {
      super(props);
    }

    handleDLClick(csvFile) {
        const NUM_ITERATIONS = 100;
        let formData = new FormData();
        let a = document.createElement('a');
        formData.append('file', csvFile);
        fetch(`/api/v1/recrudescences?iterations=${NUM_ITERATIONS}`, {
            method: 'POST',
            body: formData
        })
        .then(response => {
            response.json().then(jsonData => {
              console.log(jsonData);
              // TODO: Does this link-clicking actually do anything?
              a.download = jsonData;
              a.click();
            });
        });
    }

    render() {
      return (
        <table className="board">
          <tr>
            <td>
              <div className='status'>Welcome!</div>
            </td>
            <td colspan='3' align='right'>
              <Settings />
            </td>
          </tr>
          <tr>
            <td colspan='4'>
              <Logout />
            </td>
          </tr>
          <tr>
            <td padding='10'>
              How to use application:
            </td>
            <td>
              <Help />
            </td>
            <td>
            </td>
            <td>
            </td>
          </tr>
          <tr><td><br/></td></tr>
          <tr>
            <td>
              <RunButton handleSubmit={(csvFile, numIters) => {
                  this.handleDLClick(csvFile);
                }}
              />
            </td>
          </tr>
        </table>
      );
    }
  }
