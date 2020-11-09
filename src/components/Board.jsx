import React from 'react';

import Help from './Help';
import Logout from './Logout';
import RunButton from './RunButton';
import Settings from './Settings';
import {recrudescenceAPIRequest} from '../utils';

export default function Board(props) {
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
          <RunButton handleSubmit={(inputFile, locirepeats, numIters) => {
              recrudescenceAPIRequest(inputFile, locirepeats, numIters)
                .then(jsonData => console.log(jsonData));
            }}
          />
        </td>
      </tr>
    </table>
  );
}
