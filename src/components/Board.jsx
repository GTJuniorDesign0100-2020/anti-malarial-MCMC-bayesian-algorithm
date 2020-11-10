import React from 'react';

import Help from './Help';
import Logout from './Logout';
import RunButton from './RunButton';
import Settings from './Settings';
import DynTable from './DynTable';
import {recrudescenceAPIRequest} from '../utils';

export default function Board(props) {

  const welcomeStyle = {
    gridColumnStart: 1,
    gridRowStart: 1
  };
  const helpTextStyle = {
    gridColumnStart: 1,
    gridRowStart: 3
  };
  const helpStyle = {
    gridColumnStart: 2,
    gridRowStart: 3
  };
  const tableStyle = {
      gridColumnStart: 2,
      gridRowStart: 5
  };

  return (
    <div className="board">
      <div className='status' style={welcomeStyle}>Welcome!</div>
      <Settings />
      <div style={helpTextStyle}>How to use application:</div>
      <Help style={helpStyle} />
      <RunButton handleSubmit={(inputFile, locirepeats, numIters) => {
              recrudescenceAPIRequest(inputFile, locirepeats, numIters)
                .then(jsonData => console.log(jsonData));

            }}
          />
      <div style={tableStyle}>
        <DynTable />
      </div>
    </div>
  );
}
