import React from 'react';

import Help from './Help';
import Logout from './Logout';
import RunButton from './RunButton';
import Settings from './Settings';
import DynTable from './DynTable';
import DemoLoadingBar from './DemoLoadingBar';
import {recrudescenceAPIRequest} from '../utils';

export default class Board extends React.Component {

    constructor(props) {
        super(props);
        this.state = {
            tableData: []
        }
    }


    render() {
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
              this.state.tableData = ['11/9/2020', '21:45', inputFile, <DemoLoadingBar/>];
              recrudescenceAPIRequest(inputFile, locirepeats, numIters)
                .then(jsonData => console.log(jsonData));

            }}
          />
      <div style={tableStyle}>
        <DynTable data={this.state.tableData}/>
      </div>
    </div>
  );
}
}
