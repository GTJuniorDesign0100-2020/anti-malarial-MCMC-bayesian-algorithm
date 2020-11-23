import React from 'react';

import Help from './Help';
import RunButton from './RunButton';
import Settings from './Settings';
import DynTable from './DynTable';
import PresetLoadingBar from './PresetLoadingBar';
import {recrudescenceAPIRequest, estimateRunTime} from '../utils';

export default class MainScreen extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      tableData: {}
    }

    this.createNewAlgoRun = this.createNewAlgoRun.bind(this);
    this.updateRunValues = this.updateRunValues.bind(this);
  }

  createNewAlgoRun(inputFile, locirepeats, numIters) {
    const runDatetime = new Date();
    const runKey = runDatetime.toISOString();
    const estimatedTime = estimateRunTime(inputFile, numIters);

    // Add new run to state
    this.setState({tableData: {...this.state.tableData,
      [runKey]: {
        date: runDatetime,
        inputFilename: inputFile.name,
        status: <PresetLoadingBar timeToComplete={estimatedTime}/>,
        results: {}
      }
    }});

    // Run algorithm on API
    recrudescenceAPIRequest(inputFile, locirepeats, numIters)
      .then(jsonData => {
        // Update w/ algorithm results if successful
        this.updateRunValues(runKey, {
          results: jsonData,
          // TODO: Use HTML instead of string?
          status: `Completed (in ${jsonData.totalRunTime}s)`
        });
      }, errorJSON => {
        // Update w/ algorithm failure message
        // TODO: Use HTML instead of string?
        this.updateRunValues(runKey, {status: `ERROR: ${errorJSON.message}`});
      });
  }

  updateRunValues(runKey, newRunValues) {
    this.setState({tableData: {
      ...this.state.tableData,
      [runKey]: {
        ...this.state.tableData[runKey],
        ...newRunValues
      }
    }});
  }

  render() {
    const tableStyle = {
        gridColumnStart: 2,
        gridRowStart: 5
    };

    return (
      <div className="main-screen">
        <div className="joined-buttons top-left">
          <RunButton handleSubmit={this.createNewAlgoRun}/>
          <Help />
        </div>
        <Settings />

        <div style={tableStyle}>
          <DynTable data={this.state.tableData}/>
        </div>
      </div>
    );
  }
}
