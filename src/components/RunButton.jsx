import React from 'react';
import Popup from 'reactjs-popup';
import {estimateRunTime} from '../utils';

export default class RunButton extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            inputFile: '', locirepeatsString: '', numIters: 10000
        };

        // bind functions to this class
        this.handleSubmit = this.handleSubmit.bind(this);
        this.handleChangeFile = this.handleChangeFile.bind(this);
        this.handleChangeLocirepeats = this.handleChangeLocirepeats.bind(this);
        this.handleChangeIters = this.handleChangeIters.bind(this);
    }

    handleSubmit(event) {
        event.preventDefault();
        // Convert string into a list of locirepeats
        let locirepeats = this.state.locirepeatsString.split(',')
          .map(substring => parseInt(substring.trim(), 10))
          .filter(substring => substring); // Eliminate empty/non-numerical strings

        this.props.handleSubmit(this.state.inputFile, locirepeats, this.state.numIters);
    }

    handleChangeFile(event) {
        console.log(`Got file ${event.target.files[0].name}`);
        this.setState({inputFile: event.target.files[0]});
    }

    handleChangeLocirepeats(event) {
        this.setState({locirepeatsString: event.target.value});
    }

    handleChangeIters(event) {
        this.setState({numIters: event.target.value});
    }

  render() {
    const timeEstimationStyles = {
      marginLeft: '5px',
    }
    return (
      <Popup trigger={<button className="testRun main-button">Run Test</button>} position="bottom left">
        <div className="testPop">
          <form onSubmit={this.handleSubmit}>
            <label htmlFor="InputFile">Excel File:</label><br/>
            <input type="file" name="InputFile" onChange={this.handleChangeFile}/><br/><br/>

            <label htmlFor="locirepeats">Loci Repeats:</label><br/>
            <input type="text" value={this.state.locirepeatsString} name="locirepeats" onChange={this.handleChangeLocirepeats}/><br/><br/>

            <label htmlFor="numIts">Number of Iterations:</label><br/>
            <input type="number" value={this.state.numIters} name="numIts" onChange={this.handleChangeIters}/>

            <hr/>

            <input type="submit" value="Run Test"/>
            {this.state.inputFile &&
              <span style={timeEstimationStyles}>(Estimated Time: <b>{estimateRunTime(this.state.inputFile, this.state.numIters).toFixed(2)}s</b>)</span>
            }
          </form>
        </div>
      </Popup>
    );
  }
}
