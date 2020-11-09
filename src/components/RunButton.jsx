import React from 'react';
import Popup from 'reactjs-popup';

export default class RunButton extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            csvFile: '', locirepeatsString: '', numIters: 10000
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
        console.log(locirepeats);

        this.props.handleSubmit(this.state.csvFile, this.state.numIters);
    }

    handleChangeFile(event) {
        console.log(`Got file ${event.target.files[0].name}`);
        this.setState({csvFile: event.target.files[0]});
    }

    handleChangeLocirepeats(event) {
        this.setState({locirepeatsString: event.target.value});
    }

    handleChangeIters(event) {
        this.setState({numIters: event.target.value});
    }

  render() {
    return (
      <Popup trigger={<button className="testRun"> Run Test</button>} position="bottom left">
      <div className="testPop">
        <form onSubmit={this.handleSubmit}>
            <label>
                CSV File: <br/>
                <input type="file" name="CSVFile" onChange={this.handleChangeFile}/><br/><br/>
                Loci Repeats: <br/>
                <input type="text" value={this.state.locirepeatsString} name="locirepeats" onChange={this.handleChangeLocirepeats}/><br/><br/>
                Number of Iterations: <br/>
                <input type="number" value={this.state.numIters} name="numIts" onChange={this.handleChangeIters}/><br/><br/>
            </label>
            <input type="submit" value="Run Test"/>
        </form>
        </div>
      </Popup>
    );
  }
}
