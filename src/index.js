import React from 'react';
import ReactDOM from 'react-dom';
import Popup from 'reactjs-popup';
import './index.css';


class Help extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            value: 'Help', clicked: false,
        };
    }

  render() {
    return (
      <button className="help" onClick={() => this.props.onClick()}>
        {this.props.value}
      </button>
    );
  }
}

class Settings extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            value: 'Settings', clicked: false,
        };
    }

  render() {
    return (
      <Popup trigger={<button className="settings"> Settings</button>} position="bottom right">
        <div className="settingsData">This is filler content for the settings tab!
        <br/><br/>Username: FakeName@FakeISP.FakeDomain<br/><br/>Password: ********
        <br/><br/><div align="center" padding="10"><button>Save</button></div></div>
      </Popup>
    );
  }
}

class RunButton extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            value: 'Settings', clicked: false,
        };
    }

  render() {
    return (
      <Popup trigger={<button className="testRun"> Run Test</button>} position="bottom left">
        <div className="testPop">Please upload data to test!<br/><br/>
        <input type="file"></input><br/><br/>Convergence Precision: +/-
        <input type="number"></input>%<br/><br/>Number of Chr Fragments:
        <input type="number"></input><br/><br/><div align="center"><button>Run Test</button><br/>
        </div></div>
      </Popup>
    );
  }
}

class HelpData extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            value: 'Help', clicked: false,
        };
    }

  render() {
    return (
      <button className="helpData" onClick={() => this.props.onClick()}>
        {this.props.value}
      </button>
    );
  }
}

class Board extends React.Component {

  constructor(props) {
      super(props);
      this.state = {
          help: Array(5).fill(null),
      };
      this.state.help[0] = "Help";
      this.state.help[2] = "Settings";
      this.state.help[3] = "Run Test";
  }

  handleHelpClick(i) {
      const help = this.state.help.slice();
      if (help[1] == null) {
          help[1] = <div className="helpdata">"This will be a full help document eventually! Enjoy this placeholder!"</div>;
      } else {
          help[1] = null;
      }
      this.setState({help: help});

  }

  renderHelp(i) {
      return (
          <Help
            value={this.state.help[i]}
            onClick={() => this.handleHelpClick(i)}
            />
        );
  }

  renderTestRun(i) {
      return (
          <RunButton
            value={this.state.help[i]}
          />
      );
  }

  renderSettings(i) {
      return (
          <Settings
            value={this.state.help[i]}
            />
        );
  }

  renderHelpData(i) {
      return (
          <HelpData
            value={this.state.help[i]}
            />
        );
  }


  render() {
    const status = 'Welcome!';

    return (
      <div>
        <div className="status">{status}</div>
        <div className="row">
          {this.renderHelp(0)}
          {this.renderHelpData(1)}
          {this.renderSettings(2)}
        </div>
        <div className="row">
          {this.renderTestRun(3)}
        </div>
      </div>
    );
  }
}

class MainScreen extends React.Component {
  render() {
    return (
      <div className="mainscreen">
        <div className="main-board">
          <Board />
        </div>
        <div className="main-info">
          <div></div>
          <ol></ol>
        </div>
      </div>
    );
  }
}

ReactDOM.render(
  <MainScreen />,
  document.getElementById('root')
);
