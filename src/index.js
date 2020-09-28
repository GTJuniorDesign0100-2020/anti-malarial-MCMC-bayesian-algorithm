import React from 'react';
import ReactDOM from 'react-dom';
import Popup from 'reactjs-popup';
import './index.css';

/*
Creates a settings button that toggles the visibility of a popup allowing the user to view
associated email account, and eventually change language preferences and password.

TODO: add functionality for password change (after creating user), make text abstract
to allow for language changes
*/

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

/*
Creates a button that directs to help document.

TODO: evaluate whether help document should be displayed on page, as popup, or
as separate page depending on length & content of document. Make text abstract
to allow for language changes.
*/

class Help extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            value: '?', clicked: false,
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

class RunButton extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            csvFile: '', numIters: ''
        };
        this.handleChangeFile = this.handleChangeFile.bind(this);
        this.handleSubmit = this.handleSubmit.bind(this);
    }

    handleSubmit(event) {
        event.preventDefault();
        return this.state.csvFile;
    }

    handleChangeFile(event) {
        alert('File Submitted for Processing');
        this.setState({csvFile: event.target.value});
    }

  render() {
    return (
      <Popup trigger={<button className="testRun"> Run Test</button>} position="bottom left">
      <div className="testPop">
        <form onsubmit={this.handleSubmit}>
            <label>
                CSV File: <br/>
                <input type="file" name="CSVFile" onChange={this.handleChangeFile}/><br/><br/>
                Convergence Precision +/-: <br/>
                <input type="number" name="convPrecision"/><br/><br/>
                Number of Chr Fragments: <br/>
                <input type="number" name="numFrags"/><br/><br/>
                Number of Iterations: <br/>
                <input type="number" name="numIts"/><br/><br/>
            </label>
            <input type="submit" value="Run Test"/>
        </form>
        </div>
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
      <button className="helpData">
        {this.props.value}
      </button>
    );
  }
}


class Logout extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            value: 'Logout', clicked: false,
        };
    }

  render() {
    return (
      <Popup trigger={<button className="logout"> Log Out</button>} position="bottom right">
        <div className="logoutData">This is filler content for the Logout Button!
        <br/><br/><br/><br/>
        <br/><br/><div align="center" padding="10"></div></div>
      </Popup>
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
      this.state.help[4] = "Logout";
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

  handleDLClick(fileInputElement) {
      const NUM_ITERATIONS = 100
      let formData = new FormData();
      let a = document.createElement('a');
      formData.append('file', fileInputElement.files[0])
      fetch(`/api/v1/recrudescences?iterations=${NUM_ITERATIONS}`, {
          method: 'POST',
          body: formData
      })
      .then(response => {
          a.download = response.json();
          a.click();
      });
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
            onClick={() => this.handleDLClick(i)}
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

  renderLogout(i) {
      return (
          <Logout
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
    const helpText = 'How to use application:'

    return (
      <div>
        <div className="status">{status}</div>
        <div className="helpText">{helpText}</div>
        <br></br>
        <div className="row">
          {this.renderSettings(2)}
          {this.renderLogout(4)}
          {this.renderHelp(0)}
          {this.renderHelpData(1)}
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
