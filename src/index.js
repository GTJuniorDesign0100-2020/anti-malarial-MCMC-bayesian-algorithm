import React from 'react';
import ReactDOM from 'react-dom';
import Popup from 'reactjs-popup';
import DemoLoadingBar from './components/DemoLoadingBar';
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
      <Popup trigger={<button className="settings" style={{display:'flex', justifyContent:'center'}}>Settings</button>} position="bottom right">
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
            csvFile: '', numIters: 10000
        };

        // bind functions to this class
        this.handleChangeFile = this.handleChangeFile.bind(this);
        this.handleSubmit = this.handleSubmit.bind(this);
        this.handleChangeIters = this.handleChangeIters.bind(this);
    }

    handleSubmit(event) {
        event.preventDefault();
        this.props.handleSubmit(this.state.csvFile, this.state.numIters);
    }

    handleChangeFile(event) {
        console.log(`Got CSV file ${event.target.files[0].name}`);
        this.setState({csvFile: event.target.files[0]});
        alert('Got CSV file');
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
                Convergence Precision +/-: <br/>
                <input type="number" name="convPrecision"/><br/><br/>
                Loci Repeats: <br/>
                <input type="text" name="locirepeats"/><br/><br/>
                Number of Chr Fragments: <br/>
                <input type="number" name="numFrags"/><br/><br/>
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
      <div className="logOut">
      <Popup trigger={<button className="logout" style={{display:'flex', justifyContent:'right'}}>Log Out</button>} position="bottom right">
        <div className="logoutData">This is filler content for the Logout Button!

        <br/><br/><div align="center" padding="10"></div></div>
      </Popup>
      </div>
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
      this.state.help[2] = "Help";
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
            handleSubmit={(csvFile, numIters) => {
              this.handleDLClick(csvFile);
            }}
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
    const helpText = '  How to use application:'

    return (
      <table>
        <tr>
          <td>
            <div className='status'>{status}</div>
          </td>
          <td colspan='3' align='right'>
            {this.renderSettings(0)}
          </td>
        </tr>
        <tr>
          <td colspan='4'>
            {this.renderLogout(3)}
          </td>
        </tr>
        <tr>
          <td padding='10'>
            {helpText}
          </td>
          <td>
            {this.renderHelp(2)}
          </td>
          <td>
            {this.renderHelpData(1)}
          </td>
          <td>
          </td>
        </tr>
        <tr><td><br/></td></tr>
        <tr>
          <td>
            {this.renderTestRun(4)}
          </td>
        </tr>
      </table>
    );
  }
}

class MainScreen extends React.Component {
  render() {
    return (
      <div className="mainscreen">
        <div className="main-board">
          <Board />
          <DemoLoadingBar />
        </div>
      </div>
    );
  }
}

ReactDOM.render(
  <MainScreen />,
  document.getElementById('root')
);
