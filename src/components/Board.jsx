import React from 'react';

import Help from './Help';
import HelpData from './HelpData';
import Logout from './Logout';
import RunButton from './RunButton';
import Settings from './Settings';

export default class Board extends React.Component {

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
