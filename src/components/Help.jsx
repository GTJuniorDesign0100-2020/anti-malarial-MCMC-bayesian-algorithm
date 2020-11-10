import React from 'react';

/*
Creates a button that directs to help document.

TODO: evaluate whether help document should be displayed on page, as popup, or
as separate page depending on length & content of document. Make text abstract
to allow for language changes.
*/

export default class Help extends React.Component {
  constructor(props) {
      super(props);
      this.state = { help: null };

      this.handleHelpClick = this.handleHelpClick.bind(this);
  }

  handleHelpClick() {
    let helpElement = null
    if (!this.state.help) {
      helpElement = <span className="helpData">"This will be a full help document eventually! Enjoy this placeholder!"</span>;
    }
    this.setState({help: helpElement});
  }

  render() {
    const helpContainerStyles = {
      display: 'flex',
      width: '100%'
    };

    return (
      <div className="help" style={helpContainerStyles}>
        <button className="help main-button" onClick={() => this.handleHelpClick()}>
          Help
        </button>
        {this.state.help}
      </div>
    );
  }
}
