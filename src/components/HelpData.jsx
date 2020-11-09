import React from 'react';

/*
Creates a button that directs to help document.

TODO: evaluate whether help document should be displayed on page, as popup, or
as separate page depending on length & content of document. Make text abstract
to allow for language changes.
*/

export default class HelpData extends React.Component {
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
