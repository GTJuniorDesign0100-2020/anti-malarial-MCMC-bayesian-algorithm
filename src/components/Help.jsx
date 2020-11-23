import React from 'react';
import Popup from 'reactjs-popup';

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
  }

  render() {
    return (
      <Popup
        trigger={
          <button className="help main-button">
            Help
          </button>
        }
        position="right">
        <div className="helpData">If you need help with how to use the application, check out the <u>FAQ</u>!</div>
      </Popup>
    );
  }
}
