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
          <button className="help main-button shadow">
            Help
          </button>
        }
        position="right">
        <div className="help-data">If you need help with how to use the application, check out the <a target="_blank" rel="noopener noreferrer" href="https://docs.google.com/document/d/14xnfxBzDkTYQqryIv3YDtga34cPTaUkFUND-bL-N7UY/edit?usp=sharing">FAQ Page</a>!</div>
      </Popup>
    );
  }
}
