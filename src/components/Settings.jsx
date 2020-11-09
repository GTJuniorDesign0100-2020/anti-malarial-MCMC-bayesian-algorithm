import React from 'react';
import Popup from 'reactjs-popup';

/*
Creates a settings button that toggles the visibility of a popup allowing the user to view
associated email account, and eventually change language preferences and password.

TODO: add functionality for password change (after creating user), make text abstract
to allow for language changes
*/

export default class Settings extends React.Component {
    constructor(props) {
        super(props);
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
