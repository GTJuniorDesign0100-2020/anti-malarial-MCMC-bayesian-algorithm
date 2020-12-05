import React from 'react';
import Popup from 'reactjs-popup';

/*
Creates a settings button that toggles the visibility of a popup allowing the user to view
associated email account, and eventually change language preferences and password.

TODO: add functionality for password change (after creating user), make text abstract
to allow for language changes
*/

export default function Settings(props) {
  return (
    <Popup trigger={<button className="settings main-button shadow" style={{display:'flex', justifyContent:'center'}}>Settings</button>} position="bottom right">
      <div className="settingsData">Select a language: <br/><select name='language'><option value='english'>English</option></select></div>
    </Popup>
  );
}
