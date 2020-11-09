import React from 'react';
import Popup from 'reactjs-popup';

export default class Logout extends React.Component {
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
