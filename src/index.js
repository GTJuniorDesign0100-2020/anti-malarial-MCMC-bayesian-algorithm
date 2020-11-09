import React from 'react';
import ReactDOM from 'react-dom';

import Board from './components/Board';
import DemoLoadingBar from './components/DemoLoadingBar';
import './index.css';

class MainScreen extends React.Component {
  render() {
    return (
      <div className="mainscreen">
          <Board />
          <DemoLoadingBar />
      </div>
    );
  }
}

ReactDOM.render(
  <MainScreen />,
  document.getElementById('root')
);
